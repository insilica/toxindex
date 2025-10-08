import json
import os
import uuid
import logging
import tempfile
import requests
import openai
from workflows.celery_app import celery
from webserver.model.message import MessageSchema
from webserver.model.task import Task
from pathlib import Path
from workflows.utils import (
    emit_status,
    upload_local_file_to_gcs,
    emit_task_message,
    emit_task_file,
)
logging.getLogger().info("chemical_fetcher_task.py module loaded")
logger = logging.getLogger(__name__)

def extract_limit_with_ai(user_query):
    """Extract the number of chemicals to fetch from user query using AI"""
    try:
        # Use AI to extract the limit number (new OpenAI API format)
        client = openai.OpenAI(api_key=os.environ.get('OPENAI_API_KEY'))
        response = client.chat.completions.create(
            model="gpt-3.5-turbo",
            messages=[
                {
                    "role": "system", 
                    "content": "You are a helpful assistant that extracts numbers from user queries about fetching chemical data. Return ONLY the number as an integer. If no number is found, return 200. If the number is greater than 1000, return 1000. Examples: 'fetch 200 chemicals' -> 200, 'get 50 sample chemicals' -> 50, 'show me some chemicals' -> 200"
                },
                {"role": "user", "content": user_query},
            ],
            max_tokens=10,
            temperature=0.1,
        )
        
        # Extract the number from AI response (new API format)
        ai_response = response.choices[0].message.content.strip()
        
        # Try to extract number from AI response
        import re
        numbers = re.findall(r'\d+', ai_response)
        if numbers:
            limit = min(int(numbers[0]), 1000)  # Cap at 1000 for safety
            logger.info(f"AI extracted limit: {limit} from query: '{user_query}'")
            return limit
        else:
            logger.warning(f"AI response '{ai_response}' contained no numbers, using default 200")
            return 200
            
    except Exception as e:
        logger.warning(f"AI limit extraction failed: {e}, using default 200")
        return 200

def fetch_chemicals_sparql(limit=200):
    """Fetch chemical data from SPARQL endpoint"""
    endpoint = "https://neptune.biobricks.ai/sparql"
    
    sparql_query = f"""
    PREFIX rdfs: <http://www.w3.org/2000/01/rdf-schema#>
    PREFIX rdf: <http://www.w3.org/1999/02/22-rdf-syntax-ns#>
    PREFIX CHEMINF: <http://purl.obolibrary.org/obo/CHEMINF_>
    PREFIX CAS: <http://identifiers.org/cas/>
    PREFIX EDAM: <http://edamontology.org/>

    SELECT ?cas ?dsstox ?name
    WHERE {{
      ?s a CHEMINF:000000 ;
         rdfs:label ?name .
      OPTIONAL {{
        ?s EDAM:has_identifier ?dsstox .
        ?dsstox a CHEMINF:000568 . # DSSTOX SID
      }}
      OPTIONAL {{
        ?s EDAM:has_identifier ?cas .
        ?cas a CHEMINF:000446 . # CAS RN
      }}
    }}
    LIMIT {limit}
    """
    
    headers = {
        "Accept": "application/sparql-results+json",
        "Content-Type": "application/sparql-query"
    }
    
    try:
        response = requests.post(
            f"{endpoint}?timeout=60000",
            data=sparql_query,
            headers=headers,
            timeout=120  # 2 minute timeout
        )
        response.raise_for_status()
        return response.json()
    except requests.exceptions.RequestException as e:
        logger.error(f"SPARQL request failed: {e}")
        raise Exception(f"Failed to fetch chemical data: {e}")

@celery.task(bind=True, queue='SPARQL')
def sparql_task(self, payload):
    """Fetch chemical data from SPARQL endpoint and upload as JSON file to GCS."""
    logger.info(f"=== TASK STARTED: sparql_task ===")
    logger.info(f"Task ID: {self.request.id}")
    logger.info(f"Payload: {payload}")
    
    try:
        task_id = payload.get("task_id")
        user_id = payload.get("user_id")
        user_query = payload.get("user_query", "")

        if not all([task_id, user_id]):
            raise ValueError(f"Missing required fields. task_id={task_id}, user_id={user_id}")

        logger.info(f"Starting chemical fetcher task for user {user_id}")
        emit_status(task_id, "parsing query")
        
        # Parse user query to extract limit using AI
        limit = extract_limit_with_ai(user_query)
        
        emit_status(task_id, f"fetching {limit} chemicals from database")
        
        # Fetch chemical data
        chemical_data = fetch_chemicals_sparql(limit)
        
        emit_status(task_id, "processing results")
        
        # Process the results for better formatting
        processed_results = {
            "metadata": {
                "query": user_query,
                "limit_requested": limit,
                "total_results": len(chemical_data.get("results", {}).get("bindings", [])),
                "endpoint": "https://neptune.biobricks.ai/sparql",
                "timestamp": str(uuid.uuid4())
            },
            "chemicals": []
        }
        
        # Extract and format chemical data
        for binding in chemical_data.get("results", {}).get("bindings", []):
            chemical = {
                "name": binding.get("name", {}).get("value", "Unknown"),
                "cas_number": binding.get("cas", {}).get("value", "").split("/")[-1] if binding.get("cas") else None,
                "dsstox_id": binding.get("dsstox", {}).get("value", "").split("/")[-1] if binding.get("dsstox") else None,
                "dsstox_url": binding.get("dsstox", {}).get("value", "") if binding.get("dsstox") else None,
                "cas_url": binding.get("cas", {}).get("value", "") if binding.get("cas") else None
            }
            processed_results["chemicals"].append(chemical)
        
        emit_status(task_id, "sending results message")
        
        # Create summary message for user
        summary_message = f"**Query:** {user_query}\n\n"
        summary_message += f"**Results:** {len(processed_results['chemicals'])} chemicals found from the BioBricks database\n\n\n"
        summary_message += f"**Data includes:** Chemical names, CAS numbers, and DSSTOX identifiers\n\n\n"
        
        # Add preview of first 5 chemicals with improved formatting
        summary_message += "**Fetched Data:**\n\n"
        for i, chemical in enumerate(processed_results['chemicals'][:limit]):
            summary_message += f"**{i+1}.** *{chemical['name']}*\n"
            if chemical['cas_number']:
                summary_message += f"   • **CAS:** `{chemical['cas_number']}`\n"
            if chemical['dsstox_id']:
                summary_message += f"   • **DSSTOX:** `{chemical['dsstox_id']}`\n"
            summary_message += "\n"
        
        summary_message += f"The complete dataset ({len(processed_results['chemicals'])} chemicals) has been saved as a JSON file for download."

        emit_status(task_id, "uploading JSON file to GCS")

        # Create JSON file
        json_filename = f"chemical_data_{uuid.uuid4().hex}.json"
        
        # Upload JSON file to GCS
        with tempfile.NamedTemporaryFile(mode='w', suffix='.json', delete=False, encoding='utf-8') as temp_json_file:
            temp_json_path = temp_json_file.name
            json.dump(processed_results, temp_json_file, indent=2, ensure_ascii=False)
        
        try:
            # Upload JSON file to GCS
            json_gcs_path = f"tasks/{task_id}/{json_filename}"
            upload_local_file_to_gcs(Path(temp_json_path), json_gcs_path, content_type='application/json')
            emit_status(task_id, "files uploaded")

            # Emit JSON file event
            json_file_data = {
                "user_id": user_id,
                "filename": json_filename,
                "filepath": json_gcs_path,
                "file_type": "json",
                "content_type": "application/json"
            }
            emit_task_file(task_id, json_file_data)
            
            # Add download link to the main message and send it
            # summary_message += f"\n\n**Download Link:** [Download JSON File]({json_gcs_path})"
            message = MessageSchema(role="assistant", content=summary_message)
            emit_task_message(task_id, message.model_dump())
            
        finally:
            # Clean up temporary files
            try:
                os.unlink(temp_json_path)
            except Exception as e:
                logger.warning(f"Failed to clean up temp files: {e}")

        logger.info(f"Task {task_id} completed successfully")
        finished_at = Task.mark_finished(task_id)
        emit_status(task_id, "done")
        return {"done": True, "finished_at": finished_at}

    except Exception as e:
        logger.error(f"Error in sparql_task: {str(e)}", exc_info=True)
        emit_status(task_id, "error")
        raise  # Re-raise the exception so Celery knows the task failed

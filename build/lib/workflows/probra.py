import redis
import json
import os
import uuid
import logging
import tempfile
from datetime import datetime
from workflows.celery_worker import celery
from webserver.model.message import MessageSchema
import hashlib
from webserver.model.task import Task
from webserver.storage import GCSFileStorage
from webserver.cache_manager import cache_manager

from webserver.tools.deeptox_agent import deeptox_agent
from webserver.tools.toxicity_models import ChemicalToxicityAssessment
from webserver.ai_service import convert_pydantic_to_markdown

logging.getLogger().info("probra_gcs.py module loaded")
logger = logging.getLogger(__name__)



def emit_status(task_id, status):
    logger.info(f"[emit_status] {task_id} -> {status}")
    Task.set_status(task_id, status)
    r = redis.Redis(
        host=os.environ.get("REDIS_HOST", "localhost"),
        port=int(os.environ.get("REDIS_PORT", "6379"))
    )
    task = Task.get_task(task_id)
    event = {
        "type": "task_status_update",
        "task_id": task.task_id,
        "data": task.to_dict(),
    }
    r.publish("celery_updates", json.dumps(event, default=str))

@celery.task(bind=True)
def probra_task(self, payload):
    """GCS-enabled background task that emits progress messages and uploads files to GCS."""
    try:
        logger.info(f"Starting probra task with payload: {payload}")
        r = redis.Redis(
            host=os.environ.get("REDIS_HOST", "localhost"),
            port=int(os.environ.get("REDIS_PORT", "6379"))
        )
        task_id = payload.get("task_id")
        user_id = payload.get("user_id")

        if not all([task_id, user_id]):
            raise ValueError(f"Missing required fields. task_id={task_id}, user_id={user_id}")

        emit_status(task_id, "starting")
        chemprop_query = payload.get("payload", "Is Gentamicin nephrotoxic?")
        logger.info(f"User query for deeptox_agent: {chemprop_query}")

        emit_status(task_id, "checking cache")
        # Caching: use hash of query as key for both structured and markdown content
        cache_key_base = f"probra_cache:{hashlib.sha256(chemprop_query.encode()).hexdigest()}"
        structured_cache_key = f"{cache_key_base}:structured"
        markdown_cache_key = f"{cache_key_base}:markdown"
        
        cached_structured = r.get(structured_cache_key)
        cached_markdown = r.get(markdown_cache_key)
        
        if cached_structured and cached_markdown:
            # Use cached data
            original_response = json.loads(cached_structured.decode())
            response_content = cached_markdown.decode()
            emit_status(task_id, "using cache")
        else:
            # Run agent and generate both formats
            emit_status(task_id, "running agent")
            response = deeptox_agent.run(chemprop_query)
            
            # Get original structured response
            if isinstance(response.content, ChemicalToxicityAssessment):
                original_response = response.content.model_dump()
            elif isinstance(response.content, dict):
                original_response = response.content
            else:
                # Fallback for string response
                original_response = {"content": response.content}
            
            # Convert to markdown
            emit_status(task_id, "converting to markdown")
            response_content = convert_pydantic_to_markdown(
                original_response, 
                chemprop_query
            )
            
            # Cache both formats
            r.set(structured_cache_key, json.dumps(original_response), ex=60*60*24)
            r.set(markdown_cache_key, response_content, ex=60*60*24)
            emit_status(task_id, "agent complete")

        emit_status(task_id, "sending message")
        # display raw markdown content directly to user
        message = MessageSchema(role="assistant", content=response_content)
        event = {
            "type": "task_message",
            "data": message.model_dump(),
            "task_id": task_id,
        }
        r.publish("celery_updates", json.dumps(event, default=str))
        logger.info(f"Published task_message event to Redis for task_id={task_id}")

        emit_status(task_id, "uploading files to GCS")
        
        # Get the original structured data for JSON file
        # original_response is already set from the cache or agent response above
        
        # Create JSON file
        json_filename = f"probra_result_{uuid.uuid4().hex}.json"
        json_content = ""
        if isinstance(original_response, ChemicalToxicityAssessment):
            json_content = json.dumps(original_response.model_dump(), indent=2, ensure_ascii=False)
        elif isinstance(original_response, dict):
            json_content = json.dumps(original_response, indent=2, ensure_ascii=False)
        else:
            # Fallback: create minimal JSON
            json_content = json.dumps({"message": "Analysis completed", "content": str(original_response)}, indent=2)
        
        # Create Markdown file
        md_filename = f"probra_result_{uuid.uuid4().hex}.md"
        
        # Upload JSON file to GCS
        with tempfile.NamedTemporaryFile(mode='w', suffix='.json', delete=False, encoding='utf-8') as temp_json_file:
            temp_json_path = temp_json_file.name
            temp_json_file.write(json_content)
        
        # Upload Markdown file to GCS
        with tempfile.NamedTemporaryFile(mode='w', suffix='.md', delete=False, encoding='utf-8') as temp_md_file:
            temp_md_path = temp_md_file.name
            temp_md_file.write(response_content)
        
        try:
            # Upload to GCS
            gcs_storage = GCSFileStorage()
            
            # Upload JSON file
            json_gcs_path = f"tasks/{task_id}/{json_filename}"
            gcs_storage.upload_file(temp_json_path, json_gcs_path, content_type='application/json')
            cache_manager.cache_file_content(json_gcs_path, json_content)
            
            # Upload Markdown file
            md_gcs_path = f"tasks/{task_id}/{md_filename}"
            gcs_storage.upload_file(temp_md_path, md_gcs_path, content_type='text/markdown')
            cache_manager.cache_file_content(md_gcs_path, response_content)
            
            logger.info(f"Files uploaded to GCS: {json_gcs_path}, {md_gcs_path}")
            
            emit_status(task_id, "files uploaded")
            
            # Publish JSON file event
            json_file_event = {
                "type": "task_file",
                "task_id": task_id,
                "data": {
                    "user_id": user_id,
                    "filename": json_filename,
                    "filepath": json_gcs_path,
                    "file_type": "json",
                    "content_type": "application/json"
                },
            }
            r.publish("celery_updates", json.dumps(json_file_event, default=str))
            
            # Publish Markdown file event
            md_file_event = {
                "type": "task_file",
                "task_id": task_id,
                "data": {
                    "user_id": user_id,
                    "filename": md_filename,
                    "filepath": md_gcs_path,
                    "file_type": "markdown",
                    "content_type": "text/markdown"
                },
            }
            r.publish("celery_updates", json.dumps(md_file_event, default=str))
            
        finally:
            # Clean up temporary files
            try:
                os.unlink(temp_json_path)
                os.unlink(temp_md_path)
            except Exception as e:
                logger.warning(f"Failed to clean up temp files: {e}")

        finished_at = Task.mark_finished(task_id)
        emit_status(task_id, "done")
        return {"done": True, "finished_at": finished_at}

    except Exception as e:
        logger.error(f"Error in probra_task: {str(e)}", exc_info=True)
        emit_status(task_id, "error")
        raise  # Re-raise the exception so Celery knows the task failed 
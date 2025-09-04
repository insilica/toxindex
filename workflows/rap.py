import redis
import json
import os
import uuid
import logging
import tempfile
from datetime import datetime
from workflows.celery_app import celery
from webserver.model.message import MessageSchema
import hashlib
from webserver.model.task import Task
from webserver.storage import GCSFileStorage
from webserver.cache_manager import cache_manager

from webserver.tools.deeptox_agent import deeptox_agent
from webserver.tools.toxicity_models_simple import ChemicalToxicityAssessment
from webserver.ai_service import convert_pydantic_to_markdown
from workflows.utils import emit_status, publish_to_celery_updates, publish_to_socketio, get_redis_connection

logging.getLogger().info("rap.py module loaded")
logger = logging.getLogger(__name__)




def emit_task_message(task_id, message_data):
    """Emit task message to both database and real-time channels"""
    # Publish to celery_updates for database processing
    publish_to_celery_updates("task_message", task_id, message_data)
    
    # Publish to Socket.IO for real-time updates
    task = Task.get_task(task_id)
    if task and getattr(task, 'session_id', None):
        # Emit to chat session room
        publish_to_socketio("new_message", f"chat_session_{task.session_id}", message_data)
    
    # Emit to task room
    publish_to_socketio("task_message", f"task_{task_id}", {
        "type": "task_message",
        "data": message_data,
        "task_id": str(task_id),  # Convert UUID to string for JSON serialization
    })


def emit_task_file(task_id, file_data):
    """Emit task file to both database and real-time channels"""
    # Publish to celery_updates for database processing
    publish_to_celery_updates("task_file", task_id, file_data)
    
    # Publish to Socket.IO for real-time updates
    publish_to_socketio("task_file", f"task_{task_id}", file_data)


@celery.task(bind=True, queue='rap')
def rap_task(self, payload):
    """GCS-enabled background task that emits progress messages and uploads files to GCS."""
    # Add detailed task start logging
    logger.info(f"=== TASK STARTED: rap_task ===")
    logger.info(f"Task ID: {self.request.id}")
    logger.info(f"Payload: {payload}")
    
    try:
        logger.info(f"Starting rap task with payload: {payload}")
        r = get_redis_connection()
        task_id = payload.get("task_id")
        user_id = payload.get("user_id")

        if not all([task_id, user_id]):
            raise ValueError(f"Missing required fields. task_id={task_id}, user_id={user_id}")

        logger.info(f"Processing task {task_id} for user {user_id}")
        emit_status(task_id, "starting")
        chemprop_query = payload.get("user_query", "Is Gentamicin nephrotoxic?")
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
            logger.info(f"Using cached data for query: {chemprop_query}")
            try:
                original_response = json.loads(cached_structured.decode())
                response_content = cached_markdown.decode()
                logger.info(f"Cached data loaded successfully. Type: {type(original_response)}")
                emit_status(task_id, "using cache")
            except Exception as e:
                logger.error(f"Error loading cached data: {e}")
                # Clear corrupted cache and regenerate
                r.delete(structured_cache_key)
                r.delete(markdown_cache_key)
                cached_structured = None
                cached_markdown = None
        else:
            # Run agent and generate both formats
            logger.info(f"Running agent for query: {chemprop_query}")
            emit_status(task_id, "running agent")
            response = deeptox_agent.run(chemprop_query)
            
            # Get original structured response and ensure it's always a dict
            if isinstance(response.content, ChemicalToxicityAssessment):
                original_response = response.content.model_dump()
            elif isinstance(response.content, dict):
                original_response = response.content
            else:
                # Fallback for string response
                original_response = {"content": response.content}
            
            # Convert to markdown
            logger.info(f"Converting response to markdown")
            emit_status(task_id, "converting to markdown")
            
            # original_response is already a dict from earlier processing
            logger.info(f"Converting to markdown. Input type: {type(original_response)}")
            response_content = convert_pydantic_to_markdown(
                original_response, 
                chemprop_query
            )
            
            # Cache both formats
            logger.info(f"Caching results for future use")
            r.set(structured_cache_key, json.dumps(original_response, default=str), ex=60*60*24)
            r.set(markdown_cache_key, response_content, ex=60*60*24)
            emit_status(task_id, "agent complete")

        logger.info(f"Sending message to user")
        emit_status(task_id, "sending message")
        # display raw markdown content directly to user
        message = MessageSchema(role="assistant", content=response_content)
        emit_task_message(task_id, message.model_dump())

        logger.info(f"Uploading files to GCS")
        emit_status(task_id, "uploading files to GCS")
        
        # Get the original structured data for JSON file
        # original_response is already set from the cache or agent response above
        
        # Create JSON file
        json_filename = f"probra_result_{uuid.uuid4().hex}.json"
        json_content = ""
        
        # original_response is already a dict from earlier processing
        json_content = json.dumps(original_response, indent=2, ensure_ascii=False, default=str)
        
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
            
            # Emit JSON file event
            json_file_data = {
                "user_id": user_id,
                "filename": json_filename,
                "filepath": json_gcs_path,
                "file_type": "json",
                "content_type": "application/json"
            }
            emit_task_file(task_id, json_file_data)
            
            # Emit Markdown file event
            md_file_data = {
                "user_id": user_id,
                "filename": md_filename,
                "filepath": md_gcs_path,
                "file_type": "markdown",
                "content_type": "text/markdown"
            }
            emit_task_file(task_id, md_file_data)
            
        finally:
            # Clean up temporary files
            try:
                os.unlink(temp_json_path)
                os.unlink(temp_md_path)
            except Exception as e:
                logger.warning(f"Failed to clean up temp files: {e}")

        logger.info(f"Task {task_id} completed successfully")
        finished_at = Task.mark_finished(task_id)
        emit_status(task_id, "done")
        return {"done": True, "finished_at": finished_at}

    except Exception as e:
        logger.error(f"Error in rap_task: {str(e)}", exc_info=True)
        emit_status(task_id, "error")
        raise  # Re-raise the exception so Celery knows the task failed 
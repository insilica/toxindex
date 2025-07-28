import redis
import json
import os
import uuid
import logging
import pydantic
import tempfile
from pathlib import Path
from workflows.celery_worker import celery
from webserver.model.message import MessageSchema
import hashlib
from webserver.model.task import Task
from webserver.storage import GCSFileStorage
from webserver.cache_manager import cache_manager

from RAP.tool_deeptox import deeptox_agent

logging.getLogger().info("probra_gcs.py module loaded")
logger = logging.getLogger(__name__)

def get_pydantic_serializer(obj):
    """Get the appropriate serialization method based on Pydantic version."""
    if pydantic.__version__.startswith('2'):
        return obj.model_dump()
    return obj.dict()

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
        # Caching: use hash of query as key
        cache_key = f"probra_cache:{hashlib.sha256(chemprop_query.encode()).hexdigest()}"
        cached = r.get(cache_key)
        if cached:
            response_content = cached.decode()
            emit_status(task_id, "using cache")
        else:
            emit_status(task_id, "running agent")
            response = deeptox_agent.run(chemprop_query, history=True)
            response_content = response.content
            r.set(cache_key, response_content, ex=60*60*24)
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

        emit_status(task_id, "uploading file to GCS")
        # Create temporary file and upload to GCS
        tmp_filename = f"probra_result_{uuid.uuid4().hex}.md"
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.md', delete=False, encoding='utf-8') as temp_file:
            temp_path = temp_file.name
            temp_file.write(response_content)
        
        try:
            # Upload to GCS
            gcs_storage = GCSFileStorage()
            gcs_path = f"tasks/{task_id}/{tmp_filename}"
            
            # Upload file to GCS
            gcs_storage.upload_file(temp_path, gcs_path, content_type='text/markdown')
            
            # Cache the content for future access
            cache_manager.cache_file_content(gcs_path, response_content)
            
            logger.info(f"File uploaded to GCS: {gcs_path}")
            
            emit_status(task_id, "file uploaded")
            # Publish file event with GCS path
            file_event = {
                "type": "task_file",
                "task_id": task_id,
                "data": {
                    "user_id": user_id,
                    "filename": tmp_filename,
                    "filepath": gcs_path,  # Use GCS path instead of local path
                },
            }
            
            r.publish("celery_updates", json.dumps(file_event, default=str))
            
        finally:
            # Clean up temporary file
            try:
                os.unlink(temp_path)
            except Exception as e:
                logger.warning(f"Failed to clean up temp file {temp_path}: {e}")

        finished_at = Task.mark_finished(task_id)
        emit_status(task_id, "done")
        return {"done": True, "finished_at": finished_at}

    except Exception as e:
        logger.error(f"Error in probra_task: {str(e)}", exc_info=True)
        emit_status(task_id, "error")
        raise  # Re-raise the exception so Celery knows the task failed 
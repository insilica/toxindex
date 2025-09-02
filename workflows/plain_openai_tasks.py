import redis
import json
import os
import uuid
import logging
import openai
import hashlib
import tempfile
# from pathlib import Path
from workflows.celery_app import celery
from webserver.model.message import MessageSchema
from webserver.tools.toxicity_schema import TOXICITY_SCHEMA
from webserver.model.task import Task
from webserver.storage import GCSFileStorage
from webserver.cache_manager import cache_manager

logger = logging.getLogger(__name__)

def emit_status(task_id, status):
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

@celery.task(bind=True, queue='openai')
def plain_openai_task(self, payload):
    """GCS-enabled OpenAI task that uploads results to GCS."""
    try:
        logger.info(f"Running plain_openai_task for task_id={payload.get('task_id')}, user_id={payload.get('user_id')}, payload={payload}")
        r = redis.Redis(
            host=os.environ.get("REDIS_HOST", "localhost"),
            port=int(os.environ.get("REDIS_PORT", "6379"))
        )
        task_id = payload.get("task_id")
        user_id = payload.get("user_id")
        user_query = payload.get("payload", "Is Gentamicin nephrotoxic?")
        model = "gpt-4"
        emit_status(task_id, "checking cache")
        cache_key = f"openai_cache:{model}:{hashlib.sha256(user_query.encode()).hexdigest()}"
        cached = r.get(cache_key)
        if cached:
            response_data = json.loads(cached)
            content = response_data["content"]
            emit_status(task_id, "using cache")
        else:
            emit_status(task_id, "calling openai")
            client = openai.OpenAI()
            response = client.chat.completions.create(
                model=model,
                messages=[{"role": "user", "content": user_query}],
                temperature=0.7,
            )
            content = response.choices[0].message.content
            response_data = {"content": content}
            r.set(cache_key, json.dumps(response_data), ex=60*60*24)
            emit_status(task_id, "openai complete")

        emit_status(task_id, "publishing message")
        # display raw markdown content directly to user
        message = MessageSchema(role="assistant", content=content)
        event = {
            "type": "task_message",
            "data": message.model_dump(),
            "task_id": task_id,
        }
        logger.info(f"[plain_openai_task] Publishing message event: {event}")
        r.publish("celery_updates", json.dumps(event, default=str))

        emit_status(task_id, "uploading file to GCS")
        # Create temporary file and upload to GCS
        tmp_filename = f"plain_openai_result_{uuid.uuid4().hex}.md"
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.md', delete=False, encoding='utf-8') as temp_file:
            temp_path = temp_file.name
            temp_file.write(content)
        
        try:
            # Upload to GCS
            gcs_storage = GCSFileStorage()
            gcs_path = f"tasks/{task_id}/{tmp_filename}"
            
            # Upload file to GCS
            gcs_storage.upload_file(temp_path, gcs_path, content_type='text/markdown')
            
            # Cache the content for future access
            cache_manager.cache_file_content(gcs_path, content)
            
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
            logger.info(f"[plain_openai_task] Publishing file event: {file_event}")
            r.publish("celery_updates", json.dumps(file_event, default=str))
            
        finally:
            # Clean up temporary file
            try:
                os.unlink(temp_path)
            except Exception as e:
                logger.warning(f"Failed to clean up temp file {temp_path}: {e}")
        
        logger.info("plain_openai_task completed successfully")

        finished_at = Task.mark_finished(task_id)
        emit_status(task_id, "done")
        return {
            "done": True,
            "finished_at": finished_at,
        }
    except Exception as e:
        logger.error(f"Error in plain_openai_task: {str(e)}", exc_info=True)
        emit_status(task_id, "error")
        raise

@celery.task(bind=True)
def openai_json_schema_task(self, payload):
    """GCS-enabled OpenAI JSON schema task that uploads results to GCS."""
    try:
        logger.info(f"[ToxJson GCS] Running openai_json_schema_task for task_id={payload.get('task_id')}, user_id={payload.get('user_id')}, payload={payload}")
        r = redis.Redis(
            host=os.environ.get("REDIS_HOST", "localhost"),
            port=int(os.environ.get("REDIS_PORT", "6379"))
        )
        task_id = payload.get("task_id")
        user_id = payload.get("user_id")
        user_query = payload.get("payload", "Is Gentamicin nephrotoxic?")
        model = "gpt-4"
        schema_str = json.dumps(TOXICITY_SCHEMA, indent=2)
        system_prompt = (
            "You are a toxicology research assistant. "
            "Answer the user's question about chemical toxicity in JSON format matching this schema:\n"
            f"{schema_str}\n\n"
            "Provide a detailed analysis with evidence from scientific literature."
        )
        
        emit_status(task_id, "checking cache")
        cache_key = f"openai_json_cache:{model}:{hashlib.sha256(user_query.encode()).hexdigest()}"
        cached = r.get(cache_key)
        if cached:
            response_data = json.loads(cached)
            content = response_data["content"]
            emit_status(task_id, "using cache")
        else:
            emit_status(task_id, "calling openai")
            client = openai.OpenAI()
            response = client.chat.completions.create(
                model=model,
                messages=[
                    {"role": "system", "content": system_prompt},
                    {"role": "user", "content": user_query}
                ],
                temperature=0.7,
            )
            content = response.choices[0].message.content
            response_data = {"content": content}
            r.set(cache_key, json.dumps(response_data), ex=60*60*24)
            emit_status(task_id, "openai complete")

        emit_status(task_id, "publishing message")
        message = MessageSchema(role="assistant", content=content)
        event = {
            "type": "task_message",
            "data": message.model_dump(),
            "task_id": task_id,
        }
        logger.info(f"[openai_json_schema_task] Publishing message event: {event}")
        r.publish("celery_updates", json.dumps(event, default=str))

        emit_status(task_id, "uploading file to GCS")
        # Create temporary file and upload to GCS
        tmp_filename = f"openai_json_result_{uuid.uuid4().hex}.json"
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.json', delete=False, encoding='utf-8') as temp_file:
            temp_path = temp_file.name
            temp_file.write(content)
        
        try:
            # Upload to GCS
            gcs_storage = GCSFileStorage()
            gcs_path = f"tasks/{task_id}/{tmp_filename}"
            
            # Upload file to GCS
            gcs_storage.upload_file(temp_path, gcs_path, content_type='application/json')
            
            # Cache the content for future access
            cache_manager.cache_file_content(gcs_path, content)
            
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
            logger.info(f"[openai_json_schema_task] Publishing file event: {file_event}")
            r.publish("celery_updates", json.dumps(file_event, default=str))
            
        finally:
            # Clean up temporary file
            try:
                os.unlink(temp_path)
            except Exception as e:
                logger.warning(f"Failed to clean up temp file {temp_path}: {e}")
        
        logger.info("openai_json_schema_task completed successfully")

        finished_at = Task.mark_finished(task_id)
        emit_status(task_id, "done")
        return {
            "done": True,
            "finished_at": finished_at,
        }
    except Exception as e:
        logger.error(f"Error in openai_json_schema_task: {str(e)}", exc_info=True)
        emit_status(task_id, "error")
        raise 
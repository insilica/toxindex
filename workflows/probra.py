import redis
import json
import os
import uuid
import logging
import pydantic
from workflows.celery_worker import celery
from webserver.model.message import MessageSchema
import hashlib
from webserver.model.task import Task
from webserver.data_paths import CHATS_ROOT

from RAP.tool_deeptox import deeptox_agent

logging.getLogger().info("probra.py module loaded")
logger = logging.getLogger(__name__)

def get_pydantic_serializer(obj):
    """Get the appropriate serialization method based on Pydantic version."""
    if pydantic.__version__.startswith('2'):
        return obj.model_dump()
    return obj.dict()

def emit_status(task_id, status):
    logger.info(f"[emit_status] {task_id} -> {status}")
    Task.set_status(task_id, status)
    r = redis.Redis()
    task = Task.get_task(task_id)
    event = {
        "type": "task_status_update",
        "task_id": task.task_id,
        "data": task.to_dict(),
    }
    r.publish("celery_updates", json.dumps(event, default=str))

@celery.task(bind=True)
# async def probra_task(self, payload):
def probra_task(self, payload):
    """Example background task that emits progress messages and uploads a file."""
    try:
        logger.info(f"Starting probra task with payload: {payload}")
        r = redis.Redis()
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

        tmp_filename = f"probra_result_{uuid.uuid4().hex}.md"
        project_tmp_dir = CHATS_ROOT()
        os.makedirs(project_tmp_dir, exist_ok=True)
        tmp_path = os.path.join(project_tmp_dir, tmp_filename)

        logger.info(f"File creating in tmp: {tmp_path} for task {task_id}")
        with open(tmp_path, 'w', encoding='utf-8') as f:
            f.write(response_content)
        logger.info(f"File created in tmp: {tmp_path} for task {task_id}")

        emit_status(task_id, "file created")
        # S3 upload removed; only local file creation remains
        file_event = {
            "type": "task_file",
            "task_id": task_id,
            "data": {
                "user_id": user_id,
                "filename": tmp_filename,
                "filepath": tmp_path,  # Use filepath for local file path
            },
        }
        
        r.publish("celery_updates", json.dumps(file_event, default=str))

        finished_at = Task.mark_finished(task_id)
        emit_status(task_id, "done")
        return {"done": True, "finished_at": finished_at}

    except Exception as e:
        logger.error(f"Error in probra_task: {str(e)}", exc_info=True)
        emit_status(task_id, "error")
        raise  # Re-raise the exception so Celery knows the task failed

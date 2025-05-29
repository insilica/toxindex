import redis
import json
import time
import os
import uuid
import logging
import pydantic
from workflows.celery_worker import celery
from webserver.model.message import MessageSchema
from webserver import S3FileStorage

logger = logging.getLogger(__name__)

def get_pydantic_serializer(obj):
    """Get the appropriate serialization method based on Pydantic version."""
    if pydantic.__version__.startswith('2'):
        return obj.model_dump()
    return obj.dict()

@celery.task(bind=True)
def probra_task(self, payload):
    """Example background task that emits progress messages and uploads a file."""
    try:
        logger.info(f"Starting probra task with payload: {payload}")
        r = redis.Redis()
        task_id = payload.get("task_id")
        user_id = payload.get("user_id")

        if not all([task_id, user_id]):
            raise ValueError(f"Missing required fields. task_id={task_id}, user_id={user_id}")

        for i in range(5):
            logger.debug(f"Processing step {i}")
            message = MessageSchema(role="assistant", content=f"Step {i}")
            event = {
                "type": "task_message",
                "data": get_pydantic_serializer(message),
                "task_id": task_id,
            }
            r.publish("celery_updates", json.dumps(event))
            time.sleep(1)

        # Create a simple result file and upload it to S3
        tmp_filename = f"probra_result_{uuid.uuid4().hex}.txt"
        tmp_path = os.path.join("/tmp", tmp_filename)
        logger.info(f"Creating temporary file: {tmp_path}")
        print(f"Creating temporary file: {tmp_path}")
        with open(tmp_path, "w", encoding="utf-8") as f:
            f.write("ProbRA task complete.\n")

        storage = S3FileStorage()
        s3_key = storage.upload_file(tmp_path, f"{task_id}/{tmp_filename}", "text/plain")
        download_url = storage.generate_download_url(s3_key)
        logger.info(f"File uploaded to S3: {s3_key}")

        file_event = {
            "type": "task_file",
            "task_id": task_id,
            "data": {
                "user_id": user_id,
                "filename": tmp_filename,
                "filepath": s3_key,
                "s3_url": download_url,
            },
        }
        r.publish("celery_updates", json.dumps(file_event))
        logger.info("Task completed successfully")
        return {"done": True}

    except Exception as e:
        logger.error(f"Error in probra_task: {str(e)}", exc_info=True)
        raise  # Re-raise the exception so Celery knows the task failed

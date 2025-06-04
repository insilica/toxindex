import redis
import json
import time
import os
import uuid
import logging
import pydantic
from workflows.celery_worker import celery
from webserver.model.message import MessageSchema
from webserver.storage import S3FileStorage

from RAP.tool_deeptox import deeptox_agent
from RAP.toxicity_schema import TOXICITY_SCHEMA

logger = logging.getLogger(__name__)

def get_pydantic_serializer(obj):
    """Get the appropriate serialization method based on Pydantic version."""
    if pydantic.__version__.startswith('2'):
        return obj.model_dump()
    return obj.dict()

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

        chemprop_query = payload.get("payload", "Is Gentamicin nephrotoxic?")
        logger.info(f"User query for deeptox_agent: {chemprop_query}")

        # async for chunk in await deeptox_agent.arun(chemprop_query, stream=True):
        #     message = MessageSchema(role="assistant", content=chunk)
        #     event = {
        #         "type": "task_message",
        #         "data": message.model_dump(),
        #         "task_id": task_id,
        #     }
        #     r.publish("celery_updates", json.dumps(event))

        response = deeptox_agent.run(chemprop_query, history=True)
        # display raw markdown content directly to user
        message = MessageSchema(role="assistant", content=response.content)
        event = {
            "type": "task_message",
            "data": message.model_dump(),
            "task_id": task_id,
        }
        r.publish("celery_updates", json.dumps(event))

        
        # the tmp path is within nix environment, not local project folder.
        tmp_filename = f"probra_result_{uuid.uuid4().hex}.md"
        project_tmp_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'tmp'))
        os.makedirs(project_tmp_dir, exist_ok=True)
        tmp_path = os.path.join(project_tmp_dir, tmp_filename)

        # Save the response based on its type
        with open(tmp_path, 'w', encoding='utf-8') as f:
            f.write(response.content)

        # Create a simple result file and upload it to S3
        logger.info(f"Creating temporary file: {tmp_path}")
        

        # storage = S3FileStorage()
        # s3_key = storage.upload_file(tmp_path, f"{task_id}/{tmp_filename}", "text/plain")
        # download_url = storage.generate_download_url(s3_key)
        # logger.info(f"File uploaded to S3: {s3_key}")

        # file_event = {
        #     "type": "task_file",
        #     "task_id": task_id,
        #     "data": {
        #         "user_id": user_id,
        #         "filename": tmp_filename,
        #         "filepath": s3_key,
        #         "s3_url": download_url,
        #     },
        # }
        # r.publish("celery_updates", json.dumps(file_event))
        logger.info("Task completed successfully")
        return {"done": True}

    except Exception as e:
        logger.error(f"Error in probra_task: {str(e)}", exc_info=True)
        raise  # Re-raise the exception so Celery knows the task failed

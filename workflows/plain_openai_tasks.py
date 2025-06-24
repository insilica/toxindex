import redis
import json
import os
import uuid
import logging
import openai
import hashlib
from workflows.celery_worker import celery
from webserver.model.message import MessageSchema
from webserver.storage import S3FileStorage
from RAP.toxicity_schema import TOXICITY_SCHEMA
# Update finished_at in the database
import datetime
import webserver.datastore as ds
from webserver.model.task import Task

logger = logging.getLogger(__name__)

@celery.task(bind=True)
def plain_openai_task(self, payload):
    """Call OpenAI with a plain prompt and publish the result."""
    try:
        logger.info(f"[ToxDirect] Running plain_openai_task for task_id={payload.get('task_id')}, user_id={payload.get('user_id')}, payload={payload}")
        r = redis.Redis()
        task_id = payload.get("task_id")
        user_id = payload.get("user_id")
        user_query = payload.get("payload", "Is Gentamicin nephrotoxic?")
        model = "gpt-4"
        cache_key = f"openai_cache:{model}:{hashlib.sha256(user_query.encode()).hexdigest()}"
        cached = r.get(cache_key)
        if cached:
            response_data = json.loads(cached)
            content = response_data["content"]
        else:
            client = openai.OpenAI()
            response = client.chat.completions.create(
                model=model,
                messages=[{"role": "user", "content": user_query}],
                temperature=0.7,
            )
            content = response.choices[0].message.content
            response_data = {"content": content}
            r.set(cache_key, json.dumps(response_data), ex=60*60*24)

        message = MessageSchema(role="assistant", content=content)
        event = {
            "type": "task_message",
            "data": message.model_dump(),
            "task_id": task_id,
        }
        logger.info(f"[plain_openai_task] Publishing message event: {event}")
        r.publish("celery_updates", json.dumps(event, default=str))

        tmp_filename = f"plain_openai_result_{uuid.uuid4().hex}.md"
        project_tmp_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'tmp'))
        os.makedirs(project_tmp_dir, exist_ok=True)
        tmp_path = os.path.join(project_tmp_dir, tmp_filename)
        with open(tmp_path, 'w', encoding='utf-8') as f:
            f.write(content)
        storage = S3FileStorage()
        s3_key = storage.upload_file(tmp_path, f"{task_id}/{tmp_filename}", "text/plain")
        download_url = storage.generate_download_url(s3_key)
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
        logger.info(f"[plain_openai_task] Publishing file event: {file_event}")
        r.publish("celery_updates", json.dumps(file_event, default=str))
        logger.info("plain_openai_task completed successfully")

        finished_at = Task.mark_finished(task_id)

        return {
            "done": True,
            "finished_at": finished_at,
        }
    except Exception as e:
        logger.error(f"Error in plain_openai_task: {str(e)}", exc_info=True)
        raise

@celery.task(bind=True)
def openai_json_schema_task(self, payload):
    """Call OpenAI with a prompt instructing output as JSON matching TOXICITY_SCHEMA."""
    try:
        logger.info(f"[ToxJson] Running openai_json_schema_task for task_id={payload.get('task_id')}, user_id={payload.get('user_id')}, payload={payload}")
        r = redis.Redis()
        task_id = payload.get("task_id")
        user_id = payload.get("user_id")
        user_query = payload.get("payload", "Is Gentamicin nephrotoxic?")
        model = "gpt-4"
        schema_str = json.dumps(TOXICITY_SCHEMA, indent=2)
        system_prompt = (
            "You are a toxicology research assistant. "
            "Answer the user's question as a JSON object matching the following schema. "
            "Strictly adhere to the structure and required fields.\n\nSchema:\n" + schema_str
        )
        cache_key = f"openai_json_schema_cache:{model}:{hashlib.sha256((user_query+system_prompt).encode()).hexdigest()}"
        cached = r.get(cache_key)
        if cached:
            response_data = json.loads(cached)
            content = response_data["content"]
        else:
            client = openai.OpenAI()
            response = client.chat.completions.create(
                model=model,
                messages=[
                    {"role": "system", "content": system_prompt},
                    {"role": "user", "content": user_query},
                ],
                temperature=0.0,
            )
            content = response.choices[0].message.content
            response_data = {"content": content}
            r.set(cache_key, json.dumps(response_data), ex=60*60*24)

        message = MessageSchema(role="assistant", content=content)
        event = {
            "type": "task_message",
            "data": message.model_dump(),
            "task_id": task_id,
        }
        r.publish("celery_updates", json.dumps(event, default=str))

        tmp_filename = f"openai_json_schema_result_{uuid.uuid4().hex}.json"
        project_tmp_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'tmp'))
        os.makedirs(project_tmp_dir, exist_ok=True)
        tmp_path = os.path.join(project_tmp_dir, tmp_filename)
        with open(tmp_path, 'w', encoding='utf-8') as f:
            f.write(content)
        storage = S3FileStorage()
        s3_key = storage.upload_file(tmp_path, f"{task_id}/{tmp_filename}", "application/json")
        download_url = storage.generate_download_url(s3_key)
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
        r.publish("celery_updates", json.dumps(file_event, default=str))
        logger.info("openai_json_schema_task completed successfully")
        return {"done": True}
    except Exception as e:
        logger.error(f"Error in openai_json_schema_task: {str(e)}", exc_info=True)
        raise 
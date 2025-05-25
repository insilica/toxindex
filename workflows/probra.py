import redis
import json
import time
import os
import uuid
from workflows.celery_worker import celery
from webserver.model.message import MessageSchema
from webserver import S3FileStorage

@celery.task(bind=True)
def probra_task(self, payload):
    """Example background task that emits progress messages and uploads a file."""
    r = redis.Redis()
    task_id = payload.get("task_id")
    sid = payload.get("sid")
    user_id = payload.get("user_id")

    for i in range(5):
        message = MessageSchema(role="assistant", content=f"Step {i}")
        event = {
            "type": "task_message",
            "data": message.model_dump(),
            "task_id": task_id,
            "sid": sid,
        }
        r.publish("celery_updates", json.dumps(event))
        time.sleep(1)

    # ------------------------------------------------------------------
    # Create a simple result file and upload it to S3
    tmp_filename = f"probra_result_{uuid.uuid4().hex}.txt"
    tmp_path = os.path.join("/tmp", tmp_filename)
    with open(tmp_path, "w", encoding="utf-8") as f:
        f.write("ProbRA task complete.\n")

    storage = S3FileStorage()
    s3_key = storage.upload_file(tmp_path, f"{task_id}/{tmp_filename}", "text/plain")
    download_url = storage.generate_download_url(s3_key)

    file_event = {
        "type": "task_file",
        "task_id": task_id,
        "sid": sid,
        "data": {
            "user_id": user_id,
            "filename": tmp_filename,
            "filepath": s3_key,
            "s3_url": download_url,
        },
    }
    r.publish("celery_updates", json.dumps(file_event))

    return {"done": True}

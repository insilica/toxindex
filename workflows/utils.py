import redis
import json
import os
import uuid
import logging
import tempfile

from pathlib import Path

# For data manipulation
import pandas as pd

# Webserver models and storage
from webserver.model.message import MessageSchema
from webserver.model.task import Task
from webserver.model.file import File
from webserver.storage import GCSFileStorage

# <----- Utility functions (for workflows) ----->
def emit_status(task_id, status):
    Task.set_status(task_id, status)
    r = redis.Redis(
        host=os.environ.get("REDIS_HOST", "localhost"),
        port=int(os.environ.get("REDIS_PORT", "6379"))
    )
    task = Task.get_task(task_id)

    # If no task exists, publish a minimal event
    if not task:
        event = {
            "type": "task_status_update",
            "task_id": task_id,
            "data": {"status": status},
        }
        r.publish("celery_updates", json.dumps(event, default=str))
        return

    event = {
        "type": "task_status_update",
        "task_id": task.task_id,
        "data": task.to_dict(),
    }
    r.publish("celery_updates", json.dumps(event, default=str))

def download_gcs_file_to_temp(gcs_path: str, temp_dir: Path) -> Path:
    """Download a file from GCS to a temporary local path with caching."""
    from webserver.cache_manager import cache_manager
    
    # Try to get from cache first
    cached_content = cache_manager.get_file_content(gcs_path)
    if cached_content:
        # Write cached content to temp file
        local_path = temp_dir / f"{uuid.uuid4().hex}_{Path(gcs_path).name}"
        with open(local_path, 'w', encoding='utf-8') as f:
            f.write(cached_content)
        return local_path
    
    # Cache miss - download from GCS
    gcs_storage = GCSFileStorage()
    local_path = temp_dir / f"{uuid.uuid4().hex}_{Path(gcs_path).name}"
    gcs_storage.download_file(gcs_path, str(local_path))
    
    # Cache the content for future use
    with open(local_path, 'r', encoding='utf-8') as f:
        content = f.read()
    cache_manager.cache_file_content(gcs_path, content)
    
    return local_path

def upload_local_file_to_gcs(local_path: Path, gcs_path: str) -> str:
    """Upload a local file to GCS and return the GCS path."""
    gcs_storage = GCSFileStorage()
    gcs_storage.upload_file(str(local_path), gcs_path)
    return gcs_path

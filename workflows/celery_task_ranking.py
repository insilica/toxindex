import redis
import json
import os
import uuid
import logging
import tempfile
import subprocess
import sys
from workflows.celery_worker_ranking import celery
from webserver.model.message import MessageSchema
from webserver.model.task import Task
from webserver.storage import GCSFileStorage
from webserver.model.file import File
from pathlib import Path

logging.getLogger().info("celery_task_ranking.py module loaded")
logger = logging.getLogger(__name__)


def get_redis_connection():
    """Get Redis connection with consistent configuration"""
    return redis.Redis(
        host=os.environ.get("REDIS_HOST", "localhost"),
        port=int(os.environ.get("REDIS_PORT", "6379"))
    )


def publish_to_celery_updates(event_type, task_id, data):
    """Publish event to celery_updates channel for database processing"""
    r = get_redis_connection()
    event = {
        "type": event_type,
        "task_id": task_id,
        "data": data,
    }
    r.publish("celery_updates", json.dumps(event, default=str))
    logger.info(f"Published {event_type} to celery_updates for task {task_id}")


def publish_to_socketio(event_name, room, data):
    """Publish event to Socket.IO Redis channel for real-time updates"""
    r = get_redis_connection()
    socketio_event = {
        "method": "emit",
        "event": event_name,
        "room": room,
        "data": data
    }
    r.publish("socketio", json.dumps(socketio_event, default=str))
    logger.info(f"Published {event_name} to Socket.IO for room {room}")


def emit_status(task_id, status):
    """Emit task status update to both database and real-time channels"""
    logger.info(f"[emit_status] {task_id} -> {status}")
    
    # Update database directly
    Task.set_status(task_id, status)
    task = Task.get_task(task_id)
    
    # Publish to celery_updates for any additional database processing
    publish_to_celery_updates("task_status_update", task.task_id, task.to_dict())
    
    # Publish to Socket.IO for real-time updates
    publish_to_socketio("task_status_update", f"task_{task_id}", task.to_dict())


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

@celery.task(bind=True, queue='ranking')
def ranking_task(self, payload):
    """Run the ranking workflow using ranking_workflow/src/main.py and upload results to GCS."""
    logger.info(f"=== TASK STARTED: Ranking Workflow ===")
    logger.info(f"Task ID: {self.request.id}")
    logger.info(f"Payload: {payload}")

    try:
        logger.info(f"Starting Ranking Workflow with payload: {payload}")
        r = get_redis_connection()
        task_id = payload.get("task_id")
        user_id = payload.get("user_id")
        file_id = payload.get("payload")

        if not all([task_id, user_id, file_id]):
            raise ValueError(f"Missing required fields. task_id={task_id}, user_id={user_id}, file_id={file_id}")

        emit_status(task_id, "fetching file from GCS")
        # Get file info from DB
        file_obj = File.get_file(file_id)
        if not file_obj or not file_obj.filepath:
            raise FileNotFoundError(f"Input file not found for file_id={file_id}")

        # Create temporary directory for processing
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)

            # Download input file from GCS
            input_path = download_gcs_file_to_temp(file_obj.filepath, temp_path)
            logger.info(f"Downloaded GCS file {file_obj.filepath} to {input_path}")

            logger.info(f"Processing task {task_id} for user {user_id}")
            emit_status(task_id, "starting")

            # Prepare to run ranking workflow
            score_type = (payload.get("score_type") or "PubChem").strip()
            endpoint = (payload.get("endpoint") or "endocrine disruption").strip()
            emit_status(task_id, f"running workflow ({score_type} / {endpoint})")
            main_py = (Path(__file__).resolve().parent.parent / 'ranking-workflow' / 'ranking_workflow' / 'src' / 'main.py').resolve()
            if not main_py.exists():
                raise FileNotFoundError(f"Ranking workflow entrypoint not found at {main_py}")

            output_dir = temp_path / 'output'
            output_dir.mkdir(parents=True, exist_ok=True)

            # main.py expects: input_file, output_dir, score_type, endpoint
            cmd = [sys.executable, str(main_py), str(input_path), str(output_dir), score_type, endpoint]
            logger.info(f"Running command: {' '.join(cmd)}")
            proc = subprocess.run(cmd, capture_output=True, text=True)
            logger.info(f"Workflow stdout:\n{proc.stdout}")
            if proc.returncode != 0:
                logger.error(f"Workflow stderr:\n{proc.stderr}")
                raise RuntimeError(f"Ranking workflow failed with exit code {proc.returncode}")

            # Find generated output file
            emit_status(task_id, "collecting results")
            ranked_files = list(output_dir.glob('ranked_*.txt'))
            if not ranked_files:
                raise FileNotFoundError("No ranked_*.txt file produced by workflow")
            ranked_file_path = sorted(ranked_files)[0]

            # Optionally send a brief message
            try:
                head_lines = []
                with open(ranked_file_path, 'r', encoding='utf-8') as f:
                    for _ in range(10):
                        line = f.readline()
                        if not line:
                            break
                        head_lines.append(line.rstrip('\n'))
                preview = "\n".join(head_lines)
                message = MessageSchema(
                    role="assistant",
                    content=f"Ranking completed. Preview of results (first lines):\n\n" + "\n".join([f"`{l}`" for l in head_lines])
                )
                emit_task_message(task_id, message.model_dump())
            except Exception as _:
                # If preview fails, proceed without blocking the task
                pass

            # Upload results to GCS
            emit_status(task_id, "uploading files to GCS")
            gcs_storage = GCSFileStorage()
            ranked_filename = ranked_file_path.name
            ranked_gcs_path = f"tasks/{task_id}/{ranked_filename}"
            gcs_storage.upload_file(str(ranked_file_path), ranked_gcs_path, content_type='text/plain')

            emit_status(task_id, "files uploaded")

            # Emit file event
            txt_file_data = {
                "user_id": user_id,
                "filename": ranked_filename,
                "filepath": ranked_gcs_path,
                "file_type": "txt",
                "content_type": "text/plain"
            }
            emit_task_file(task_id, txt_file_data)

        logger.info(f"Task {task_id} completed successfully")
        finished_at = Task.mark_finished(task_id)
        emit_status(task_id, "done")
        return {"done": True, "finished_at": finished_at}

    except Exception as e:
        logger.error(f"Error in ranking_task: {str(e)}", exc_info=True)
        emit_status(payload.get("task_id"), "error")
        raise
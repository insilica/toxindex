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
from workflows.utils import (
    get_redis_connection,
    publish_to_celery_updates,
    publish_to_socketio,
    emit_status,
    download_gcs_file_to_temp,
)

logging.getLogger().info("celery_task_ranking.py module loaded")
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
        file_id = payload.get("file_id")

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

            # Prepare to run ranking workflow take them from query
            score_type = (payload.get("score_type") or "GHS").strip()
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
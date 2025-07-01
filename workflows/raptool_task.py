import redis
import json
import os
import uuid
import logging
# import openai
# import hashlib
from workflows.celery_worker import celery
from webserver.model.message import MessageSchema
# from webserver.storage import S3FileStorage
# from RAP.toxicity_schema import TOXICITY_SCHEMA
# Update finished_at in the database
# import datetime
# import webserver.datastore as ds
from webserver.model.task import Task
import pandas as pd
from webserver.model.file import File
from raptool.parse_chemicals import parse_chemicals

logger = logging.getLogger(__name__)

@celery.task(bind=True)
def raptool_task(self, payload):
    """Run RAPtool parse_chemicals on an uploaded file and publish the result."""
    try:
        logger.info(f"[RAPtool] Running raptool_task for task_id={payload.get('task_id')}, user_id={payload.get('user_id')}, payload={payload}")
        r = redis.Redis()
        task_id = payload.get("task_id")
        user_id = payload.get("user_id")
        file_id = payload.get("file_id")
        if not all([task_id, user_id, file_id]):
            raise ValueError(f"Missing required fields. task_id={task_id}, user_id={user_id}, file_id={file_id}")

        # Get file info from DB
        file_obj = File.get_file(file_id)
        if not file_obj or not file_obj.filepath or not os.path.exists(file_obj.filepath):
            raise FileNotFoundError(f"Input file not found for file_id={file_id}")
        input_path = file_obj.filepath
        ext = os.path.splitext(input_path)[1].lower()

        # Prepare TXT if needed
        if ext == ".csv":
            txt_path = os.path.join("outputs", f"{os.path.splitext(os.path.basename(input_path))[0]}_{uuid.uuid4().hex}.txt")
            os.makedirs(os.path.dirname(txt_path), exist_ok=True)
            csv_to_txt_robust(input_path, txt_path)
            parse_input_path = txt_path
        else:
            parse_input_path = input_path

        # Prepare output CSV path
        output_filename = f"raptool_result_{uuid.uuid4().hex}.csv"
        output_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'outputs'))
        os.makedirs(output_dir, exist_ok=True)
        output_path = os.path.join(output_dir, output_filename)

        # Run RAPtool parse_chemicals
        parse_chemicals(parse_input_path, output_path)

        # Publish message event (simple notification)
        message = MessageSchema(role="assistant", content=f"RAPtool parse_chemicals completed. Output saved to {output_filename}.")
        event = {
            "type": "task_message",
            "data": message.model_dump(),
            "task_id": task_id,
        }
        r.publish("celery_updates", json.dumps(event, default=str))

        # Publish file event (for download)
        file_event = {
            "type": "task_file",
            "task_id": task_id,
            "data": {
                "user_id": user_id,
                "filename": output_filename,
                "filepath": output_path,
            },
        }
        r.publish("celery_updates", json.dumps(file_event, default=str))
        logger.info("raptool_task completed successfully")

        finished_at = Task.mark_finished(task_id)
        return {"done": True, "finished_at": finished_at}
    except Exception as e:
        logger.error(f"Error in raptool_task: {str(e)}", exc_info=True)
        raise

def csv_to_txt_robust(csv_path, txt_path):
    """
    Convert a CSV file to a plain text file with one chemical name per line.
    Prefer columns named 'name' or 'chemical', else use the first column.
    Handles missing files, empty files, encoding issues, and missing columns.
    """
    try:
        df = pd.read_csv(csv_path, encoding='utf-8')
    except UnicodeDecodeError:
        try:
            df = pd.read_csv(csv_path, encoding='latin1')
        except Exception as e:
            logging.error(f"Failed to read CSV file {csv_path} with utf-8 and latin1: {e}")
            raise
    except Exception as e:
        logging.error(f"Failed to read CSV file {csv_path}: {e}")
        raise

    if df.empty:
        logging.error(f"CSV file {csv_path} is empty.")
        raise ValueError(f"CSV file {csv_path} is empty.")

    # Prefer 'name' or 'chemical' columns
    col = None
    for candidate in ['name', 'chemical', 'Name', 'Chemical']:
        if candidate in df.columns:
            col = candidate
            break
    if col is None:
        col = df.columns[0]  # Use first column
        logging.warning(f"No 'name' or 'chemical' column found in {csv_path}. Using first column: {col}")

    names = df[col].dropna().astype(str).str.strip().unique()
    names = [n for n in names if n]
    if not names:
        logging.error(f"No valid chemical names found in column '{col}' of {csv_path}.")
        raise ValueError(f"No valid chemical names found in column '{col}' of {csv_path}.")

    with open(txt_path, 'w', encoding='utf-8') as f:
        for name in names:
            f.write(f"{name}\n")
    logging.info(f"Converted CSV {csv_path} to TXT {txt_path} with {len(names)} names.") 
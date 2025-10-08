import redis
import json
import os
import uuid
import logging
import tempfile
from pathlib import Path
from workflows.celery_app import celery
from webserver.model.message import MessageSchema
from webserver.model.task import Task
import pandas as pd
from webserver.model.file import File
from webserver.storage import GCSFileStorage
from RAPtool.parse_chemicals import parse_chemicals
from RAPtool.categorize_chemicals import categorize_chemicals
from RAPtool.predict_chemicals import predict_chemicals
from RAPtool.select_feature import select_feature
from RAPtool.build_heatmap import build_heatmap
from RAPtool.build_stripchart import build_stripchart
from workflows.utils import emit_status, download_gcs_file_to_temp, upload_local_file_to_gcs, publish_to_socketio

logger = logging.getLogger(__name__)


def emit_task_message(task_id, message_data):
    """Emit task message to both database and real-time channels"""
    r = redis.Redis(
        host=os.environ.get("REDIS_HOST", "localhost"),
        port=int(os.environ.get("REDIS_PORT", "6379"))
    )
    
    # Publish to Redis for database processing
    event = {
        "type": "task_message",
        "data": message_data,
        "task_id": task_id,
    }
    r.publish("celery_updates", json.dumps(event, default=str))
    
    # Emit to Socket.IO chat room
    task = Task.get_task(task_id)
    if task and getattr(task, 'session_id', None):
        room = f"chat_session_{task.session_id}"
        logger.info(f"[raptool_task] Emitting to chat room: {room}")
        try:
            publish_to_socketio("new_message", room, message_data)
            logger.info(f"[raptool_task] Successfully emitted to socket")
        except Exception as e:
            logger.error(f"[raptool_task] Failed to emit to socket: {e}")
    else:
        logger.warning(f"[raptool_task] No session_id found for task {task_id}")


@celery.task(bind=True, queue='raptool')
def raptool_task(self, payload):
    """Run RAPtool parse_chemicals on a GCS file and publish results back to GCS."""
    try:
        logger.info(f"[RAPtool GCS] Running raptool_task for task_id={payload.get('task_id')}, user_id={payload.get('user_id')}, payload={payload}")
        r = redis.Redis(
            host=os.environ.get("REDIS_HOST", "localhost"),
            port=int(os.environ.get("REDIS_PORT", "6379"))
        )
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
            
            ext = input_path.suffix.lower()

            # Prepare TXT if needed
            if ext == ".csv":
                emit_status(task_id, "converting csv")
                txt_path = temp_path / f"{input_path.stem}_{uuid.uuid4().hex}.txt"
                csv_to_txt_robust(input_path, txt_path)
                parse_input_path = txt_path
            else:
                parse_input_path = input_path

            emit_status(task_id, "parsing chemicals")
            # Prepare output CSV path in temp
            output_filename = "parsed_chemicals.csv"
            parsed_output_path = temp_path / output_filename
            
            # Run RAPtool parse_chemicals
            parse_chemicals(parse_input_path, parsed_output_path)

            # Upload parsed results to GCS
            gcs_parsed_path = f"tasks/{task_id}/parsed_chemicals.csv"
            upload_local_file_to_gcs(parsed_output_path, gcs_parsed_path)

            emit_status(task_id, "publishing parse result")
            # Publish message event
            message = MessageSchema(role="assistant", content=f"Chemicals parsed")
            emit_task_message(task_id, message.model_dump())

            # Publish file event (for download)
            file_event = {
                "type": "task_file",
                "task_id": task_id,
                "data": {
                    "user_id": user_id,
                    "filename": output_filename,
                    "filepath": gcs_parsed_path,  # GCS path instead of local
                },
            }
            r.publish("celery_updates", json.dumps(file_event, default=str))

            emit_status(task_id, "categorizing chemicals")
            # Run RAPtool categorize_chemicals on the parsed output
            categorize_chemicals(parsed_output_path, temp_path)
            categorized_file = temp_path / 'classified_chemicals.csv'
            
            # Upload categorized results to GCS
            gcs_categorized_path = f"tasks/{task_id}/classified_chemicals.csv"
            upload_local_file_to_gcs(categorized_file, gcs_categorized_path)
        
            emit_status(task_id, "publishing categorize result")
            # Publish message event
            message = MessageSchema(role="assistant", content=f"Chemicals categorized")
            emit_task_message(task_id, message.model_dump())

            # Publish categorized file event (for download)
            categorized_file_event = {
                "type": "task_file",
                "task_id": task_id,
                "data": {
                    "user_id": user_id,
                    "filename": categorized_file.name,
                    "filepath": gcs_categorized_path,  # GCS path instead of local
                },
            }
            r.publish("celery_updates", json.dumps(categorized_file_event, default=str))

            # Run RAPtool predict_chemicals on the categorized output
            predictions_file = temp_path / 'predicted_chemicals.parquet'
            emit_status(task_id, "predicting chemicals")
            predict_chemicals(categorized_file, predictions_file)
            
            # Upload predictions to GCS
            gcs_predictions_path = f"tasks/{task_id}/predicted_chemicals.parquet"
            upload_local_file_to_gcs(predictions_file, gcs_predictions_path)
            
            # Publish message event
            message = MessageSchema(role="assistant", content=f"ToxTransformer predictions generated")
            emit_task_message(task_id, message.model_dump())

            # Publish predictions file event (for download)
            predictions_file_event = {
                "type": "task_file",
                "task_id": task_id,
                "data": {
                    "user_id": user_id,
                    "filename": predictions_file.name,
                    "filepath": gcs_predictions_path,  # GCS path instead of local
                },
            }
            r.publish("celery_updates", json.dumps(predictions_file_event, default=str))

            # Run select_feature on the predictions output
            emit_status(task_id, "selecting features")
            
            feature_output_dir = temp_path / 'selected_properties'
            feature_output_dir.mkdir(exist_ok=True)
            select_feature(predictions_file, feature_output_dir)
            methods = ["lasso", "random_forest", "mutual_info", "rfe"]
            
            # Upload feature selection results to GCS
            for method in methods:
                csv_file = feature_output_dir / f"{method}_selected_properties.csv"
                if csv_file.exists():
                    gcs_csv_path = f"tasks/{task_id}/selected_properties/{method}_selected_properties.csv"
                    upload_local_file_to_gcs(csv_file, gcs_csv_path)
                    
                    feature_file_event = {
                        "type": "task_file",
                        "task_id": task_id,
                        "data": {
                            "user_id": user_id,
                            "filename": csv_file.name,
                            "filepath": gcs_csv_path,
                        },
                    }
                    r.publish("celery_updates", json.dumps(feature_file_event, default=str))
                
                txt_file = feature_output_dir / f"{method}_selected_properties.txt"
                if txt_file.exists():
                    gcs_txt_path = f"tasks/{task_id}/selected_properties/{method}_selected_properties.txt"
                    upload_local_file_to_gcs(txt_file, gcs_txt_path)
                    
                    feature_txt_file_event = {
                        "type": "task_file",
                        "task_id": task_id,
                        "data": {
                            "user_id": user_id,
                            "filename": txt_file.name,
                            "filepath": gcs_txt_path,
                        },
                    }
                    r.publish("celery_updates", json.dumps(feature_txt_file_event, default=str))
            
            # Publish message event for select_feature
            message = MessageSchema(
                role="assistant",
                content=f"features selected (methods used: lasso, random_forest, mutual_info, rfe)."
            )
            emit_task_message(task_id, message.model_dump())

            # Run build_heatmap with mutual_info-selected features
            emit_status(task_id, "building heatmap")
            
            heatmap_output = temp_path / 'heatmap_mutual_info.png'
            build_heatmap(predictions_file, heatmap_output, feature_selection_method='mutual_info')
            
            # Upload heatmap to GCS
            gcs_heatmap_path = f"tasks/{task_id}/heatmap_mutual_info.png"
            upload_local_file_to_gcs(heatmap_output, gcs_heatmap_path)
            
            # Publish message event for build_heatmap
            message = MessageSchema(role="assistant", content=f"heatmap generated")
            emit_task_message(task_id, message.model_dump())
            
            # Publish file event for heatmap image
            heatmap_file_event = {
                "type": "task_file",
                "task_id": task_id,
                "data": {
                    "user_id": user_id,
                    "filename": heatmap_output.name,
                    "filepath": gcs_heatmap_path,
                },
            }
            r.publish("celery_updates", json.dumps(heatmap_file_event, default=str))

            # Run build_stripchart with mutual_info-selected features
            emit_status(task_id, "building stripchart")
            
            stripchart_output = temp_path / 'stripchart_mutual_info.png'
            build_stripchart(predictions_file, stripchart_output, agg_func='median', feature_selection_method='mutual_info')
            
            # Upload stripchart to GCS
            gcs_stripchart_path = f"tasks/{task_id}/stripchart_mutual_info.png"
            upload_local_file_to_gcs(stripchart_output, gcs_stripchart_path)
            
            # Publish message event for build_stripchart
            message = MessageSchema(role="assistant", content=f"stripchart generated")
            emit_task_message(task_id, message.model_dump())
            
            # Publish file event for stripchart image
            stripchart_file_event = {
                "type": "task_file",
                "task_id": task_id,
                "data": {
                    "user_id": user_id,
                    "filename": stripchart_output.name,
                    "filepath": gcs_stripchart_path,
                },
            }
            r.publish("celery_updates", json.dumps(stripchart_file_event, default=str))

        logger.info("raptool_task completed successfully")
        finished_at = Task.mark_finished(task_id)
        emit_status(task_id, "done")
        return {"done": True, "finished_at": finished_at}
    except Exception as e:
        logger.error(f"Error in raptool_task: {str(e)}", exc_info=True)
        emit_status(task_id, "error")
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
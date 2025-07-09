import redis
import json
import os
import uuid
import logging
from workflows.celery_worker import celery
from webserver.model.message import MessageSchema
from webserver.model.task import Task
import pandas as pd
from webserver.model.file import File
from raptool.parse_chemicals import parse_chemicals
from raptool.categorize_chemicals import categorize_chemicals
from raptool.predict_chemicals import predict_chemicals
from raptool.select_feature import select_feature
from raptool.build_heatmap import build_heatmap
from raptool.build_stripchart import build_stripchart
import pathlib

logger = logging.getLogger(__name__)

# Always resolve outputs relative to project root
PROJECT_ROOT = pathlib.Path(__file__).resolve().parent.parent
OUTPUTS_ROOT = PROJECT_ROOT / "outputs"

def emit_status(task_id, status):
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
def raptool_task(self, payload):
    """Run RAPtool parse_chemicals on an uploaded file and publish the result."""
    try:
        logger.info(f"[RAPtool] Running raptool_task for task_id={payload.get('task_id')}, user_id={payload.get('user_id')}, payload={payload}")
        r = redis.Redis()
        task_id = payload.get("task_id")
        user_id = payload.get("user_id")
        file_id = payload.get("payload")
        if not all([task_id, user_id, file_id]):
            raise ValueError(f"Missing required fields. task_id={task_id}, user_id={user_id}, file_id={file_id}")

        emit_status(task_id, "fetching file")
        # Get file info from DB
        file_obj = File.get_file(file_id)
        if not file_obj or not file_obj.filepath or not os.path.exists(file_obj.filepath):
            raise FileNotFoundError(f"Input file not found for file_id={file_id}")
        input_path = file_obj.filepath
        ext = os.path.splitext(input_path)[1].lower()

        # Prepare TXT if needed
        if ext == ".csv":
            emit_status(task_id, "converting csv")
            txt_path = OUTPUTS_ROOT / f"{os.path.splitext(os.path.basename(input_path))[0]}_{uuid.uuid4().hex}.txt"
            txt_path.parent.mkdir(parents=True, exist_ok=True)
            csv_to_txt_robust(input_path, txt_path)
            parse_input_path = txt_path
        else:
            parse_input_path = input_path

        emit_status(task_id, "parsing chemicals")
        # Prepare unique output directory per task
        output_dir = OUTPUTS_ROOT / str(task_id)
        output_dir.mkdir(parents=True, exist_ok=True)
        # Prepare output CSV path
        output_filename = "parsed_chemicals.csv"
        parsed_output_path = output_dir / output_filename
        
        # Run RAPtool parse_chemicals
        parse_chemicals(parse_input_path, parsed_output_path)

        emit_status(task_id, "publishing parse result")
        # Publish message event (simple notification)
        message = MessageSchema(role="assistant", content=f"Chemicals parsed")
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
                "filepath": str(parsed_output_path),
            },
        }
        r.publish("celery_updates", json.dumps(file_event, default=str))

        emit_status(task_id, "categorizing chemicals")
        # Run RAPtool categorize_chemicals on the parsed output
        categorize_chemicals(parsed_output_path, output_dir)
        categorized_file = output_dir / 'classified_chemicals.csv'
    
        emit_status(task_id, "publishing categorize result")
        # Publish message event (simple notification)
        message = MessageSchema(role="assistant", content=f"Chemicals categorized")
        event = {
            "type": "task_message",
            "data": message.model_dump(),
            "task_id": task_id,
        }
        r.publish("celery_updates", json.dumps(event, default=str))

        # Publish categorized file event (for download)
        categorized_file_event = {
            "type": "task_file",
            "task_id": task_id,
            "data": {
                "user_id": user_id,
                "filename": categorized_file.name,
                "filepath": str(categorized_file),
            },
        }
        r.publish("celery_updates", json.dumps(categorized_file_event, default=str))

        # Run RAPtool predict_chemicals on the categorized output
        predictions_file = output_dir / 'predicted_chemicals.parquet'
        emit_status(task_id, "predicting chemicals")
        predict_chemicals(categorized_file, predictions_file)
        # Publish message event (simple notification)
        message = MessageSchema(role="assistant", content=f"ToxTransformer predictions generated")
        event = {
            "type": "task_message",
            "data": message.model_dump(),
            "task_id": task_id,
        }
        r.publish("celery_updates", json.dumps(event, default=str))

        # Publish predictions file event (for download)
        predictions_file_event = {
            "type": "task_file",
            "task_id": task_id,
            "data": {
                "user_id": user_id,
                "filename": predictions_file.name,
                "filepath": str(predictions_file),
            },
        }
        r.publish("celery_updates", json.dumps(predictions_file_event, default=str))

        # Run select_feature on the predictions output
        emit_status(task_id, "selecting features")
        
        feature_output_dir = output_dir / 'selected_properties' #note this is hardcoded in build_feature.py
        feature_output_dir.mkdir(exist_ok=True)
        select_feature(predictions_file, feature_output_dir)
        methods = ["lasso", "random_forest", "mutual_info", "rfe"]
        # Publish message event for select_feature
        message = MessageSchema(
            role="assistant",
            content=f"features selected (methods used: lasso, random_forest, mutual_info, rfe)."
        )
        event = {
            "type": "task_message",
            "data": message.model_dump(),
            "task_id": task_id,
        }
        r.publish("celery_updates", json.dumps(event, default=str))
        # Publish file events for each method's CSV and TXT files
        for method in methods:
            # Publish CSV file
            csv_file = feature_output_dir / f"{method}_selected_properties.csv"
            if csv_file.exists():
                feature_file_event = {
                    "type": "task_file",
                    "task_id": task_id,
                    "data": {
                        "user_id": user_id,
                        "filename": csv_file.name,
                        "filepath": str(csv_file),
                    },
                }
                r.publish("celery_updates", json.dumps(feature_file_event, default=str))
            
            # Publish TXT file
            txt_file = feature_output_dir / f"{method}_selected_properties.txt"
            if txt_file.exists():
                feature_txt_file_event = {
                    "type": "task_file",
                    "task_id": task_id,
                    "data": {
                        "user_id": user_id,
                        "filename": txt_file.name,
                        "filepath": str(txt_file),
                    },
                }
                r.publish("celery_updates", json.dumps(feature_txt_file_event, default=str))

        # Run build_heatmap with mutual_info-selected features
        emit_status(task_id, "building heatmap")
        
        heatmap_output = output_dir / 'heatmap_mutual_info.png'
        build_heatmap(predictions_file, heatmap_output, feature_selection_method='mutual_info')
        # Publish message event for build_heatmap
        message = MessageSchema(role="assistant", content=f"heatmap generated")
        event = {
            "type": "task_message",
            "data": message.model_dump(),
            "task_id": task_id,
        }
        r.publish("celery_updates", json.dumps(event, default=str))
        # Publish file event for heatmap image
        heatmap_file_event = {
            "type": "task_file",
            "task_id": task_id,
            "data": {
                "user_id": user_id,
                "filename": heatmap_output.name,
                "filepath": str(heatmap_output),
            },
        }
        r.publish("celery_updates", json.dumps(heatmap_file_event, default=str))

        # Run build_stripchart with mutual_info-selected features
        emit_status(task_id, "building stripchart")
        
        stripchart_output = output_dir / 'stripchart_mutual_info.png'
        build_stripchart(predictions_file, stripchart_output, agg_func='median', feature_selection_method='mutual_info')
        # Publish message event for build_stripchart
        message = MessageSchema(role="assistant", content=f"stripchart generated")
        event = {
            "type": "task_message",
            "data": message.model_dump(),
            "task_id": task_id,
        }
        r.publish("celery_updates", json.dumps(event, default=str))
        # Publish file event for stripchart image
        stripchart_file_event = {
            "type": "task_file",
            "task_id": task_id,
            "data": {
                "user_id": user_id,
                "filename": stripchart_output.name,
                "filepath": str(stripchart_output),
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
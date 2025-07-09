import redis
import json
import os
import uuid
import logging
import pathlib
from workflows.celery_worker import celery
from webserver.model.message import MessageSchema
from webserver.model.task import Task
from webserver.model.file import File
from pathway_analysis_tool.annotate_pathway import save_pathway
import pandas as pd

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
def pathway_analysis_task(self, payload):
    """Run pathway analysis using pathway_analysis_tool on uploaded data."""
    try:
        logger.info(f"[Pathway Analysis] Running pathway_analysis_task for task_id={payload.get('task_id')}, user_id={payload.get('user_id')}, payload={payload}")
        r = redis.Redis()
        task_id = payload.get("task_id")
        user_id = payload.get("user_id")
        pathway_id = payload.get("pathway_id", "WP3657")  # Default to a common pathway
        property_name = payload.get("property_name", "value")  # Default property column
        entity_kind = payload.get("entity_kind", "gene")  # Default to gene entities
        file_id = payload.get("file_id")
        
        if not all([task_id, user_id]):
            raise ValueError(f"Missing required fields. task_id={task_id}, user_id={user_id}")

        emit_status(task_id, "fetching file")
        
        # Get file info from DB if file_id provided, otherwise use default data
        if file_id:
            file_obj = File.get_file(file_id)
            if not file_obj or not file_obj.filepath or not os.path.exists(file_obj.filepath):
                raise FileNotFoundError(f"Input file not found for file_id={file_id}")
            input_path = file_obj.filepath
        else:
            # Create minimal default data if no file provided
            emit_status(task_id, "creating default data")
            output_dir = OUTPUTS_ROOT / str(task_id)
            output_dir.mkdir(parents=True, exist_ok=True)
            input_path = output_dir / "default_data.csv"
            
            # Create minimal sample data with entity and value columns
            sample_data = pd.DataFrame({
                'entity': ['GENE1', 'GENE2', 'GENE3'],
                property_name: [0.5, 0.8, 0.3]
            })
            sample_data.to_csv(input_path, index=False)

        emit_status(task_id, "preparing output directory")
        # Prepare unique output directory per task
        output_dir = OUTPUTS_ROOT / str(task_id)
        output_dir.mkdir(parents=True, exist_ok=True)

        emit_status(task_id, "running pathway analysis")
        # Run pathway analysis using the save_pathway function
        save_pathway(
            pathway=pathway_id,
            prop=property_name,
            kind=entity_kind,
            data=input_path,
            outdir=output_dir
        )

        emit_status(task_id, "publishing results")
        # Publish message event
        message = MessageSchema(role="assistant", content=f"Pathway analysis completed for {pathway_id}")
        event = {
            "type": "task_message",
            "data": message.model_dump(),
            "task_id": task_id,
        }
        r.publish("celery_updates", json.dumps(event, default=str))

        # Publish image file event
        image_path = output_dir / "images" / f"{pathway_id}_{property_name}.png"
        if image_path.exists():
            image_file_event = {
                "type": "task_file",
                "task_id": task_id,
                "data": {
                    "user_id": user_id,
                    "filename": image_path.name,
                    "filepath": str(image_path),
                },
            }
            r.publish("celery_updates", json.dumps(image_file_event, default=str))

        logger.info("pathway_analysis_task completed successfully")
        finished_at = Task.mark_finished(task_id)
        emit_status(task_id, "done")
        return {"done": True, "finished_at": finished_at}
        
    except Exception as e:
        logger.error(f"Error in pathway_analysis_task: {str(e)}", exc_info=True)
        emit_status(task_id, "error")
        raise 
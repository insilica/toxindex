import redis
import json
import os
import logging
import tempfile
from pathlib import Path
from workflows.celery_app import celery
from webserver.model.message import MessageSchema
from webserver.model.task import Task
from webserver.model.file import File
from pathway_analysis_tool import save_pathway
import pandas as pd
from webserver.storage import GCSFileStorage
from webserver.cache_manager import cache_manager
from workflows.utils import emit_status, download_gcs_file_to_temp, upload_local_file_to_gcs

logger = logging.getLogger(__name__)


@celery.task(bind=True, queue='pathway')
def pathway_analysis_task(self, payload):
    """GCS-enabled pathway analysis task that downloads input files from GCS and uploads outputs to GCS."""
    try:
        logger.info(f"[Pathway Analysis GCS] Running pathway_analysis_task for task_id={payload.get('task_id')}, user_id={payload.get('user_id')}, payload={payload}")
        r = redis.Redis(
            host=os.environ.get("REDIS_HOST", "localhost"),
            port=int(os.environ.get("REDIS_PORT", "6379"))
        )
        task_id = payload.get("task_id")
        user_id = payload.get("user_id")
        pathway_id = payload.get("pathway_id", "WP3657")  # Default to a common pathway
        property_name = payload.get("property_name", "value")  # Default property column
        entity_kind = payload.get("entity_kind", "gene")  # Default to gene entities
        file_id = payload.get("file_id")
        
        if not all([task_id, user_id]):
            raise ValueError(f"Missing required fields. task_id={task_id}, user_id={user_id}")

        emit_status(task_id, "fetching file")
        
        with tempfile.TemporaryDirectory() as temp_dir:
            temp_path = Path(temp_dir)
            
            # Get file info from DB if file_id provided, otherwise use default data
            if file_id:
                file_obj = File.get_file(file_id)
                if not file_obj or not file_obj.filepath:
                    raise FileNotFoundError(f"Input file not found for file_id={file_id}")
                
                # Download from GCS
                input_path = download_gcs_file_to_temp(file_obj.filepath, temp_path)
            else:
                # Create minimal default data if no file provided
                emit_status(task_id, "creating default data")
                input_path = temp_path / "default_data.csv"
                
                # Create minimal sample data with entity and value columns
                sample_data = pd.DataFrame({
                    'entity': ['GENE1', 'GENE2', 'GENE3'],
                    property_name: [0.5, 0.8, 0.3]
                })
                sample_data.to_csv(input_path, index=False)

            emit_status(task_id, "running pathway analysis")
            # Run pathway analysis using the save_pathway function
            save_pathway(
                pathway=pathway_id,
                prop=property_name,
                kind=entity_kind,
                data=input_path,
                outdir=temp_path
            )

            emit_status(task_id, "uploading results to GCS")
            # Upload generated files to GCS
            gcs_storage = GCSFileStorage()
            
            # Upload image file if it exists
            image_path = temp_path / "images" / f"{pathway_id}_{property_name}.png"
            if image_path.exists():
                gcs_image_path = f"tasks/{task_id}/pathway_{pathway_id}_{property_name}.png"
                upload_local_file_to_gcs(image_path, gcs_image_path, content_type='image/png')
                
                # Publish image file event
                image_file_event = {
                    "type": "task_file",
                    "task_id": task_id,
                    "data": {
                        "user_id": user_id,
                        "filename": image_path.name,
                        "filepath": gcs_image_path,
                    },
                }
                r.publish("celery_updates", json.dumps(image_file_event, default=str))
                logger.info(f"Uploaded pathway image to GCS: {gcs_image_path}")

            # Upload any other generated files
            for file_path in temp_path.rglob("*"):
                if file_path.is_file() and file_path.suffix in ['.csv', '.json', '.txt']:
                    gcs_file_path = f"tasks/{task_id}/pathway_{file_path.name}"
                    upload_local_file_to_gcs(file_path, gcs_file_path)
                    
                    # Publish file event
                    file_event = {
                        "type": "task_file",
                        "task_id": task_id,
                        "data": {
                            "user_id": user_id,
                            "filename": file_path.name,
                            "filepath": gcs_file_path,
                        },
                    }
                    r.publish("celery_updates", json.dumps(file_event, default=str))
                    logger.info(f"Uploaded pathway file to GCS: {gcs_file_path}")

        emit_status(task_id, "publishing results")
        # Publish message event
        message = MessageSchema(role="assistant", content=f"Pathway analysis completed for {pathway_id}")
        event = {
            "type": "task_message",
            "data": message.model_dump(),
            "task_id": task_id,
        }
        r.publish("celery_updates", json.dumps(event, default=str))

        logger.info("pathway_analysis_task completed successfully")
        finished_at = Task.mark_finished(task_id)
        emit_status(task_id, "done")
        return {"done": True, "finished_at": finished_at}
        
    except Exception as e:
        logger.error(f"Error in pathway_analysis_task: {str(e)}", exc_info=True)
        emit_status(task_id, "error")
        raise 
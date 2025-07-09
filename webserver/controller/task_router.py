"""
Dynamic task router that reads workflows from JSON to determine which Celery task to run.
This eliminates hardcoded if/elif statements and keeps task routing in sync with workflow definitions.
"""
import json
import logging
from pathlib import Path
from typing import Dict, Any, Optional

# Import all available Celery tasks
from workflows.probra import probra_task
from workflows.plain_openai_tasks import plain_openai_task, openai_json_schema_task
from workflows.raptool_task import raptool_task
from workflows.pathway_analysis_task import pathway_analysis_task

logger = logging.getLogger(__name__)

# Task mapping - maps celery_task names to actual task functions
TASK_MAPPING = {
    'probra_task': probra_task,
    'plain_openai_task': plain_openai_task,
    'openai_json_schema_task': openai_json_schema_task,
    'raptool_task': raptool_task,
    'pathway_analysis_task': pathway_analysis_task,
}

def get_workflows_file_path() -> Path:
    """Get the path to the workflows JSON file"""
    return Path(__file__).parent.parent.parent / "resources" / "default_workflows.json"

def load_workflows() -> Dict[int, Dict[str, Any]]:
    """Load workflows from JSON file"""
    try:
        workflows_file = get_workflows_file_path()
        if not workflows_file.exists():
            logger.error(f"Workflows file not found: {workflows_file}")
            return {}
            
        with open(workflows_file, 'r') as f:
            data = json.load(f)
        
        # Create a mapping by workflow_id for easy lookup
        workflows = {}
        for workflow in data['workflows']:
            workflows[workflow['workflow_id']] = workflow
        
        return workflows
    except Exception as e:
        logger.error(f"Failed to load workflows: {e}")
        return {}

def get_task_for_workflow(workflow_id: int) -> Optional[Any]:
    """Get the Celery task function for a given workflow_id"""
    workflows = load_workflows()
    workflow = workflows.get(workflow_id)
    
    if not workflow:
        logger.error(f"No workflow found for workflow_id: {workflow_id}")
        return None
    
    celery_task_name = workflow.get('celery_task')
    if not celery_task_name:
        logger.error(f"No celery_task specified for workflow_id: {workflow_id}")
        return None
    
    task_function = TASK_MAPPING.get(celery_task_name)
    if not task_function:
        logger.error(f"No task function found for celery_task: {celery_task_name}")
        return None
    
    return task_function

def create_task_payload(workflow_id: int, task_id: str, user_id: str, message: str = "", file_id: str = None) -> Dict[str, Any]:
    """Create the appropriate payload for a given workflow"""
    workflows = load_workflows()
    workflow = workflows.get(workflow_id)
    
    if not workflow:
        logger.error(f"No workflow found for workflow_id: {workflow_id}")
        return {}
    
    # Base payload
    payload = {
        "task_id": task_id,
        "user_id": user_id,
    }
    
    # Add workflow-specific payload logic
    celery_task_name = workflow.get('celery_task')
    
    if celery_task_name == 'probra_task':
        payload["payload"] = message
    elif celery_task_name == 'plain_openai_task':
        payload["payload"] = message
    elif celery_task_name == 'openai_json_schema_task':
        payload["payload"] = message
    elif celery_task_name == 'raptool_task':
        payload["payload"] = file_id
    elif celery_task_name == 'pathway_analysis_task':
        # Extract pathway parameters from message or use defaults
        payload.update({
            "pathway_id": "WP3657",  # Default pathway
            "property_name": "value",  # Default property column
            "entity_kind": "gene",  # Default entity type
            "file_id": file_id,
        })
        
        # If message contains pathway info, try to parse it
        if message and "WP" in message.upper():
            import re
            pathway_match = re.search(r'WP\d+', message.upper())
            if pathway_match:
                payload["pathway_id"] = pathway_match.group()
    
    return payload

def route_task(workflow_id: int, task_id: str, user_id: str, message: str = "", file_id: str = None):
    """Route a task to the appropriate Celery task based on workflow_id"""
    task_function = get_task_for_workflow(workflow_id)
    
    if not task_function:
        logger.error(f"Could not determine task function for workflow_id: {workflow_id}")
        return None
    
    payload = create_task_payload(workflow_id, task_id, user_id, message, file_id)
    
    if not payload:
        logger.error(f"Could not create payload for workflow_id: {workflow_id}")
        return None
    
    logger.info(f"Routing workflow_id {workflow_id} to {task_function.__name__} with payload: {payload}")
    
    # Execute the task
    return task_function.delay(payload) 
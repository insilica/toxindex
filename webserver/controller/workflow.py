from flask import Blueprint, jsonify
import flask_login
import json
import logging
from webserver.controller.task_router import get_workflows_file_path
from webserver.model.user_group import UserGroup

workflow_bp = Blueprint('workflows', __name__, url_prefix='/api/workflows')

@workflow_bp.route('/config', methods=['GET'])
@flask_login.login_required
def get_workflows_config():
    """Serve the workflows configuration from resources directory, filtered by user group"""
    try:
        workflows_file = get_workflows_file_path()
        if not workflows_file.exists():
            return jsonify({"error": "Workflows configuration not found"}), 404
        
        with open(workflows_file, 'r') as f:
            data = json.load(f)
        
        # Filter workflows based on user's group access
        user_id = flask_login.current_user.user_id
        accessible_workflows = UserGroup.get_accessible_workflows(user_id)
        accessible_workflow_ids = {w['workflow_id'] for w in accessible_workflows}
        
        # Filter the workflows list
        if 'workflows' in data:
            data['workflows'] = [
                w for w in data['workflows'] 
                if w.get('workflow_id') in accessible_workflow_ids
            ]
        
        return jsonify(data)
    except Exception as e:
        logging.error(f"Error reading workflows config: {e}")
        return jsonify({"error": "Failed to load workflows configuration"}), 500

@workflow_bp.route('/list', methods=['GET'])
@flask_login.login_required
def list_workflows():
    """List available workflows for the current user"""
    try:
        from webserver.model.workflow import Workflow
        user_id = flask_login.current_user.user_id
        workflows = Workflow.get_workflows_by_user(user_id)
        return jsonify({"workflows": [w.to_dict() for w in workflows]})
    except Exception as e:
        logging.error(f"Error listing workflows: {e}")
        return jsonify({"error": "Failed to list workflows"}), 500 
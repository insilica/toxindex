from flask import Blueprint, request, jsonify
import flask_login
import logging
from webserver.model.user_group import UserGroup
from webserver.util import is_valid_uuid
import webserver.datastore as ds
from webserver.csrf import csrf

admin_bp = Blueprint('admin', __name__, url_prefix='/api/admin')

def require_admin():
    """Decorator to require admin access"""
    def decorator(f):
        def wrapper(*args, **kwargs):
            if not flask_login.current_user.is_authenticated:
                return jsonify({"error": "Authentication required"}), 401
            
            user_group = UserGroup.get_user_group(flask_login.current_user.user_id)
            if not user_group or user_group.name != 'admin':
                return jsonify({"error": "Admin access required"}), 403
            
            return f(*args, **kwargs)
        wrapper.__name__ = f.__name__
        return wrapper
    return decorator

# ============================================================================
# USER MANAGEMENT ENDPOINTS
# ============================================================================

@admin_bp.route('/users', methods=['GET'])
@flask_login.login_required
@require_admin()
def get_all_users():
    """Get all users with their group information (Admin only)"""
    try:
        users = UserGroup.get_all_users_with_groups()
        return jsonify({
            "users": [
                {
                    "user_id": str(user['user_id']),
                    "email": user['email'],
                    "group_name": user.get('group_name', 'No Group'),
                    "group_description": user.get('group_description'),
                    "created_at": user['created_at'].strftime('%Y-%m-%d %H:%M:%S') if user['created_at'] else None
                }
                for user in users
            ]
        })
    except Exception as e:
        logging.error(f"Error getting users: {e}")
        return jsonify({"error": "Failed to get users"}), 500

@admin_bp.route('/users/<user_id>/group', methods=['PUT'])
@flask_login.login_required
@require_admin()
@csrf.exempt
def update_user_group(user_id):
    """Update a user's group (Admin only)"""
    if not is_valid_uuid(user_id):
        return jsonify({"error": "Invalid user ID"}), 400
    
    try:
        data = request.get_json()
        group_id = data.get('group_id')
        
        if not group_id or not is_valid_uuid(group_id):
            return jsonify({"error": "Valid group_id required"}), 400
        
        # Verify the group exists
        group = UserGroup.get_group(group_id)
        if not group:
            return jsonify({"error": "Group not found"}), 404
        
        UserGroup.set_user_group(user_id, group_id)
        
        return jsonify({
            "success": True,
            "message": f"User {user_id} assigned to group {group.name}"
        })
    except Exception as e:
        logging.error(f"Error updating user group: {e}")
        return jsonify({"error": "Failed to update user group"}), 500

@admin_bp.route('/groups/<group_id>/users', methods=['GET'])
@flask_login.login_required
@require_admin()
def get_users_in_group(group_id):
    """Get all users in a specific group (Admin only)"""
    if not is_valid_uuid(group_id):
        return jsonify({"error": "Invalid group ID"}), 400
    
    try:
        users = UserGroup.get_users_in_group(group_id)
        return jsonify({
            "users": [
                {
                    "user_id": str(user['user_id']),
                    "email": user['email'],
                    "group_name": user.get('group_name'),
                    "created_at": user['created_at'].strftime('%Y-%m-%d %H:%M:%S') if user['created_at'] else None
                }
                for user in users
            ]
        })
    except Exception as e:
        logging.error(f"Error getting users in group: {e}")
        return jsonify({"error": "Failed to get users in group"}), 500

# ============================================================================
# GROUP MANAGEMENT ENDPOINTS
# ============================================================================

@admin_bp.route('/groups', methods=['GET'])
@flask_login.login_required
@require_admin()
def get_all_groups():
    """Get all user groups (Admin only)"""
    try:
        groups = UserGroup.get_all_groups()
        return jsonify({
            "groups": [group.to_dict() for group in groups]
        })
    except Exception as e:
        logging.error(f"Error getting groups: {e}")
        return jsonify({"error": "Failed to get groups"}), 500

@admin_bp.route('/groups', methods=['POST'])
@flask_login.login_required
@require_admin()
@csrf.exempt
def create_group():
    """Create a new user group (Admin only)"""
    try:
        data = request.get_json()
        name = data.get('name')
        description = data.get('description')
        
        if not name:
            return jsonify({"error": "Group name required"}), 400
        
        # Check if group already exists
        existing = UserGroup.get_group_by_name(name)
        if existing:
            return jsonify({"error": "Group with this name already exists"}), 409
        
        group = UserGroup.create_group(name, description)
        if group:
            return jsonify({
                "success": True,
                "group": group.to_dict()
            })
        else:
            return jsonify({"error": "Failed to create group"}), 500
    except Exception as e:
        logging.error(f"Error creating group: {e}")
        return jsonify({"error": "Failed to create group"}), 500

@admin_bp.route('/groups/<group_id>', methods=['PUT'])
@flask_login.login_required
@require_admin()
@csrf.exempt
def update_group(group_id):
    """Update a user group (Admin only)"""
    if not is_valid_uuid(group_id):
        return jsonify({"error": "Invalid group ID"}), 400
    
    try:
        data = request.get_json()
        name = data.get('name')
        description = data.get('description')
        
        group = UserGroup.update_group(group_id, name, description)
        if group:
            return jsonify({
                "success": True,
                "group": group.to_dict()
            })
        else:
            return jsonify({"error": "Group not found"}), 404
    except Exception as e:
        logging.error(f"Error updating group: {e}")
        return jsonify({"error": "Failed to update group"}), 500

@admin_bp.route('/groups/<group_id>', methods=['DELETE'])
@flask_login.login_required
@require_admin()
@csrf.exempt
def delete_group(group_id):
    """Delete a user group (Admin only)"""
    if not is_valid_uuid(group_id):
        return jsonify({"error": "Invalid group ID"}), 400
    
    try:
        group = UserGroup.get_group(group_id)
        if not group:
            return jsonify({"error": "Group not found"}), 404
        
        # Prevent deletion of admin group
        if group.name == 'admin':
            return jsonify({"error": "Cannot delete admin group"}), 403
        
        UserGroup.delete_group(group_id)
        return jsonify({
            "success": True,
            "message": f"Group {group.name} deleted"
        })
    except Exception as e:
        logging.error(f"Error deleting group: {e}")
        return jsonify({"error": "Failed to delete group"}), 500

# ============================================================================
# WORKFLOW ACCESS MANAGEMENT ENDPOINTS
# ============================================================================

@admin_bp.route('/workflow-access', methods=['GET'])
@flask_login.login_required
@require_admin()
def get_workflow_access():
    """Get workflow access status for all groups and workflows (Admin only)"""
    try:
        # Get all workflows
        from webserver.model.workflow import Workflow
        workflows = Workflow.get_workflows_by_user(None)  # Get all workflows
        
        # Get all groups
        groups = UserGroup.get_all_groups()
        
        # Get all workflow access records
        access_records = ds.find_all("""
            SELECT group_id, workflow_id FROM workflow_group_access
        """)
        
        # Create a set of existing access records for quick lookup
        existing_access = {(str(record['group_id']), record['workflow_id']) for record in access_records}
        
        # Build the response
        access_matrix = []
        for group in groups:
            for workflow in workflows:
                access_matrix.append({
                    'group_id': str(group.group_id),
                    'group_name': group.name,
                    'workflow_id': workflow.workflow_id,
                    'workflow_title': workflow.title,
                    'has_access': (str(group.group_id), workflow.workflow_id) in existing_access
                })
        
        return jsonify({
            "access_matrix": access_matrix,
            "groups": [group.to_dict() for group in groups],
            "workflows": [workflow.to_dict() for workflow in workflows]
        })
    except Exception as e:
        logging.error(f"Error getting workflow access: {e}")
        return jsonify({"error": "Failed to get workflow access"}), 500

@admin_bp.route('/groups/<group_id>/workflows/<int:workflow_id>/access', methods=['POST'])
@flask_login.login_required
@require_admin()
@csrf.exempt
def grant_workflow_access(group_id, workflow_id):
    """Grant a group access to a workflow (Admin only)"""
    if not is_valid_uuid(group_id):
        return jsonify({"error": "Invalid group ID"}), 400
    
    try:
        UserGroup.grant_workflow_access(group_id, workflow_id)
        return jsonify({
            "success": True,
            "message": f"Access granted to workflow {workflow_id}"
        })
    except Exception as e:
        logging.error(f"Error granting workflow access: {e}")
        return jsonify({"error": "Failed to grant workflow access"}), 500

@admin_bp.route('/groups/<group_id>/workflows/<int:workflow_id>/access', methods=['DELETE'])
@flask_login.login_required
@require_admin()
@csrf.exempt
def revoke_workflow_access(group_id, workflow_id):
    """Revoke a group's access to a workflow (Admin only)"""
    if not is_valid_uuid(group_id):
        return jsonify({"error": "Invalid group ID"}), 400
    
    try:
        UserGroup.revoke_workflow_access(group_id, workflow_id)
        return jsonify({
            "success": True,
            "message": f"Access revoked from workflow {workflow_id}"
        })
    except Exception as e:
        logging.error(f"Error revoking workflow access: {e}")
        return jsonify({"error": "Failed to revoke workflow access"}), 500 
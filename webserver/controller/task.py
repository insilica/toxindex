from flask import Blueprint, request, jsonify, send_file, abort
import flask_login
from webserver.model import Task
from webserver.model import ChatSession
from webserver.csrf import csrf
import logging, datetime
from webserver.util import is_valid_uuid
from webserver import datastore as ds
from webserver.ai_service import generate_title
from webserver.controller.task_router import route_task
from webserver.model.file import File
import os

task_bp = Blueprint('tasks', __name__, url_prefix='/api/tasks')

@task_bp.route('', methods=['GET'])
@flask_login.login_required
def get_user_tasks():
    try:
        user_id = flask_login.current_user.user_id
        env_id = request.args.get("environment_id")
        if env_id:
            tasks = Task.get_tasks_by_environment(env_id, user_id)
        else:
            tasks = Task.get_tasks_by_user(user_id)
        active_tasks = [t for t in tasks if not t.archived]
        archived_tasks = [t for t in tasks if t.archived]
        def sort_key(t):
            return t.last_accessed or t.created_at or datetime.datetime.min
        active_tasks = sorted(active_tasks, key=sort_key, reverse=True)
        archived_tasks = sorted(archived_tasks, key=sort_key, reverse=True)
        return jsonify({
            "active_tasks": [t.to_dict() for t in active_tasks],
            "archived_tasks": [t.to_dict() for t in archived_tasks],
        })
    except Exception as e:
        logging.error(f"Error retrieving tasks: {str(e)}")
        return jsonify({"error": "Failed to retrieve tasks"}), 500

@csrf.exempt
@task_bp.route('', methods=['POST'])
@flask_login.login_required
def create_task():
    task_data = request.get_json()
    message = task_data.get("message", "")
    title = generate_title(message)
    user_id = flask_login.current_user.user_id
    workflow_id = int(task_data.get("workflow", 1))
    environment_id = task_data.get("environment_id")
    sid = task_data.get("sid")
    file_id = task_data.get("file_id")
    created_at = datetime.datetime.now(datetime.timezone.utc)

    # If no session id, create a new session
    if not sid:
        session = ChatSession.create_session(environment_id, user_id, title=title)
        sid = session.session_id if session else None

    logging.info(f"Controller: about to create task with created_at={created_at} (type: {type(created_at)})")
    task = Task.create_task(
        title=title,
        user_id=user_id,
        workflow_id=workflow_id,
        environment_id=environment_id,
        session_id=sid,
        created_at=created_at,
    )
    Task.add_message(task.task_id, flask_login.current_user.user_id, "user", message, session_id=sid)
    celery_task = route_task(workflow_id, task.task_id, str(user_id), message, file_id)
    Task.update_celery_task_id(task.task_id, celery_task.id)
    return jsonify({
        "task_id": task.task_id,
        "celery_id": celery_task.id,
        "session_id": sid,
        "created_at": task.created_at.strftime("%Y-%m-%d %H:%M:%S") if task.created_at else None,
        "finished_at": task.finished_at.strftime("%Y-%m-%d %H:%M:%S") if hasattr(task, 'finished_at') and task.finished_at else None
    })

@csrf.exempt
@task_bp.route('/<task_id>', methods=['GET'])
@flask_login.login_required
def get_task(task_id):
    try:
        if not is_valid_uuid(task_id):
            return jsonify({"error": "Invalid task ID"}), 400
        task = Task.get_task(task_id)
        if not task or task.user_id != flask_login.current_user.user_id:
            return jsonify({"error": "Task not found"}), 404
        return jsonify(task.to_dict())
    except Exception as e:
        logging.error(f"/api/tasks/{task_id} error: {e}", exc_info=True)
        return jsonify({"error": "Internal server error"}), 500

@csrf.exempt
@task_bp.route('/<task_id>/archive', methods=['POST'])
@flask_login.login_required
def archive_task(task_id):
    if not is_valid_uuid(task_id):
        return jsonify({"success": False, "error": "Invalid task ID"}), 400
    user_id = flask_login.current_user.user_id
    ds.execute("UPDATE tasks SET archived = TRUE WHERE task_id = %s AND user_id = %s", (task_id, user_id))
    return jsonify({"success": True})

@csrf.exempt
@task_bp.route('/<task_id>/unarchive', methods=['POST'])
@flask_login.login_required
def unarchive_task(task_id):
    if not is_valid_uuid(task_id):
        return jsonify({"success": False, "error": "Invalid task ID"}), 400
    user_id = flask_login.current_user.user_id
    ds.execute("UPDATE tasks SET archived = FALSE WHERE task_id = %s AND user_id = %s", (task_id, user_id))
    return jsonify({"success": True})

@task_bp.route('/<task_id>/files', methods=['GET'])
@flask_login.login_required
def get_task_files(task_id):
    files = File.get_files_by_task(task_id)
    return jsonify({'files': [f.to_dict() for f in files]})

@task_bp.route('/<task_id>/messages', methods=['GET'])
@flask_login.login_required
def get_task_messages(task_id):
    user_id = flask_login.current_user.user_id
    messages = Task.get_messages(task_id, user_id)
    return jsonify({'messages': messages})

@task_bp.route('/<task_id>/files/<file_id>/download', methods=['GET'])
@flask_login.login_required
def download_task_file(task_id, file_id):
    if not is_valid_uuid(task_id) or not is_valid_uuid(file_id):
        return jsonify({'error': 'Invalid task or file ID'}), 400
    file = File.get_file(file_id)
    if not file or not file.filepath or not os.path.exists(file.filepath) or str(file.task_id) != str(task_id):
        return abort(404)
    return send_file(file.filepath, as_attachment=True, download_name=file.filename) 
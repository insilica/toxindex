import eventlet
eventlet.monkey_patch()

# Standard library imports
import os
import uuid
import base64
import mimetypes
import threading
import logging
import json
import re
from datetime import datetime

# Third-party imports
import dotenv
from werkzeug.utils import secure_filename
import flask
import flask_login
from flask_socketio import SocketIO, emit, join_room
from flask_wtf.csrf import CSRFError
from flask_wtf import CSRFProtect
from celery.result import AsyncResult
from flask import request, send_from_directory, jsonify, send_file, abort

# Local imports
import webserver.datastore as ds
from webserver import login_manager as LM
from webserver.controller import login
from webserver.model import Task, Workflow, Message, File, Environment, ChatSession, User
from webserver.ai_service import generate_title
from workflows.plain_openai_tasks import plain_openai_task, openai_json_schema_task
from workflows.probra import probra_task
from workflows.celery_worker import celery
from workflows.chat_response_task import chat_response_task

print("RUNNING WEBSERVER/APP.PY")

dotenv.load_dotenv()

import redis

# FLASK APP ===================================================================
static_folder_path = os.path.join(os.path.dirname(__file__), "webserver", "static")
app = flask.Flask(__name__, static_folder=static_folder_path)
# app.config["SERVER_NAME"] = os.environ.get("SERVER_NAME")
app.config["PREFERRED_URL_SCHEME"] = os.environ.get("PREFERRED_URL_SCHEME")
app.config["TEMPLATES_AUTO_RELOAD"] = True
app.secret_key = os.environ.get("FLASK_APP_SECRET_KEY")

logging.info(f"SECRET_KEY: {app.secret_key}")
logging.info(f"WTF_CSRF_ENABLED: {app.config.get('WTF_CSRF_ENABLED', 'not set')}")

socketio = SocketIO(
    app,
    cors_allowed_origins="*",
    message_queue="redis://localhost:6379/0",
    manage_session=False,
    # async_mode='gevent',
    logger=True,
    engineio_logger=True
)

# Generate a unique log filename with date and time
os.makedirs('logs', exist_ok=True)
log_filename = f"logs/app_{datetime.now().strftime('%Y-%m-%d_%H')}.log"

# Configure logging to output to file with detailed formatting
logging.basicConfig(
    level=logging.DEBUG,
    format='%(asctime)s - %(name)s - %(levelname)s - %(filename)s:%(lineno)d - %(message)s',
    handlers=[
        logging.FileHandler(log_filename),
        logging.StreamHandler()  # Also keep console output
    ]
)

# Set Flask app logger to use same configuration
app.logger.setLevel(logging.INFO)
for handler in logging.getLogger().handlers:
    app.logger.addHandler(handler)
LM.init(app)
Workflow.load_default_workflows()

csrf = CSRFProtect(app)

logging.info(f"DB ENV (Flask startup): PGHOST={os.getenv('PGHOST')}, PGPORT={os.getenv('PGPORT')}, PGDATABASE={os.getenv('PGDATABASE')}, PGUSER={os.getenv('PGUSER')}, PGPASSWORD={os.getenv('PGPASSWORD')}")

# REDIS PUB/SUB BACKGROUND LISTENER ==========================================
def redis_listener(name):
    import uuid
    listener_id = uuid.uuid4().hex[:8]
    r = redis.Redis()
    pubsub = r.pubsub()
    pubsub.subscribe("celery_updates")
    logging.info(f"[redis_listener] ({listener_id}) Started Redis listener thread: {name}")
    for redis_message in pubsub.listen():
        logging.debug(f"[redis_listener] ({listener_id}) Raw redis_message: {redis_message}")
        if redis_message["type"] != "message":
            continue
        try:
            raw_data = redis_message["data"]
            if isinstance(raw_data, bytes):
                raw_data = raw_data.decode()
            event = json.loads(raw_data)
            logging.info(f"[redis_listener] ({listener_id}) Redis event: {event}")
            event_type = event.get("type")
            event_task_id = event.get("task_id")
            event_data = event.get("data")
            if event_task_id is None or event_data is None:
                logging.warning(f"[redis_listener] ({listener_id}) Redis event missing required fields: event_type={event_type}, task_id={event_task_id}, data={event_data}")
                continue
            if event_type not in ("task_message", "task_file"):
                logging.debug(f"[redis_listener] ({listener_id}) Ignoring event_type={event_type}")
                continue
            task = Task.get_task(event_task_id)
            if not task:
                logging.warning(f"[redis_listener] ({listener_id}) No task found for task_id={event_task_id}")
                continue
            if event_type == "task_message":
                logging.info(f"[redis_listener] ({listener_id}) Processing task_message for task_id={event_task_id}")
                Message.process_event(task, event_data)
            elif event_type == "task_file":
                logging.info(f"[redis_listener] ({listener_id}) Processing task_file for task_id={event_task_id}")
                File.process_event(task, event_data)
            room = f"task_{task.task_id}"
            logging.info(f"[redis_listener] ({listener_id}) Emitting {event_type} to room {room} with data: {event_data}")
            socketio.emit(event_type, event_data, to=room)
            logging.info(f"[redis_listener] ({listener_id}) {event_type} sent to {room}")
        except json.JSONDecodeError:
            logging.error(f"[redis_listener] ({listener_id}) Failed to decode Redis message as JSON.")
        except Exception as e:
            logging.error(f"[redis_listener] ({listener_id}) Redis listener error: {e}", exc_info=True)


# INDEX / ROOT ===============================================================
# @app.route("/", methods=["GET"])
# def index():
#     ...


# SOCKETIO HANDLERS ==========================================================
@socketio.on("connect")
def handle_connect(auth):
    logging.info(f"[socketio] Client connected: {request.sid}")
    emit("connected", {"sid": request.sid})


# TODO right now I think any user can join any task room
@socketio.on("join_task_room")
def handle_join_task_room(data):
    task_id = data.get("task_id")
    room = f"task_{task_id}"
    logging.info(f"[socketio] {request.sid} joining room: {room}")
    join_room(room)
    logging.info(f"[socketio] {request.sid} joined room: {room}")
    emit("joined_task_room", {"task_id": task_id})


# TASK MANAGEMENT ============================================================
@app.route("/task/new", methods=["POST"])
def task_create():
    logging.info(f"[task_create] Called with {request.method}")
    task_data = request.get_json()
    message = task_data.get("message", "")
    logging.info(f"[task_create] Received message: {message}")
    title = generate_title(message)
    user_id = flask_login.current_user.user_id
    workflow_id = int(task_data.get("workflow", 1))
    environment_id = task_data.get("environment_id")
    sid = task_data.get("sid")
    logging.info(f"[task_create] Creating task for user_id={user_id}, workflow_id={workflow_id}, environment_id={environment_id}, sid={sid}")
    task = Task.create_task(
        title=title,
        user_id=user_id,
        workflow_id=workflow_id,
        environment_id=environment_id,
    )
    logging.info(f"[task_create] Created task: {task.to_dict() if task else None}")
    Task.add_message(task.task_id, flask_login.current_user.user_id, "user", message)
    
    # Select the appropriate task based on workflow_id
    if workflow_id == 1:
        logging.info(f"[task_create] Executing ToxRAP (probra_task) for task_id={task.task_id}")
        celery_task = probra_task.delay(
            {
                "payload": message,
                "sid": sid,
                "task_id": task.task_id,
                "user_id": str(user_id),
            }
        )
    elif workflow_id == 2:
        logging.info(f"[task_create] Executing ToxDirect (plain_openai_task) for task_id={task.task_id}")
        celery_task = plain_openai_task.delay(
            {
                "payload": message,
                "sid": sid,
                "task_id": task.task_id,
                "user_id": str(user_id),
            }
        )
    elif workflow_id == 3:
        logging.info(f"[task_create] Executing ToxJson (openai_json_schema_task) for task_id={task.task_id}")
        celery_task = openai_json_schema_task.delay(
            {
                "payload": message,
                "sid": sid,
                "task_id": task.task_id,
                "user_id": str(user_id),
            }
        )
    else:
        logging.info(f"[task_create] Executing default (probra_task) for task_id={task.task_id}")
        celery_task = probra_task.delay(  # Default to probra_task for unknown workflows
            {
                "payload": message,
                "sid": sid,
                "task_id": task.task_id,
                "user_id": str(user_id),
            }
        )
    
    Task.update_celery_task_id(task.task_id, celery_task.id)

    return flask.jsonify({"task_id": task.task_id, "celery_id": celery_task.id})


@app.route("/task/<task_id>/status", methods=["GET"])
def task_status(task_id):
    task = Task.get_task(task_id, flask_login.current_user.user_id)
    celery_id = task.celery_task_id
    if not celery_id:
        return flask.jsonify({"error": "Missing celery_id"}), 400

    result = AsyncResult(celery_id, app=celery)
    return flask.jsonify(
        {
            "state": result.state,
            "info": result.info if result.info else {},
            "ready": result.ready(),
            "result": result.result if result.ready() else None,
        }
    )


@app.route("/task/<int:task_id>/archive", methods=["POST"])
@flask_login.login_required
def archive_task(task_id):
    user_id = flask_login.current_user.user_id
    ds.execute("UPDATE tasks SET archived = TRUE WHERE task_id = %s AND user_id = %s", (task_id, user_id))
    return flask.jsonify({"success": True})


@app.route("/task/<int:task_id>/unarchive", methods=["POST"])
@flask_login.login_required
def unarchive_task(task_id):
    user_id = flask_login.current_user.user_id
    ds.execute("UPDATE tasks SET archived = FALSE WHERE task_id = %s AND user_id = %s", (task_id, user_id))
    return flask.jsonify({"success": True})


@app.route("/tasks", methods=["GET"])
@flask_login.login_required
def get_user_tasks():
    try:
        user_id = flask_login.current_user.user_id
        env_id = request.args.get("environment_id")
        if env_id:
            tasks = Task.get_tasks_by_environment(env_id, user_id)
        else:
            tasks = Task.get_tasks_by_user(user_id)
        # Split into active and archived, sort by last_accessed or created_at
        active_tasks = [t for t in tasks if not t.archived]
        archived_tasks = [t for t in tasks if t.archived]
        def sort_key(t):
            return t.last_accessed or t.created_at or datetime.datetime.min
        active_tasks = sorted(active_tasks, key=sort_key, reverse=True)
        archived_tasks = sorted(archived_tasks, key=sort_key, reverse=True)
        return flask.jsonify({
            "active_tasks": [t.to_dict() for t in active_tasks],
            "archived_tasks": [t.to_dict() for t in archived_tasks],
        })
    except Exception as e:
        logging.error(f"Error retrieving tasks: {str(e)}")
        return flask.jsonify({"error": "Failed to retrieve tasks"}), 500


# @app.route("/task/<task_id>", methods=["GET"])
# def task(task_id):
#     messages = Message.get_messages(task_id)
#     files = File.get_files(task_id)

#     # Find the latest assistant message
#     assistant_messages = [m for m in messages if m.role == "assistant"]
#     latest_assistant_message = assistant_messages[-1].content if assistant_messages else None

#     # Pre-render markdown content for convenience
#     rendered_files = []
#     for f in files:
#         file_info = f.to_dict()
#         if f.filename.lower().endswith((".md", ".markdown")):
#             file_info["html"] = f.markdown_to_html()
#         else:
#             file_info["html"] = ""
#         rendered_files.append(file_info)

#     # Render latest assistant message as HTML (if present)
#     if latest_assistant_message:
#         import markdown
#         latest_assistant_message_html = markdown.markdown(latest_assistant_message, extensions=["extra", "tables"])
#     else:
#         latest_assistant_message_html = None

#     return flask.render_template(
#         "task.html", 
#         task_id=task_id, 
#         messages=messages, 
#         files=rendered_files, 
#         latest_assistant_message=latest_assistant_message,
#         latest_assistant_message_html=latest_assistant_message_html
#     )


# MESSAGE MANAGEMENT =========================================================
@app.route("/message/new", methods=["POST"])
def message_new():
    message_data = request.get_json()
    task_id = message_data.get('task_id')
    content = message_data.get('content')
    sid = message_data.get('sid')
    
    if not all([task_id, content]):
        return flask.jsonify({"error": "Missing required fields"}), 400
        
    # Get the task to check its workflow
    task = Task.get_task(task_id)
    if not task:
        return flask.jsonify({"error": "Task not found"}), 404
        
    # Store the user message
    Message.create_message(task_id, None, "user", content)
    
    # For interactive echo task (workflow_id 3), publish to Redis channel
    if task.workflow_id == 3:
        r = redis.Redis()
        input_channel = f"task:{task_id}:input"
        r.publish(input_channel, json.dumps({"content": content}))
        logging.info(f"Published message to {input_channel}: {content}")
    
    return flask.jsonify({"message_id": message_data.get("message_id")})


# WORKFLOW MANAGEMENT ========================================================
# @app.route("/workflow/new", methods=["GET", "POST"])
# def workflow_new():
#     if request.method == "POST":
#         workflow_data = request.get_json()
#         title = workflow_data.get("title", "New workflow")
#         description = workflow_data.get("description", "")
#         workflow = Workflow.create_workflow(
#             title=title,
#             description=description,
#             user_id=flask_login.current_user.user_id,
#         )
#         logging.info(f"created workflow {workflow.workflow_id}")
#         return flask.jsonify({"workflow_id": workflow.workflow_id})
#     return flask.render_template("create_workflow.html")


# @app.route("/workflow/<workflow_id>", methods=["GET"])
# def workflow(workflow_id):
#     workflow = Workflow.get_workflow(workflow_id)
#     return flask.render_template("workflow_run.html", workflow=workflow)


# ENVIRONMENT MANAGEMENT =====================================================
@csrf.exempt
@app.route("/environment/new", methods=["POST"])
@flask_login.login_required
def environment_new():
    env_data = flask.request.get_json()
    title = env_data.get("title", "New Environment")
    description = env_data.get("description")
    env = Environment.create_environment(
        title, flask_login.current_user.user_id, description
    )
    return flask.jsonify({"environment_id": env.environment_id})

@csrf.exempt
@app.route("/api/environment/new", methods=["POST"])
@flask_login.login_required
def api_environment_new():
    env_data = flask.request.get_json()
    title = env_data.get("title", "New Environment")
    description = env_data.get("description")
    env = Environment.create_environment(
        title, flask_login.current_user.user_id, description
    )
    return flask.jsonify({"environment_id": env.environment_id})

@app.route("/environments", methods=["GET"])
@flask_login.login_required
def environment_list():
    user_id = flask_login.current_user.user_id
    envs = Environment.get_environments_by_user(user_id)
    logging.info(f"[environment_list] user_id={user_id}, envs={[e.to_dict() for e in envs]}")
    return flask.jsonify({"environments": [e.to_dict() for e in envs]})


# @app.route("/environment/<env_id>", methods=["GET"])
# @flask_login.login_required
# def environment_view(env_id):
#     workflow = Workflow.get_workflow(1)
#     tasks = Task.get_tasks_by_environment(env_id, flask_login.current_user.user_id)
#     env = Environment.get_environment(env_id)
#     user_id = flask_login.current_user.user_id
#     user_workflows = Workflow.get_workflows_by_user(user_id)
#     global_workflows = Workflow.get_workflows_by_user(None)
#     all_workflows = {str(w.workflow_id): w.title for w in (user_workflows + global_workflows)}
#     logging.info(f"[environment_view] workflow_titles mapping: {all_workflows}")
#     return flask.render_template(
#         "environment.html", tasks=tasks, workflow=workflow, environment=env, workflow_titles=all_workflows
#     )


@app.route("/environment/<env_id>", methods=["DELETE"])
@flask_login.login_required
def environment_delete(env_id):
    if not is_valid_uuid(env_id):
        logging.error(f"Attempted to delete environment with invalid UUID: {env_id}")
        return flask.jsonify({"success": False, "error": "Invalid environment ID"}), 400
    try:
        # Disk cleanup: remove files from disk before deleting environment
        files = File.get_files_by_environment(env_id)
        for file in files:
            try:
                if file.filepath and os.path.exists(file.filepath):
                    os.remove(file.filepath)
                    logging.info(f"Deleted file from disk: {file.filepath}")
            except Exception as e:
                logging.warning(f"Failed to remove file from disk: {file.filepath}, error: {e}")
        # Now delete the environment (DB rows for files/tasks will be deleted by ON DELETE CASCADE)
        Environment.delete_environment(env_id, flask_login.current_user.user_id)
        return flask.jsonify({"success": True})
    except Exception as e:
        logging.error(f"Failed to delete environment {env_id}: {e}")
        return flask.jsonify({"success": False, "error": str(e)}), 500


# LOGIN MANAGEMENT ===========================================================
app.route("/register", methods=["GET", "POST"])(login.register)
app.route("/api/register", methods=["GET", "POST"])(csrf.exempt(login.register))
app.route("/verify", methods=["GET"])(login.verify_message)
app.route("/api/verify", methods=["GET"])(csrf.exempt(login.verify_message))
app.route("/verification/<token>", methods=["GET", "POST"])(login.verification)
app.route("/api/verification/<token>", methods=["GET", "POST"])(csrf.exempt(login.verification))
app.route("/login", methods=["GET", "POST"])(login.login)
app.route("/api/login", methods=["GET", "POST"])(csrf.exempt(login.login))
app.route("/logout", methods=["GET"])(login.logout)
app.route("/api/logout", methods=["GET"])(csrf.exempt(login.logout))
app.route("/forgot_password", methods=["GET", "POST"])(login.forgot_password)
app.route("/api/forgot_password", methods=["GET", "POST"])(csrf.exempt(login.forgot_password))
app.route("/reset_password/<token>", methods=["GET", "POST"])(login.reset_password)
app.route("/api/reset_password/<token>", methods=["GET", "POST"])(csrf.exempt(login.reset_password))


# ICONS ======================================================================
@app.route("/favicon.ico")
def favicon():
    return app.send_static_file("favicon.png")


@app.route("/icons/<path:filename>")
def serve_icon(filename):
    icons_dir = os.path.join(os.path.dirname(__file__), "icons")
    return send_from_directory(icons_dir, filename)

@app.route('/log_tab_switch', methods=['POST'])
def log_tab_switch():
    data = flask.request.get_json()
    tab = data.get('tab')
    task_id = data.get('task_id')
    timestamp = data.get('timestamp')
    logging.info(f"[Tab Switch] User switched to tab '{tab}' for task_id={task_id} at {timestamp}")
    return flask.jsonify({'status': 'ok'})

@app.route("/api/me", methods=["GET"])
def api_me():
    import flask_login
    from flask import make_response, jsonify
    logging.warning(f"current_user: {flask_login.current_user}")
    logging.warning(f"is_authenticated: {flask_login.current_user.is_authenticated}")
    if not flask_login.current_user.is_authenticated:
        response = make_response(jsonify({"error": "Not authenticated"}), 401)
    else:
        response = make_response(jsonify({
            "user_id": flask_login.current_user.user_id,
            "email": flask_login.current_user.email,
        }))
    # Prevent caching
    response.headers["Cache-Control"] = "no-store, no-cache, must-revalidate, max-age=0"
    response.headers["Pragma"] = "no-cache"
    response.headers["Expires"] = "0"
    return response

@app.route("/test-alive")
def test_alive():
    return "ALIVE"

@app.route("/api/environments", methods=["GET"])
@flask_login.login_required
def api_environments():
    try:
        user_id = flask_login.current_user.user_id
        envs = Environment.get_environments_by_user(user_id)
        if not envs:
            Environment.create_environment("Base environment", user_id, description="Your starter workspace")
            envs = Environment.get_environments_by_user(user_id)
        return flask.jsonify({
            "environments": [
                {
                    "environment_id": e.environment_id,
                    "title": e.title,
                    "user_id": flask_login.current_user.email,  # Use email for creator
                    "created_at": e.created_at.strftime("%Y-%m-%dT%H:%M:%S") if e.created_at else ""
                }
                for e in envs
            ]
        })
    except Exception as e:
        logging.error(f"/api/environments error: {e}", exc_info=True)
        return flask.jsonify({"error": "Internal server error"}), 500

# FILE UPLOAD ENDPOINT ======================================================
@csrf.exempt
@app.route('/api/upload-file', methods=['POST'])
@flask_login.login_required
def upload_file():
    if 'file' not in flask.request.files:
        return flask.jsonify({'error': 'No file part'}), 400
    file = flask.request.files['file']
    if file.filename == '':
        return flask.jsonify({'error': 'No selected file'}), 400
    if not file.filename.lower().endswith('.csv'):
        return flask.jsonify({'error': 'Only CSV files are allowed'}), 400
    filename = secure_filename(file.filename)
    upload_dir = os.path.join(os.path.dirname(__file__), '..', 'uploads')
    os.makedirs(upload_dir, exist_ok=True)
    file_path = os.path.join(upload_dir, filename)
    file.save(file_path)
    # Store file metadata in the database
    environment_id = flask.request.form.get('environment_id')
    user_id = flask_login.current_user.user_id
    File.create_file(
        task_id=None,
        user_id=user_id,
        filename=filename,
        filepath=file_path,
        s3_url='',
        environment_id=environment_id
    )
    return flask.jsonify({'success': True, 'filename': filename})

# List files for an environment
@app.route('/api/environment/<env_id>/files', methods=['GET'])
@flask_login.login_required
def environment_files(env_id):
    files = File.get_files_by_environment(env_id)
    return flask.jsonify({'files': [f.to_dict() for f in files]})

@csrf.exempt
@app.route('/api/file/<file_id>', methods=['DELETE'])
@flask_login.login_required
def delete_file(file_id):
    try:
        # Find the file
        file = File.get_file(file_id)
        if not file:
            return flask.jsonify({'success': False, 'error': 'File not found'}), 404
        # Remove from disk if present
        try:
            if file.filepath and os.path.exists(file.filepath):
                os.remove(file.filepath)
        except Exception as e:
            logging.warning(f"Failed to remove file from disk: {e}")
        # Remove from database
        File.delete_file(file_id, flask_login.current_user.user_id)
        return flask.jsonify({'success': True})
    except Exception as e:
        logging.error(f"Failed to delete file {file_id}: {e}")
        return flask.jsonify({'success': False, 'error': str(e)}), 500

@app.route('/api/file/<file_id>/inspect', methods=['GET'])
@flask_login.login_required
def inspect_file(file_id):
    import logging
    logging.info(f"[inspect_file] Called with file_id={file_id} (type: {type(file_id)})")
    file = File.get_file(file_id)
    logging.info(f"[inspect_file] File.get_file({file_id}) returned: {file}")
    if not file:
        logging.warning(f"[inspect_file] No file found in DB for file_id={file_id}")
        return flask.jsonify({'error': 'File not found'}), 404
    if not file.filepath:
        logging.warning(f"[inspect_file] File object for file_id={file_id} has no filepath")
        return flask.jsonify({'error': 'File not found'}), 404
    if not os.path.exists(file.filepath):
        logging.warning(f"[inspect_file] Filepath does not exist on disk: {file.filepath}")
        return flask.jsonify({'error': 'File not found'}), 404
    ext = os.path.splitext(file.filename)[1].lower()
    mimetype, _ = mimetypes.guess_type(file.filename)
    try:
        if ext in ['.txt', '.csv', '.json', '.md', '.markdown']:
            with open(file.filepath, 'r', encoding='utf-8') as f:
                content = f.read()
            if ext == '.json':
                import json
                try:
                    parsed = json.loads(content)
                    return flask.jsonify({'type': 'json', 'content': parsed, 'filename': file.filename, 'mimetype': mimetype})
                except Exception as e:
                    return flask.jsonify({'type': 'text', 'content': content, 'filename': file.filename, 'mimetype': mimetype, 'warning': f'Invalid JSON: {e}'})
            elif ext in ['.md', '.markdown']:
                import markdown
                html = markdown.markdown(content, extensions=["extra", "tables"])
                return flask.jsonify({'type': 'markdown', 'content': html, 'filename': file.filename, 'mimetype': mimetype})
            elif ext == '.csv':
                import csv
                import io
                reader = csv.reader(io.StringIO(content))
                rows = list(reader)
                return flask.jsonify({'type': 'csv', 'content': rows, 'filename': file.filename, 'mimetype': mimetype})
            else:
                return flask.jsonify({'type': 'text', 'content': content, 'filename': file.filename, 'mimetype': mimetype})
        elif ext in ['.xlsx', '.xls']:
            try:
                import pandas as pd
                df = pd.read_excel(file.filepath)
                preview = df.head(100).to_dict(orient='records')
                columns = list(df.columns)
                return flask.jsonify({'type': 'xlsx', 'content': preview, 'columns': columns, 'filename': file.filename, 'mimetype': mimetype})
            except Exception as e:
                return flask.jsonify({'error': f'Failed to parse Excel: {e}', 'type': 'error', 'filename': file.filename, 'mimetype': mimetype}), 400
        elif ext in ['.png', '.jpg', '.jpeg', '.gif']:
            with open(file.filepath, 'rb') as f:
                data = f.read()
            b64 = base64.b64encode(data).decode('utf-8')
            data_url = f"data:{mimetype};base64,{b64}"
            return flask.jsonify({'type': 'image', 'content': data_url, 'filename': file.filename, 'mimetype': mimetype})
        else:
            return flask.jsonify({'error': 'Preview not supported for this file type', 'type': 'unsupported', 'filename': file.filename, 'mimetype': mimetype}), 415
    except Exception as e:
        logging.error(f"[inspect_file] Exception for file_id={file_id}: {e}")
        return flask.jsonify({'error': str(e), 'type': 'error', 'filename': file.filename, 'mimetype': mimetype}), 500

@app.route('/api/file/<file_id>/download', methods=['GET'])
@flask_login.login_required
def download_file(file_id):
    file = File.get_file(file_id)
    print(f"DOWNLOAD: file_id={file_id}, file={file}, filepath={getattr(file, 'filepath', None)}")
    if not file or not file.filepath or not os.path.exists(file.filepath):
        print("DOWNLOAD: File not found or missing on disk")
        return flask.abort(404)
    return flask.send_file(file.filepath, as_attachment=True, download_name=file.filename)

# API endpoint to trigger probra_task from dashboard chatbox
@csrf.exempt
@app.route("/api/run-probra-task", methods=["POST"])
@flask_login.login_required
def api_run_probra_task():
    data = flask.request.get_json()
    prompt = data.get("prompt")
    environment_id = data.get("environment_id")
    user_id = flask_login.current_user.user_id
    if not prompt:
        return flask.jsonify({"error": "Missing prompt"}), 400
    # Always create a new chat session for dashboard submissions
    session = ChatSession.create_session(environment_id, user_id, title=prompt[:60])
    session_id = session.session_id if session else None
    # Create a new Task (workflow_id=1 for ToxIndex RAP)
    title = generate_title(prompt)
    task = Task.create_task(
        title=title,
        user_id=user_id,
        workflow_id=1,
        environment_id=environment_id,
        session_id=session_id,
    )
    Task.add_message(task.task_id, user_id, "user", prompt, session_id=session_id)
    celery_task = probra_task.delay({
        "payload": prompt,
        "task_id": task.task_id,
        "user_id": str(user_id),
    })
    Task.update_celery_task_id(task.task_id, celery_task.id)
    return flask.jsonify({"task_id": task.task_id, "celery_id": celery_task.id, "session_id": session_id})

@app.route("/api/tasks", methods=["GET"])
@flask_login.login_required
def api_get_user_tasks():
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
        return flask.jsonify({
            "active_tasks": [t.to_dict() for t in active_tasks],
            "archived_tasks": [t.to_dict() for t in archived_tasks],
        })
    except Exception as e:
        logging.error(f"Error retrieving tasks: {str(e)}")
        return flask.jsonify({"error": "Failed to retrieve tasks"}), 500

@csrf.exempt
@app.route("/api/tasks/<task_id>", methods=["GET"])
@flask_login.login_required
def api_get_task(task_id):
    try:
        if not is_valid_uuid(task_id):
            return flask.jsonify({"error": "Invalid task ID"}), 400
        task = Task.get_task(task_id)
        if not task or task.user_id != flask_login.current_user.user_id:
            return flask.jsonify({"error": "Task not found"}), 404
        return flask.jsonify(task.to_dict())
    except Exception as e:
        logging.error(f"/api/tasks/{task_id} error: {e}", exc_info=True)
        return flask.jsonify({"error": "Internal server error"}), 500

@csrf.exempt
@app.route("/api/tasks/<task_id>/archive", methods=["POST"])
@flask_login.login_required
def api_archive_task(task_id):
    if not is_valid_uuid(task_id):
        return flask.jsonify({"success": False, "error": "Invalid task ID"}), 400
    user_id = flask_login.current_user.user_id
    ds.execute("UPDATE tasks SET archived = TRUE WHERE task_id = %s AND user_id = %s", (task_id, user_id))
    return flask.jsonify({"success": True})

@csrf.exempt
@app.route("/api/tasks/<task_id>/unarchive", methods=["POST"])
@flask_login.login_required
def api_unarchive_task(task_id):
    if not is_valid_uuid(task_id):
        return flask.jsonify({"success": False, "error": "Invalid task ID"}), 400
    user_id = flask_login.current_user.user_id
    ds.execute("UPDATE tasks SET archived = FALSE WHERE task_id = %s AND user_id = %s", (task_id, user_id))
    return flask.jsonify({"success": True})

@csrf.exempt
@app.route("/api/run-vanilla-task", methods=["POST"])
@flask_login.login_required
def api_run_vanilla_task():
    data = flask.request.get_json()
    prompt = data.get("prompt")
    environment_id = data.get("environment_id")
    user_id = flask_login.current_user.user_id
    if not prompt:
        return flask.jsonify({"error": "Missing prompt"}), 400
    # Always create a new chat session for dashboard submissions
    session = ChatSession.create_session(environment_id, user_id, title=prompt[:60])
    session_id = session.session_id if session else None
    # Create a new Task (workflow_id=2 for ToxIndex Vanilla)
    title = generate_title(prompt)
    task = Task.create_task(
        title=title,
        user_id=user_id,
        workflow_id=2,
        environment_id=environment_id,
        session_id=session_id,
    )
    Task.add_message(task.task_id, user_id, "user", prompt, session_id=session_id)
    logging.info(f"[api_run_vanilla_task] Queuing plain_openai_task for task_id={task.task_id}, user_id={user_id}")
    celery_task = plain_openai_task.delay({
        "payload": prompt,
        "task_id": task.task_id,
        "user_id": str(user_id),
    })
    logging.info(f"[api_run_vanilla_task] Queued plain_openai_task with celery_id={celery_task.id}")
    Task.update_celery_task_id(task.task_id, celery_task.id)
    return flask.jsonify({"task_id": task.task_id, "celery_id": celery_task.id, "session_id": session_id})

# --- Chat Session Endpoints ---


def is_valid_uuid(val):
    try:
        uuid.UUID(str(val))
        return True
    except Exception:
        return False

@csrf.exempt
@app.route('/api/environment/<env_id>/chat_sessions', methods=['POST'])
@flask_login.login_required
def create_chat_session(env_id):
    if env_id in ("__add__", "__manage__") or not is_valid_uuid(env_id):
        return jsonify({"error": "Invalid environment ID"}), 400
    user_id = flask_login.current_user.user_id
    d = flask.request.get_json(force=True)
    title = d.get('title') or 'New chat'
    session = ChatSession.create_session(env_id, user_id, title)
    if not session:
        return jsonify({"error": "Failed to create chat session"}), 500
    return jsonify(session.to_dict())

@app.route('/api/environment/<env_id>/chat_sessions', methods=['GET'])
@flask_login.login_required
def list_chat_sessions(env_id):
    user_id = flask_login.current_user.user_id
    sessions = ChatSession.get_sessions_by_environment(env_id, user_id)
    return jsonify({'sessions': [s.to_dict() for s in sessions]})

@app.route('/api/chat_session/<session_id>/messages', methods=['GET'])
@flask_login.login_required
def get_chat_session_messages(session_id):
    messages = Message.get_messages_by_session(session_id)
    return jsonify({'messages': [m.to_dict() for m in messages]})

@csrf.exempt
@app.route('/api/chat_session/<session_id>/message', methods=['POST'])
@flask_login.login_required
def send_message_in_session(session_id):
    from webserver.model import Task
    user_id = flask_login.current_user.user_id
    data = request.get_json() or {}
    prompt = data.get('prompt')
    model = data.get('model')
    environment_id = data.get('environment_id')
    if not prompt:
        return jsonify({'error': 'Missing prompt'}), 400
    # Create a new task for this message
    title = prompt[:60]  # Or use generate_title(prompt)
    workflow_id = 1 if model == 'toxindex-rap' else 2 if model == 'toxindex-vanilla' else 1
    task = Task.create_task(
        title=title,
        user_id=user_id,
        workflow_id=workflow_id,
        environment_id=environment_id,
        session_id=session_id,
    )
    # Add user message to session
    Message.create_message(task.task_id, user_id, 'user', prompt, session_id=session_id)
    # If this is the first user message, update the session title
    messages = Message.get_messages_by_session(session_id)
    user_messages = [m for m in messages if m.role == 'user']
    if len(user_messages) == 1:
        from webserver.model.chat_session import ChatSession as CSModel
        CSModel.update_title(session_id, prompt[:60])
        import logging
        logging.info(f"[DEBUG] Updated chat session {session_id} title to: {prompt[:60]}")
    # Trigger the appropriate celery task
    if model == 'toxindex-rap':
        celery_task = probra_task.delay({
            'payload': prompt,
            'task_id': task.task_id,
            'user_id': str(user_id),
        })
    elif model == 'toxindex-vanilla':
        celery_task = plain_openai_task.delay({
            'payload': prompt,
            'task_id': task.task_id,
            'user_id': str(user_id),
        })
    else:
        celery_task = probra_task.delay({
            'payload': prompt,
            'task_id': task.task_id,
            'user_id': str(user_id),
        })
    Task.update_celery_task_id(task.task_id, celery_task.id)
    return jsonify({'task_id': task.task_id, 'celery_id': celery_task.id})

@csrf.exempt
@app.route("/api/environment/<env_id>", methods=["GET"])
@flask_login.login_required
def api_environment_view(env_id):
    env = Environment.get_environment(env_id)
    if not env:
        return flask.jsonify({"error": "Environment not found"}), 404
    tasks = Task.get_tasks_by_environment(env_id, flask_login.current_user.user_id)
    files = File.get_files_by_environment(env_id)
    return flask.jsonify({
        "environment": env.to_dict(),
        "tasks": [t.to_dict() for t in tasks],
        "files": [f.to_dict() for f in files]
    })

@csrf.exempt
@app.route("/api/environment/<env_id>", methods=["DELETE"])
@flask_login.login_required
def api_environment_delete(env_id):
    if not is_valid_uuid(env_id):
        logging.error(f"Attempted to delete environment with invalid UUID: {env_id}")
        return flask.jsonify({"success": False, "error": "Invalid environment ID"}), 400
    try:
        # Disk cleanup: remove files from disk before deleting environment
        files = File.get_files_by_environment(env_id)
        for file in files:
            try:
                if file.filepath and os.path.exists(file.filepath):
                    os.remove(file.filepath)
                    logging.info(f"Deleted file from disk: {file.filepath}")
            except Exception as e:
                logging.warning(f"Failed to remove file from disk: {file.filepath}, error: {e}")
        # Now delete the environment (DB rows for files/tasks will be deleted by ON DELETE CASCADE)
        Environment.delete_environment(env_id, flask_login.current_user.user_id)
        return flask.jsonify({"success": True})
    except Exception as e:
        logging.error(f"Failed to delete environment {env_id}: {e}")
        return flask.jsonify({"success": False, "error": str(e)}), 500

@app.route("/api/user/<user_id>", methods=["GET"])
@flask_login.login_required
def api_get_user(user_id):
    user = User.get(user_id)
    if not user:
        return flask.jsonify({"error": "User not found"}), 404
    return flask.jsonify({"user_id": user.user_id, "email": user.email})

print("Registered routes:")
for rule in app.url_map.iter_rules():
    print(rule)

# --- ADD catch-all route for React SPA ---
@app.route('/', defaults={'path': ''})
@app.route('/<path:path>')
def serve_react_app(path):
    # Only serve index.html for non-API, non-static routes
    if path.startswith('api/') or path.startswith('static/') or path.startswith('uploads/'):
        abort(404)
    return send_from_directory(app.static_folder, 'index.html')

# LAUNCH =====================================================================
# In development, start the Redis listener only in the reloader child process
if __name__ == "__main__":
    # Only start the listener in the reloader child process (not parent)
    if os.environ.get("WERKZEUG_RUN_MAIN") == "true" or not app.debug:
        thread_name = f"RedisListenerThread-{uuid.uuid4().hex[:8]}"
        threading.Thread(
            target=redis_listener,
            args=(thread_name,),
            daemon=True,
            name=thread_name,
        ).start()
    socketio.run(app, host="0.0.0.0", port=6513, debug=True, use_reloader=True) # set to true if you want to reload while editing code.

# In production (Gunicorn), run the Redis listener as a separate process using redis_listener_service.py
# See redis_listener_service.py for details.

@csrf.exempt
@app.route('/api/chat_session/<session_id>', methods=['DELETE'])
@flask_login.login_required
def delete_chat_session(session_id):
    user_id = flask_login.current_user.user_id
    from webserver.model.chat_session import ChatSession
    ChatSession.delete_session(session_id, user_id)
    return flask.jsonify({'success': True})

@csrf.exempt
@app.route('/api/chat_session/<session_id>', methods=['PATCH'])
@flask_login.login_required
def rename_chat_session(session_id):
    user_id = flask_login.current_user.user_id
    data = flask.request.get_json(force=True)
    new_title = data.get('title')
    if not new_title:
        return flask.jsonify({'success': False, 'error': 'Missing title'}), 400
    from webserver.model.chat_session import ChatSession
    ChatSession.update_title(session_id, new_title)
    return flask.jsonify({'success': True})

from gevent import monkey
monkey.patch_all()

import dotenv
import uuid

dotenv.load_dotenv()

from webserver import login_manager as LM
from webserver.controller import login
from webserver.model import Task, Workflow, Message, File, Environment
from webserver.ai_service import generate_title
from workflows.chat_response_task import chat_response_task
from workflows.interactive_echo_task import interactive_echo_task
from flask import request, Response, send_from_directory

import flask, flask_login
import os, logging, requests, threading, redis, json
from datetime import datetime

# CELERY ======================================================================
from workflows.probra import probra_task
from workflows.chat_response_task import chat_response_task

from celery.result import AsyncResult
from workflows.celery_worker import celery

# SOCKETIO ====================================================================
from flask_socketio import SocketIO, emit
from flask_socketio import join_room

# FLASK APP ===================================================================
static_folder_path = os.path.join(os.path.dirname(__file__), "webserver", "static")
app = flask.Flask(__name__, template_folder="templates")
app.config["SERVER_NAME"] = os.environ.get("SERVER_NAME")
app.config["PREFERRED_URL_SCHEME"] = os.environ.get("PREFERRED_URL_SCHEME")
app.config["TEMPLATES_AUTO_RELOAD"] = True
app.secret_key = os.environ.get("FLASK_APP_SECRET_KEY")

socketio = SocketIO(
    app,
    cors_allowed_origins="*",
    message_queue="redis://localhost:6379/0",
    manage_session=False,
    async_mode='gevent',
    logger=True,
    engineio_logger=True
)

# Configure logging to output to file with detailed formatting
logging.basicConfig(
    level=logging.DEBUG,
    format='%(asctime)s - %(name)s - %(levelname)s - %(filename)s:%(lineno)d - %(message)s',
    handlers=[
        logging.FileHandler('app.log'),
        logging.StreamHandler()  # Also keep console output
    ]
)

# Set Flask app logger to use same configuration
app.logger.setLevel(logging.INFO)
for handler in logging.getLogger().handlers:
    app.logger.addHandler(handler)
LM.init(app)
Workflow.load_default_workflows()


# REDIS PUB/SUB BACKGROUND LISTENER ==========================================
def redis_listener(name):
    r = redis.Redis()
    pubsub = r.pubsub()
    pubsub.subscribe("celery_updates")

    for redis_message in pubsub.listen():
        if redis_message["type"] != "message":
            continue

        try:
            raw_data = redis_message["data"]
            if isinstance(raw_data, bytes):
                raw_data = raw_data.decode()

            event = json.loads(raw_data)
            logging.info(f"Redis event: {event}")
            
            event_type = event.get("type")
            event_task_id = event.get("task_id")
            event_data = event.get("data")

            if event_task_id is None or event_data is None:
                logging.warning(
                    f"Redis event missing required fields: event_type={event_type}, "
                    f"task_id={event_task_id}, data={event_data}"
                )
                continue

            if event_type not in ("task_message", "task_file"):
                continue

            task = Task.get_task(event_task_id)
            if not task:
                logging.warning(f"No task found for task_id={event_task_id}")
                continue

            if event_type == "task_message":
                Message.process_event(task, event_data)
            elif event_type == "task_file":
                logging.info(f"task_file {raw_data}")
                File.process_event(task, event_data)

            room = f"task_{task.task_id}"
            socketio.emit(event_type, event_data, to=room)
            logging.info(f"{event_type} sent to {room}")

        except json.JSONDecodeError:
            logging.error("Failed to decode Redis message as JSON.")
        except Exception as e:
            logging.error(f"Redis listener error: {e}", exc_info=True)


# INDEX / ROOT ===============================================================
@app.route("/", methods=["GET"])
def index():
    logging.info(f"current user: {flask_login.current_user}")
    if flask_login.current_user.is_authenticated:
        environments = Environment.get_environments_by_user(
            flask_login.current_user.user_id
        )
        return flask.render_template("index.html", environments=environments)
    else:
        print("user is not logged in")
        return flask.render_template("login_register.html")


# SOCKETIO HANDLERS ==========================================================
@socketio.on("connect")
def handle_connect(auth):
    emit("connected", {"sid": request.sid})
    logging.info(f"user connected {request.sid}")


# TODO right now I think any user can join any task room
@socketio.on("join_task_room")
def handle_join_task_room(data):
    task_id = data.get("task_id")
    room = f"task_{task_id}"
    join_room(room)
    logging.info(f"user joined task room: {room}")
    emit("joined_task_room", {"task_id": task_id})


# TASK MANAGEMENT ============================================================
@app.route("/task/new", methods=["POST"])
def task_create():
    logging.info(f"task_create called with {request.method}")
    task_data = request.get_json()
    message = task_data.get("message", "")

    title = generate_title(message)
    user_id = flask_login.current_user.user_id
    workflow_id = int(task_data.get("workflow", 1))
    environment_id = task_data.get("environment_id")
    sid = task_data.get("sid")
    task = Task.create_task(
        title=title,
        user_id=user_id,
        workflow_id=workflow_id,
        environment_id=environment_id,
    )

    Task.add_message(task.task_id, flask_login.current_user.user_id, "user", message)
    
    # Select the appropriate task based on workflow_id
    if workflow_id == 1:
        celery_task = probra_task.delay(
            {
                "payload": message,
                "sid": sid,
                "task_id": task.task_id,
                "user_id": str(user_id),
            }
        )
    elif workflow_id == 2:
        celery_task = chat_response_task.delay(
            {
                "payload": message,
                "sid": sid,
                "task_id": task.task_id,
                "user_id": str(user_id),
            }
        )
    elif workflow_id == 3:
        celery_task = interactive_echo_task.delay(
            {
                "payload": message,
                "sid": sid,
                "task_id": task.task_id,
                "user_id": str(user_id),
            }
        )
    else:
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
        tasks_data = [
            {
                "task_id": task.task_id,
                "title": task.title,
                "created_at": (
                    task.created_at.isoformat() if hasattr(task, "created_at") else None
                ),
                "workflow_id": task.workflow_id,
            }
            for task in tasks
        ]
        return flask.jsonify({"tasks": tasks_data})
    except Exception as e:
        logging.error(f"Error retrieving tasks: {str(e)}")
        return flask.jsonify({"error": "Failed to retrieve tasks"}), 500


@app.route("/task/<task_id>", methods=["GET"])
def task(task_id):
    messages = Message.get_messages(task_id)
    files = File.get_files(task_id)

    # Pre-render markdown content for convenience
    rendered_files = []
    for f in files:
        file_info = f.to_dict()
        if f.filename.lower().endswith((".md", ".markdown")):
            file_info["html"] = f.markdown_to_html()
        else:
            file_info["html"] = ""
        rendered_files.append(file_info)

    return flask.render_template(
        "task.html", task_id=task_id, messages=messages, files=rendered_files
    )


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
@app.route("/workflow/new", methods=["GET", "POST"])
def workflow_new():
    if request.method == "POST":
        workflow_data = request.get_json()
        title = workflow_data.get("title", "New workflow")
        description = workflow_data.get("description", "")
        workflow = Workflow.create_workflow(
            title=title,
            description=description,
            user_id=flask_login.current_user.user_id,
        )
        logging.info(f"created workflow {workflow.workflow_id}")
        return flask.jsonify({"workflow_id": workflow.workflow_id})
    return flask.render_template("create_workflow.html")


@app.route("/workflow/<workflow_id>", methods=["GET"])
def workflow(workflow_id):
    workflow = Workflow.get_workflow(workflow_id)
    return flask.render_template("workflow_run.html", workflow=workflow)


# ENVIRONMENT MANAGEMENT =====================================================
@app.route("/environment/new", methods=["POST"])
@flask_login.login_required
def environment_new():
    env_data = request.get_json()
    title = env_data.get("title", "New Environment")
    description = env_data.get("description")
    env = Environment.create_environment(
        title, flask_login.current_user.user_id, description
    )
    return flask.jsonify({"environment_id": env.environment_id})


@app.route("/environments", methods=["GET"])
@flask_login.login_required
def environment_list():
    envs = Environment.get_environments_by_user(flask_login.current_user.user_id)
    return flask.jsonify({"environments": [e.to_dict() for e in envs]})


@app.route("/environment/<env_id>", methods=["GET"])
@flask_login.login_required
def environment_view(env_id):
    workflow = Workflow.get_workflow(1)
    tasks = Task.get_tasks_by_environment(env_id, flask_login.current_user.user_id)
    env = Environment.get_environment(env_id)
    return flask.render_template(
        "environment.html", tasks=tasks, workflow=workflow, environment=env
    )


# LOGIN MANAGEMENT ===========================================================
app.route("/register", methods=["GET", "POST"])(login.register)
app.route("/verify", methods=["GET"])(login.verify_message)
app.route("/verification/<token>", methods=["GET", "POST"])(login.verification)
app.route("/login", methods=["GET", "POST"])(login.login)
app.route("/logout", methods=["GET"])(login.logout)
app.route("/forgot_password", methods=["GET", "POST"])(login.forgot_password)
app.route("/reset_password/<token>", methods=["GET", "POST"])(login.reset_password)


# ICONS ======================================================================
@app.route("/favicon.ico")
def favicon():
    return app.send_static_file("favicon.png")


@app.route("/icons/<path:filename>")
def serve_icon(filename):
    icons_dir = os.path.join(os.path.dirname(__file__), "icons")
    return send_from_directory(icons_dir, filename)


# LAUNCH =====================================================================
if __name__ == "__main__":
    if os.environ.get("WERKZEUG_RUN_MAIN") == "true":
        thread_name = f"RedisListenerThread-{uuid.uuid4().hex[:8]}"
        threading.Thread(
            target=redis_listener,
            args=(thread_name,),  # <- fix here
            daemon=True,
            name=thread_name,
        ).start()

    socketio.run(app, host="0.0.0.0", port=6513, debug=True, use_reloader=True)

from gevent import monkey
monkey.patch_all()

# Standard library imports
import os
import uuid
import threading
import logging
import json
from datetime import datetime
import redis
# Third-party imports
import dotenv
import flask
import flask_login
from flask_socketio import emit, join_room
from flask import request, send_from_directory, abort, make_response, jsonify

# Local imports
from webserver.login_manager import login_manager as LM
from webserver.model import Task, Workflow, Message, File, ChatSession
from webserver.ai_service import generate_title
from workflows.plain_openai_tasks import plain_openai_task
from workflows.probra import probra_task
from webserver.controller.environment import env_bp
from webserver.controller.task import task_bp
from webserver.controller.chat import chat_bp
from webserver.controller.auth import auth_bp
from webserver.controller.file import file_bp
from webserver.controller.workflow import workflow_bp
from webserver.controller.admin import admin_bp
from webserver.csrf import csrf
from webserver.socketio import socketio
from webserver.controller.user import user_bp
from webserver.controller.schema import schema_bp
from webserver.data_paths import LOGS_ROOT
from webserver.health_checker import health_checker

dotenv.load_dotenv()

# FLASK APP ===================================================================
static_folder_path = os.path.join(os.path.dirname(__file__), "webserver", "static")
app = flask.Flask(__name__, static_folder=static_folder_path)
# app.config["SERVER_NAME"] = os.environ.get("SERVER_NAME")
app.config["PREFERRED_URL_SCHEME"] = os.environ.get("PREFERRED_URL_SCHEME")
app.config["TEMPLATES_AUTO_RELOAD"] = True
app.secret_key = os.environ.get("FLASK_APP_SECRET_KEY")

logging.info(f"SECRET_KEY: {app.secret_key}")
logging.info(f"WTF_CSRF_ENABLED: {app.config.get('WTF_CSRF_ENABLED', 'not set')}")

socketio.init_app(app)

LOGS_ROOT().mkdir(exist_ok=True)
log_filename = LOGS_ROOT() / f"app_{datetime.now().strftime('%Y-%m-%d_%H')}.log"

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
LM.init_app(app)
Workflow.load_default_workflows()

csrf.init_app(app)

logging.info(f"DB ENV (Flask startup): PGHOST={os.getenv('PGHOST')}, PGPORT={os.getenv('PGPORT')}, PGDATABASE={os.getenv('PGDATABASE')}, PGUSER={os.getenv('PGUSER')}, PGPASSWORD={os.getenv('PGPASSWORD')}")

# REDIS PUB/SUB BACKGROUND LISTENER ==========================================
def redis_listener(name):
    
    listener_id = uuid.uuid4().hex[:8]
    r = redis.Redis(
        host=os.environ.get("REDIS_HOST", "localhost"),
        port=int(os.environ.get("REDIS_PORT", "6379"))
    )
    pubsub = r.pubsub()
    pubsub.subscribe("celery_updates")
    logging.info(f"[redis_listener] ({listener_id}) Started Redis listener thread: {name}")
    for redis_message in pubsub.listen():
        logging.info(f"[redis_listener] ({listener_id}) Received message: {redis_message}")
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
            if event_type not in ("task_message", "task_file", "task_status_update"):
                logging.debug(f"[redis_listener] ({listener_id}) Ignoring event_type={event_type}")
                continue
            task = Task.get_task(event_task_id)
            if not task:
                logging.warning(f"[redis_listener] ({listener_id}) No task found for task_id={event_task_id}")
                continue
            if event_type == "task_message":
                logging.info(f"[redis_listener] ({listener_id}) Processing task_message for task_id={event_task_id}")
                msg = Message.process_event(task, event_data)
                # Emit to chat session room if possible
                if msg and getattr(task, 'session_id', None):
                    chat_room = f"chat_session_{task.session_id}"
                    logging.info(f"[redis_listener] ({listener_id}) Emitting new_message to room {chat_room} with message: {msg.to_dict()}")
                    socketio.emit('new_message', msg.to_dict(), to=chat_room)
                    logging.info(f"[redis_listener] ({listener_id}) new_message sent to {chat_room}")
            elif event_type == "task_file":
                logging.info(f"[redis_listener] ({listener_id}) Processing task_file for task_id={event_task_id}")
                File.process_event(task, event_data)
            elif event_type == "task_status_update":
                logging.info(f"[redis_listener] ({listener_id}) Processing task_status_update for task_id={event_task_id}")
                # Emit a flat task object, not the event wrapper
                room = f"task_{event_task_id}"
                logging.info(f"[redis_listener] ({listener_id}) Emitting task_status_update to room {room} with data: {task.to_dict()}")
                socketio.emit("task_status_update", task.to_dict(), to=room)
                logging.info(f"[redis_listener] ({listener_id}) task_status_update sent to {room}")
            if event_type in ("task_message", "task_file"):
                room = f"task_{task.task_id}"
                logging.info(f"[redis_listener] ({listener_id}) Emitting {event_type} to room {room} with data: {event_data}")
                socketio.emit(event_type, event_data, to=room)
                logging.info(f"[redis_listener] ({listener_id}) {event_type} sent to {room}")
        except json.JSONDecodeError:
            logging.error(f"[redis_listener] ({listener_id}) Failed to decode Redis message as JSON.")
        except Exception as e:
            logging.error(f"[redis_listener] ({listener_id}) Redis listener error: {e}", exc_info=True)

# SOCKETIO HANDLERS ==========================================================
@socketio.on("connect")
def handle_connect(auth):
    logging.info(f"[socketio] Client connected: {request.sid}")
    emit("connected", {"sid": request.sid})


# TODO right now I think any user can join any task room
@socketio.on("join_task_room")
def handle_join_task_room(data):
    logging.info(f"JOIN ROOM {data}")
    task_id = data.get("task_id")
    user = flask_login.current_user

    # Check authentication
    if not user.is_authenticated:
        logging.warning(f"[socketio] join_task_room called without authentication by {request.sid}")
        emit("error", {"error": "Authentication required"})
        return

    # Check authorization: does this user own the task?
    task = Task.get_task(task_id)
    if not task or task.user_id != user.user_id:
        logging.warning(f"[socketio] join_task_room called without authorization by {request.sid}")
        emit("error", {"error": "Not authorized to join this task room"})
        return

    room = f"task_{task_id}"
    logging.info(f"[socketio] {request.sid} joining room: {room}")
    join_room(room)
    logging.info(f"[socketio] {request.sid} joined room: {room}")
    emit("joined_task_room", {"task_id": task_id})
    # Emit current status to the joining client
    logging.info(f"[socketio] Emitting current status for task {task_id} to {request.sid}: {task.to_dict()}")
    emit("task_status_update", task.to_dict(), room=request.sid)

@socketio.on("join_chat_session")
def handle_join_chat_session(data):
    session_id = data.get("session_id")
    if not session_id:
        logging.warning(f"[socketio] join_chat_session called without session_id by {request.sid}")
        return
    join_room(f"chat_session_{session_id}")
    logging.info(f"[socketio] {request.sid} joined chat_session_{session_id}")
    emit("joined_chat_session", {"session_id": session_id}, room=request.sid)

# Register Blueprints for modularized API endpoints
app.register_blueprint(env_bp)
app.register_blueprint(task_bp)
app.register_blueprint(chat_bp)
app.register_blueprint(file_bp)
app.register_blueprint(auth_bp)
app.register_blueprint(user_bp)
app.register_blueprint(schema_bp)
app.register_blueprint(workflow_bp)
app.register_blueprint(admin_bp)

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


@app.route("/api/healthz")
def healthz():
    """Basic health check endpoint."""
    app.logger.debug("Basic health check")
    return "ok", 200

@app.route("/api/health")
def comprehensive_health():
    """Comprehensive health check endpoint."""
    try:
        # Run all health checks
        results = health_checker.run_all_checks()
        
        # Determine overall status
        overall_status, http_status = health_checker.get_overall_status(results)
        
        # Prepare response
        response = {
            'status': overall_status,
            'http_status': http_status,
            'checks': results,
            'timestamp': datetime.now().isoformat()
        }
        
        # Add cache headers for health checks
        response_obj = jsonify(response)
        response_obj.headers['Cache-Control'] = 'no-cache, no-store, must-revalidate'
        response_obj.headers['Pragma'] = 'no-cache'
        response_obj.headers['Expires'] = '0'
        
        return response_obj, http_status
        
    except Exception as e:
        app.logger.error(f"Health check failed: {e}")
        return jsonify({
            'status': 'error',
            'message': f'Health check failed: {str(e)}',
            'timestamp': datetime.now().isoformat()
        }), 503

@app.route("/api/health/ready")
def readiness_probe():
    """Kubernetes readiness probe endpoint."""
    try:
        # Only check critical services for readiness
        critical_checks = {
            'database': health_checker._check_database,
            'redis': health_checker._check_redis,
        }
        
        results = {}
        for check_name, check_func in critical_checks.items():
            try:
                result = check_func()
                results[check_name] = result
            except Exception as e:
                results[check_name] = {
                    'status': 'error',
                    'message': str(e),
                    'timestamp': datetime.now().isoformat()
                }
        
        # Check if all critical services are healthy
        all_healthy = all(
            result.get('status') == 'healthy' 
            for result in results.values()
        )
        
        if all_healthy:
            return jsonify({
                'status': 'ready',
                'message': 'All critical services are ready',
                'timestamp': datetime.now().isoformat()
            }), 200
        else:
            return jsonify({
                'status': 'not_ready',
                'message': 'Critical services are not ready',
                'checks': results,
                'timestamp': datetime.now().isoformat()
            }), 503
            
    except Exception as e:
        app.logger.error(f"Readiness probe failed: {e}")
        return jsonify({
            'status': 'not_ready',
            'message': f'Readiness probe failed: {str(e)}',
            'timestamp': datetime.now().isoformat()
        }), 503

@app.route("/api/health/live")
def liveness_probe():
    """Kubernetes liveness probe endpoint."""
    try:
        # Simple check that the application is responding
        return jsonify({
            'status': 'alive',
            'message': 'Application is alive',
            'timestamp': datetime.now().isoformat()
        }), 200
    except Exception as e:
        app.logger.error(f"Liveness probe failed: {e}")
        return jsonify({
            'status': 'dead',
            'message': f'Application is not responding: {str(e)}',
            'timestamp': datetime.now().isoformat()
        }), 503

@app.route("/test-alive")
def test_alive():
    return "ALIVE"

# --- ADD catch-all route for React SPA ---
@app.route('/', defaults={'path': ''})
@app.route('/<path:path>')
def serve_react_app(path):
    # Only serve index.html for non-API, non-static routes
    # Don't intercept API routes that are already defined
    if path.startswith('static/') or path.startswith('uploads/'):
        abort(404)
    return send_from_directory(app.static_folder, 'index.html')

# --- Start Redis listener thread on every process startup (not just when run as __main__) ---
thread_name = f"RedisListenerThread-{uuid.uuid4().hex[:8]}"
threading.Thread(
    target=redis_listener,
    args=(thread_name,),
    daemon=True,
    name=thread_name,
).start()
# This ensures the listener runs in both dev (python app.py) and prod (gunicorn/uwsgi/etc)

# LAUNCH =====================================================================
if __name__ == "__main__":
    socketio.run(app, host="0.0.0.0", port=6513, debug=False, use_reloader=False) # Always debug=False in prod

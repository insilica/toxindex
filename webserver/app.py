from gevent import monkey
monkey.patch_all()

# Standard library imports
import os
import uuid
import threading
import logging
import json
import warnings
from datetime import datetime, timedelta
import redis

# Suppress specific warnings from third-party libraries
warnings.filterwarnings("ignore", message=".*is not.*int.*literal.*", category=SyntaxWarning)
warnings.filterwarnings("ignore", message=".*invalid escape sequence.*", category=SyntaxWarning)
# Third-party imports
import dotenv
import flask
import flask_login
from flask_socketio import emit, join_room
from flask import request, send_from_directory, abort, make_response, jsonify
from flask_cors import CORS

# Local imports
from webserver.login_manager import login_manager as LM
from webserver.model import Task, Workflow
from webserver.model.system_settings import SystemSettings
from webserver.controller.environment import env_bp
from webserver.controller.task import task_bp
from webserver.controller.chat import chat_bp
from webserver.controller.auth import auth_bp
from webserver.controller.file import file_bp
from webserver.controller.workflow import workflow_bp
from webserver.controller.admin import admin_bp
from webserver.controller.health import health_bp
from webserver.csrf import csrf
from webserver.socketio import socketio
from webserver.controller.user import user_bp
from webserver.controller.schema import schema_bp
from webserver.logging_utils import setup_logging, log_service_startup, get_logger

dotenv.load_dotenv()

# Setup logging with shared utility
setup_logging("webserver", log_level=logging.DEBUG)
logger = get_logger("webserver")

# FLASK APP ===================================================================
static_folder_path = os.path.join(os.path.dirname(__file__), "webserver", "static")
app = flask.Flask(__name__, static_folder=static_folder_path)
# app.config["SERVER_NAME"] = os.environ.get("SERVER_NAME")
app.config["PREFERRED_URL_SCHEME"] = os.environ.get("PREFERRED_URL_SCHEME")
app.config["TEMPLATES_AUTO_RELOAD"] = True
app.secret_key = os.environ.get("FLASK_APP_SECRET_KEY")

# Session timeout configuration (dynamic from database)
session_timeout_minutes = SystemSettings.get_setting_int('session_timeout_minutes', 60)
app.config["PERMANENT_SESSION_LIFETIME"] = timedelta(minutes=session_timeout_minutes)
app.config["SESSION_COOKIE_SECURE"] = True
app.config["SESSION_COOKIE_HTTPONLY"] = True
app.config["SESSION_COOKIE_SAMESITE"] = "Lax"

# CSRF Configuration
# Detect environment and configure CSRF appropriately
is_development = (
    os.environ.get('FLASK_ENV') == 'development' or 
    os.environ.get('FLASK_DEBUG') == '1'
)
is_localhost = (
    os.environ.get('PREFERRED_URL_SCHEME', 'https') == 'http'
)

# Initialize CORS only in development
if is_development or is_localhost:
    # Only enable CORS in development
    try:
        from flask_cors import CORS
        CORS(app, 
             origins=[
                 "http://localhost:5173",
                 "http://localhost:3000"
             ], 
             supports_credentials=True,
             allow_headers=["Content-Type", "Authorization", "X-Requested-With"],
             methods=["GET", "POST", "PUT", "DELETE", "OPTIONS"],
             expose_headers=["Content-Type", "Authorization"]
        )
        logger.info("CORS enabled for development environment")
    except ImportError:
        logger.warning("Flask-CORS not installed, CORS disabled")
else:
    logger.info("CORS disabled for production environment")

app.config["WTF_CSRF_TIME_LIMIT"] = 3600  # 1 hour

# Configure SSL settings based on environment
if is_development or is_localhost:
    # Local development: Allow HTTP, disable SSL strict
    app.config["WTF_CSRF_ENABLED"] = False
    app.config["WTF_CSRF_SSL_STRICT"] = False
    app.config["SESSION_COOKIE_SECURE"] = False  # Allow HTTP cookies
    logging.info(f"CSRF configured for development environment (is_development={is_development}, is_localhost={is_localhost})")
else:
    # Production: Require HTTPS, enable SSL strict
    app.config["WTF_CSRF_ENABLED"] = True
    app.config["WTF_CSRF_SSL_STRICT"] = True
    app.config["SESSION_COOKIE_SECURE"] = True  # Require HTTPS cookies
    logging.info(f"CSRF configured for production environment (is_development={is_development}, is_localhost={is_localhost})")

logging.info(f"SECRET_KEY: {app.secret_key}")
logging.info(f"WTF_CSRF_ENABLED: {app.config.get('WTF_CSRF_ENABLED', 'not set')}")
logging.info(f"Session timeout: {app.config['PERMANENT_SESSION_LIFETIME']} (from database: {session_timeout_minutes} minutes)")

socketio.init_app(app)

# Log startup information
log_service_startup("webserver")

# Set Flask app logger to use same configuration
app.logger.setLevel(logging.INFO)
for handler in logging.getLogger().handlers:
    app.logger.addHandler(handler)
LM.init_app(app)
Workflow.load_default_workflows()

csrf.init_app(app)

logger.info(f"DB ENV (Flask startup): PGHOST={os.getenv('PGHOST')}, PGPORT={os.getenv('PGPORT')}, PGDATABASE={os.getenv('PGDATABASE')}, PGUSER={os.getenv('PGUSER')}, PGPASSWORD={os.getenv('PGPASSWORD')}")

# SOCKETIO HANDLERS ==========================================================
@socketio.on("connect")
def handle_connect(auth):
    try:
        logger.info(f"[socketio] Client connected: {request.sid}")
        # Don't require authentication for initial connection
        # Authentication will be checked when joining specific rooms
        emit("connected", {"sid": request.sid})
    except Exception as e:
        logger.error(f"[socketio] Error in connect handler: {e}")
        emit("error", {"error": "Connection failed"})

@socketio.on_error()
def error_handler(e):
    """Global error handler for Socket.IO events."""
    error_str = str(e)
    logger.error(f"[socketio] Socket.IO error: {error_str}")
    
    # Handle specific payload errors
    if "Too many packets in payload" in error_str:
        logger.warning(f"[socketio] Payload size exceeded for {request.sid}, reducing data size")
        try:
            emit("error", {"error": "Payload too large, reducing data size"})
        except:
            pass  # Don't emit if we can't
    elif "Payload decode error" in error_str:
        logger.warning(f"[socketio] Payload decode error for {request.sid}")
        try:
            emit("error", {"error": "Invalid payload format"})
        except:
            pass
    else:
        try:
            emit("error", {"error": "An error occurred"})
        except:
            pass  # Don't emit if we can't

@socketio.on("disconnect")
def handle_disconnect(data=None):
    try:
        logger.info(f"[socketio] Client disconnected: {request.sid}")
    except Exception as e:
        logger.error(f"[socketio] Error in disconnect handler: {e}")

# TODO right now I think any user can join any task room
@socketio.on("join_task_room")
def handle_join_task_room(data):
    try:
        logger.info(f"[socketio] join_task_room called with data: {data}")
        task_id = data.get("task_id")
        user = flask_login.current_user

        # Check authentication
        if not user.is_authenticated:
            logger.warning(f"[socketio] join_task_room called without authentication by {request.sid}")
            emit("error", {"error": "Authentication required"})
            return

        # Check authorization: does this user own the task?
        task = Task.get_task(task_id)
        if not task or task.user_id != user.user_id:
            logger.warning(f"[socketio] join_task_room called without authorization by {request.sid} for task {task_id}")
            emit("error", {"error": "Not authorized to join this task room"})
            return

        room = f"task_{task_id}"
        logger.info(f"[socketio] {request.sid} joining room: {room}")
        join_room(room)
        logger.info(f"[socketio] {request.sid} joined room: {room}")
        emit("joined_task_room", {"task_id": task_id})
        
        # Send minimal task status to avoid payload issues
        try:
            task_dict = task.to_dict()
            # Remove potentially large fields to reduce payload size
            minimal_task = {
                "task_id": task_dict.get("task_id"),
                "status": task_dict.get("status"),
                "title": task_dict.get("title"),
                "created_at": task_dict.get("created_at")
            }
            logging.info(f"[socketio] Emitting minimal task status for task {task_id} to {request.sid}")
            emit("task_status_update", minimal_task, room=request.sid)
        except Exception as e:
            logging.error(f"[socketio] Error sending task status: {e}")
            emit("error", {"error": "Failed to send task status"})
            
    except Exception as e:
        logging.error(f"[socketio] Error in join_task_room handler: {e}")
        emit("error", {"error": "Failed to join task room"})

@socketio.on("join_chat_session")
def handle_join_chat_session(data):
    try:
        logging.info(f"[socketio] join_chat_session called with data: {data}")
        session_id = data.get("session_id")
        if not session_id:
            logging.warning(f"[socketio] join_chat_session called without session_id by {request.sid}")
            emit("error", {"error": "session_id is required"})
            return
        
        # Check if user is authenticated
        user = flask_login.current_user
        if not user.is_authenticated:
            logging.warning(f"[socketio] join_chat_session called without authentication by {request.sid}")
            emit("error", {"error": "Authentication required"})
            return
        
        room = f"chat_session_{session_id}"
        logging.info(f"[socketio] {request.sid} joining room: {room}")
        join_room(room)
        logging.info(f"[socketio] {request.sid} joined chat_session_{session_id}")
        emit("joined_chat_session", {"session_id": session_id}, room=request.sid)
    except Exception as e:
        logging.error(f"[socketio] Error in join_chat_session handler: {e}")
        emit("error", {"error": "Failed to join chat session"})

@socketio.on("leave_task_room")
def handle_leave_task_room(data):
    try:
        logging.info(f"[socketio] leave_task_room called with data: {data}")
        task_id = data.get("task_id")
        if task_id:
            room = f"task_{task_id}"
            logging.info(f"[socketio] {request.sid} leaving room: {room}")
            # Note: Flask-SocketIO doesn't have a direct leave_room method for clients
            # The room will be cleaned up automatically when the client disconnects
    except Exception as e:
        logging.error(f"[socketio] Error in leave_task_room handler: {e}")

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
app.register_blueprint(health_bp)

# SESSION TIMEOUT MIDDLEWARE =================================================
@app.before_request
def check_session_timeout():
    """Check if user session has expired and log them out if necessary."""
    
    # Skip session checks for auth endpoints (login, register, etc.)
    if request.endpoint and request.endpoint.startswith('auth.'):
        return
    
    # Skip session checks for static files and other non-API routes
    if request.path.startswith('/static/') or request.path.startswith('/uploads/'):
        return
    
    if flask_login.current_user.is_authenticated:
        try:
            # Import the session management functions
            from webserver.controller.auth import is_session_expired, update_session_activity
            
            user_id = flask_login.current_user.user_id
            
            # Check if session has expired in Redis
            if is_session_expired(user_id):
                logging.info(f"[session] Session expired for user {flask_login.current_user.email}")
                flask_login.logout_user()
                flask.session.clear()
                return jsonify({'error': 'Session expired. Please log in again.'}), 401
            
            # Update activity for this request
            update_session_activity(user_id)
            
        except Exception as e:
            logging.error(f"[session] Error checking session timeout: {e}")
            # Fall back to old Flask session check for compatibility
            session_timeout_minutes = SystemSettings.get_setting_int('session_timeout_minutes', 60)
            session_lifetime = timedelta(minutes=session_timeout_minutes)
            
            if flask.session.get('_permanent', False):
                now = datetime.now()
                created = flask.session.get('_created', now)
                
                if hasattr(created, 'tzinfo') and created.tzinfo is not None:
                    created = created.replace(tzinfo=None)
                
                session_age = now - created
                if session_age > session_lifetime:
                    logging.info(f"[session] Session expired for user {flask_login.current_user.email}")
                    flask_login.logout_user()
                    flask.session.clear()
                    return jsonify({'error': 'Session expired. Please log in again.'}), 401

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
    user_id = flask_login.current_user.user_id if flask_login.current_user.is_authenticated else None
    logging.info(f"Tab switch: user={user_id}, task={task_id}, tab={tab}")
    return jsonify({"success": True})

# HEALTH ENDPOINTS ===========================================================
@app.route("/api/healthz")
def healthz():
    """Simple health check endpoint."""
    return jsonify({"status": "healthy"}), 200

# --- ADD catch-all route for React SPA ---
@app.route('/', defaults={'path': ''})
@app.route('/<path:path>')
def serve_react_app(path):
    # Only serve index.html for non-API, non-static routes
    # Don't intercept API routes that are already defined
    if path.startswith('static/') or path.startswith('uploads/'):
        abort(404)
    return send_from_directory(app.static_folder, 'index.html')

# --- Redis listener is now in a separate deployment ---
# The main app no longer starts Redis listeners
# This prevents duplicate listeners across multiple web pods

logging.info("Redis listener disabled - using separate deployment for background processing")

# LAUNCH =====================================================================
if __name__ == "__main__":
    socketio.run(app, host="0.0.0.0", port=6513, debug=False, use_reloader=False) # Always debug=False in prod

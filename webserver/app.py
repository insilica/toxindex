import eventlet
eventlet.monkey_patch()

import dotenv
dotenv.load_dotenv()

from webserver import login_manager as LM
from webserver.controller import login, stripe
from webserver.model import Task, Workflow
from webserver.ai_service import generate_title
from flask import request, Response, send_from_directory

import flask, flask_login
import os, logging, requests, threading, redis, json
from datetime import datetime

# CELERY ======================================================================
from workflows.probra import probra_task
from celery.result import AsyncResult
from workflows.celery_worker import celery

# SOCKETIO ====================================================================
from flask_socketio import SocketIO, emit

# FLASK APP ===================================================================
static_folder_path = os.path.join(os.path.dirname(__file__), 'webserver', 'static')
app = flask.Flask(__name__, template_folder="templates")
app.config['SERVER_NAME'] = os.environ.get('SERVER_NAME')
app.config['PREFERRED_URL_SCHEME'] = os.environ.get('PREFERRED_URL_SCHEME')
app.secret_key = os.environ.get('FLASK_APP_SECRET_KEY')

socketio = SocketIO(app, cors_allowed_origins="*", message_queue='redis://localhost:6379/0')

app.logger.setLevel(logging.INFO)
logging.basicConfig(level=logging.DEBUG)
LM.init(app)
Workflow.load_default_workflows()

# REDIS PUB/SUB BACKGROUND LISTENER ==========================================
def redis_listener():
    r = redis.Redis()
    pubsub = r.pubsub()
    pubsub.subscribe("celery_updates")

    for message in pubsub.listen():
        if message["type"] != "message":
            continue
        try:
            data = json.loads(message["data"])
            sid = data["sid"]
            event = data["event"]
            payload = data["data"]
            socketio.emit(event, payload, to=sid)
        except Exception as e:
            logging.error(f"Redis listener error: {e}")

threading.Thread(target=redis_listener, daemon=True).start()

# INDEX / ROOT ===============================================================
@app.route('/', methods=['GET'])
def index():
    logging.info(f"current user: {flask_login.current_user}")
    if flask_login.current_user.is_authenticated:
        tasks = Task.get_tasks_by_user(flask_login.current_user.user_id)
        workflow = Workflow.get_workflow(1)
        return flask.render_template('index.html', tasks=tasks, workflow=workflow)
    else:
        print('user is not logged in')
        return flask.render_template('login_register.html')

# SOCKETIO HANDLERS ==========================================================
@socketio.on('connect')
def handle_connect():
    logging.info(f"Client connected: {request.sid}")
    emit('connected', {'sid': request.sid})

# TASK MANAGEMENT ============================================================
@app.route('/task/new', methods=['POST'])
def task_create():
    logging.info(f"task_create called with {request.method}")

    if request.method == 'POST':
        task_data = request.get_json()
        message = task_data.get('message', '')
        workflow_id = int(task_data.get('workflow', 1))
        sid = task_data.get('sid')

        title = generate_title(message)

        celery_task = probra_task.delay({
            "payload": message,
            "sid": sid
        })

        task = Task.create_task(
            title=title,
            user_id=flask_login.current_user.user_id,
            workflow_id=workflow_id,
            celery_task_id=celery_task.id
        )

        Task.add_message(task.task_id, flask_login.current_user.user_id, 'user', message)

        return flask.jsonify({'task_id': task.task_id, 'celery_id': celery_task.id})
    
    return flask.render_template('task.html')

@app.route('/task/status/<task_id>', methods=['GET'])
def task_status(task_id):
    # TODO You might want to load your Task object here and link it to a celery_id
    celery_id = request.args.get('celery_id')  # Pass celery_id via query param or load from DB
    if not celery_id:
        return flask.jsonify({'error': 'Missing celery_id'}), 400

    result = AsyncResult(celery_id, app=celery)
    return flask.jsonify({
        "state": result.state,
        "info": result.info if result.info else {},
        "ready": result.ready(),
        "result": result.result if result.ready() else None,
    })

@app.route('/tasks', methods=['GET'])
@flask_login.login_required
def get_user_tasks():
    try:
        user_id = flask_login.current_user.user_id
        tasks = Task.get_tasks_by_user(user_id)
        tasks_data = [{
            'task_id': task.task_id,
            'title': task.title,
            'created_at': task.created_at.isoformat() if hasattr(task, 'created_at') else None,
            'workflow_id': task.workflow_id
        } for task in tasks]
        return flask.jsonify({'tasks': tasks_data})
    except Exception as e:
        logging.error(f"Error retrieving tasks: {str(e)}")
        return flask.jsonify({'error': 'Failed to retrieve tasks'}), 500

@app.route('/task/<task_id>', methods=['GET'])
def task(task_id):
    return flask.render_template('task.html', task_id=task_id)

# WORKFLOW MANAGEMENT ========================================================
@app.route('/workflow/new', methods=['GET', 'POST'])
def workflow_new():
    if request.method == 'POST':
        workflow_data = request.get_json()
        title = workflow_data.get('title', 'New workflow')
        description = workflow_data.get('description', '')
        workflow = Workflow.create_workflow(
            title=title,
            description=description,
            user_id=flask_login.current_user.user_id
        )
        logging.info(f"created workflow {workflow.workflow_id}")
        return flask.jsonify({'workflow_id': workflow.workflow_id})
    return flask.render_template('create_workflow.html')

@app.route('/workflow/<workflow_id>', methods=['GET'])
def workflow(workflow_id):
    workflow = Workflow.get_workflow(workflow_id)
    return flask.render_template('workflow_run.html', workflow=workflow)

# LOGIN MANAGEMENT ===========================================================
app.route('/register', methods=['GET','POST'])(login.register)
app.route('/verify', methods=['GET'])(login.verify_message)
app.route('/verification/<token>', methods=['GET','POST'])(login.verification)
app.route('/login', methods=['GET','POST'])(login.login)
app.route('/logout', methods=['GET'])(login.logout)
app.route('/forgot_password', methods=['GET','POST'])(login.forgot_password)
app.route('/reset_password/<token>', methods=['GET','POST'])(login.reset_password)

# ICONS ======================================================================
@app.route('/favicon.ico')
def favicon():
    return app.send_static_file('favicon.png')

@app.route('/icons/<path:filename>')
def serve_icon(filename):
    icons_dir = os.path.join(os.path.dirname(__file__), 'icons')
    return send_from_directory(icons_dir, filename)

# LAUNCH =====================================================================
if __name__ == '__main__':
    app.secret_key
    socketio.run(app, host="0.0.0.0", port=6513)

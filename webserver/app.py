from webserver import login_manager as LM
from webserver.controller import login, stripe
from flask import request, Response, send_from_directory

import flask, flask_login
import os, logging, requests
from datetime import datetime

# Defining the static folder path
static_folder_path = os.path.join(os.path.dirname(__file__), 'webserver', 'static')
app = flask.Flask(__name__, template_folder="templates")

app.config['SERVER_NAME'] = os.environ.get('SERVER_NAME')
app.config['PREFERRED_URL_SCHEME'] = os.environ.get('PREFERRED_URL_SCHEME')
app.secret_key =  os.environ.get('FLASK_APP_SECRET_KEY')
app.logger.setLevel(logging.INFO)
logging.basicConfig(level=logging.DEBUG) # Use INFO for less verbose output

LM.init(app)

tasks = [
            {"id": "1", "title": "Pesticide exposure assessment"},
            {"id": "2", "title": "Heavy metal contamination analysis"},
            {"id": "3", "title": "Drug interaction toxicity report"}
        ]
        
@app.route('/', methods=['GET'])
def index():
  logging.info(f"current user: {flask_login.current_user}")
  if flask_login.current_user.is_authenticated:
    return flask.render_template('chat.html', tasks=tasks)
  else:
    print('user is not logged in')
    return flask.render_template('landing.html')

# Tasks
@app.route('/task/<task_id>', methods=['GET', 'POST'])
@flask_login.login_required
def task(task_id):
    task_data = {
        "id": task_id,
        "title": f"Task {task_id}",
        "created_at": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        "messages": []
    }
    
    return flask.render_template('task.html', task=task_data, tasks=tasks)

@app.route('/task/new', methods=['POST'])
@flask_login.login_required
def create_task():
    data = request.get_json()
    title = data.get('title')
    agent = data.get('agent')
    
    # In a real app, create task in database
    # For now, mock creating a new task
    new_task_id = str(len(tasks) + 1)
    new_task = {
        "id": new_task_id,
        "title": title
    }
    tasks.append(new_task)
    
    return flask.jsonify({
        "task_id": new_task_id
    })

# Login and Registration
app.route('/register', methods=['GET','POST'])(login.register)
app.route('/verify', methods=['GET'])(login.verify_message)
app.route('/verification/<token>', methods=['GET','POST'])(login.verification)
app.route('/login', methods=['GET','POST'])(login.login)
app.route('/logout', methods=['GET'])(login.logout)

app.route('/forgot_password', methods=['GET','POST'])(login.forgot_password)
app.route('/reset_password/<token>', methods=['GET','POST'])(login.reset_password)
# Serving the favicon
@app.route('/favicon.ico')
def favicon():
    return app.send_static_file('favicon.png')

# Serving SVG icons from webserver/icons directory
@app.route('/icons/<path:filename>')
def serve_icon(filename):
    icons_dir = os.path.join(os.path.dirname(__file__), 'icons')
    return send_from_directory(icons_dir, filename)

from webserver import login_manager as LM
from webserver.controller import login, stripe
from webserver.model import Run, Workflow
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
Workflow.load_default_workflows()
        
@app.route('/', methods=['GET'])
def index():
  logging.info(f"current user: {flask_login.current_user}")
  runs = Run.get_runs_by_user(flask_login.current_user.user_id)
  if flask_login.current_user.is_authenticated:
    return flask.render_template('landing.html', runs=runs)
  else:
    print('user is not logged in')
    return flask.render_template('login_register.html')

@app.route('/run/new', methods=['POST'])
def run_create():
    logging.info(f"run_create called with {request.method}")
    
    if request.method == 'POST':
        run_data = request.get_json()
        title = run_data.get('title', 'New run')
        workflow_id = int(run_data.get('workflow', 1))
        
        # Create new run
        run = Run.create_run(
            title=title,
            user_id=flask_login.current_user.user_id,
            workflow_id=workflow_id
        )
        
        return flask.jsonify({'run_id': run.run_id})
    return flask.render_template('run.html')

@app.route('/run/<run_id>', methods=['GET'])
def run(run_id):
    return flask.render_template('run.html', run_id=run_id)
  
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
	return flask.render_template('workflow.html', workflow=workflow)

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

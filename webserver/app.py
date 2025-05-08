from webserver import login_manager as LM
from webserver.controller import login, stripe
from webserver.model.project import Project
from flask import request, Response

import flask, flask_login
import os, logging, requests

# Defining the static folder path
static_folder_path = os.path.join(os.path.dirname(__file__), 'webserver', 'static')
app = flask.Flask(__name__, template_folder="templates")

app.config['SERVER_NAME'] = os.environ.get('SERVER_NAME')
app.config['PREFERRED_URL_SCHEME'] = os.environ.get('PREFERRED_URL_SCHEME')
app.secret_key =  os.environ.get('FLASK_APP_SECRET_KEY')
app.logger.setLevel(logging.INFO)
logging.basicConfig(level=logging.DEBUG) # Use INFO for less verbose output

LM.init(app)

@app.route('/', methods=['GET'])
def index():
  logging.info(f"current user: {flask_login.current_user}")
  if flask_login.current_user.is_authenticated:
    active_projects = [p.to_dict() for p in Project.get_projects_by_creator(flask_login.current_user.user_id)]
    active_project = active_projects[0]["project_id"] if len(active_projects) > 0 else None
    if not active_project:
      res = Project.create_project("default", "user default project", flask_login.current_user.user_id)
      active_projects = active_projects + [res]
      active_project = res.project_id
      
    logging.info(f"active projects: {active_projects}")
    logging.info(f"active project: {active_project}")
    return flask.redirect(f'p/{active_project}/report')
    # return flask.render_template('layout.html', projects=active_projects, active_project=active_project)
  else:
    print('user is not logged in')
    return flask.render_template('landing.html')

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

# PROJECTS ====================================================================================
@app.route('/create_new_project', methods=['POST'])
def create_new_project():
    data = request.get_json()
    logging.info(f"Creating new project with name: {data['name']} for user {flask_login.current_user}")
    Project.create_project(data['name'], data.get("description"), flask_login.current_user.user_id)
    return flask.jsonify(success=True, message="Project created successfully.")

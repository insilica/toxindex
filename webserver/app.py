from webserver import login_manager as LM
from webserver.controller import login, stripe
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
    return flask.render_template('layout.html')
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

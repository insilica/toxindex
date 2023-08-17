from webserver import login_manager as LM
from webserver.controller import login

import requests
from flask import request, Response
import logging
import flask
import flask_login
import os


# Defining the static folder path
static_folder_path = os.path.join(os.path.dirname(__file__), 'webserver', 'static')
app = flask.Flask(__name__, template_folder="templates")
app.logger.setLevel(logging.INFO)
app.secret_key = 'super secret string'  # TODO Change this!
LM.init(app)

@app.route('/', methods=['GET'])
def index():
  if flask_login.current_user.is_authenticated:
    print(f'user {flask_login.current_user} is logged in')
    return flask.render_template('index.html')
  else:
    print('user is not logged in')
    return flask.redirect('login')
  

# Login and Registration
app.route('/register', methods=['GET','POST'])(login.register)
app.route('/verification/<token>', methods=['GET','POST'])(login.verification)
app.route('/login', methods=['GET','POST'])(login.login)
app.route('/logout', methods=['GET'])(login.logout)

app.route('/forgot_password', methods=['GET','POST'])(login.forgot_password)
app.route('/reset_password/<token>', methods=['GET','POST'])(login.reset_password)

# @app.route('/customer_portal', methods=['GET'])
# @flask_login.login_required
# def customer_portal():
#   customer_id = flask_login.current_user.stripe_customer_id
#   session = ToxindexStripe.create_customer_portal_session(customer_id)
#   return flask.redirect(session.url)

@app.route('/token', methods=['GET'])
@flask_login.login_required
def token():
  token = flask_login.current_user.token
  return f"your token is below. please do not share it with anyone <br><b>{token}</b>"

# Serving the favicon
@app.route('/favicon.ico')
def favicon():
    return app.send_static_file('favicon.png')

def active_projects():
  return [{"name":'test'},{"name":'test2'},{"name":"test3"}]

# OVERVIEW =====================================================================================
@app.route('/p/<project>/overview', methods=['GET'])
def overview(project):
    return flask.render_template('services/overview/overview.html', 
                                 active_project = project, 
                                 active_service = 'overview',
                                 projects = active_projects())
  
# UPLOAD =======================================================================================
# Directory to save uploaded files
@app.route('/p/<project>/upload', methods=['GET'])
def upload_get(project):
    return flask.render_template('services/upload/upload.html', 
                              active_project = project, 
                              active_service = 'upload',
                              projects = active_projects())


@app.route('/p/<project>/upload', methods=['POST'])
def upload_post(project):
    # Get the file from the request
    file = request.files['file']

    # You can now handle the file as needed (e.g., store in database, analyze, etc.)
    return flask.jsonify({"message": "File uploaded successfully", "file_path": file.filename})
  
# CHEMICALS ====================================================================================
@app.route('/p/<project>/chemicals', methods=['GET'])
def chemicals():
    template_name = 'services/upload/upload.html'
    return flask.render_template(template_name)

# REPORTS ====================================================================================
@app.route('/p/<project>/<service>', methods=['GET'])
@app.route('/p/<project>/<service>/<path:path>', methods=['GET'])
def reports(project,service,path=""):
    
    response = requests.get(f'http://{service}:6515/{path}')

    if response.headers['Content-Type'].startswith('text/html'):
        return flask.render_template('services/layout.html', 
                                     content=response.content.decode('utf-8'), 
                                     projects=active_projects(), 
                                     active_project=project,
                                     active_service='report')
    
    return Response(response.content, 
                        response.status_code,
                        dict(response.headers))

@app.route('/p/<project>/<service>', methods=['POST'])
@app.route('/p/<project>/<service>/<path:path>', methods=['POST'])
def service_post(project,service,path=""):
    
    
    service_url = f'http://{service}:6515/{path}'
    
    response = requests.request(
        method=request.method,
        url=service_url,
        headers={key: value for (key, value) in request.headers if key != 'Host'},
        data=request.get_data(),
        cookies=request.cookies,
        params=request.args
    )
    
    proxied_response = app.response_class(
        response=response.content,
        status=response.status_code,
        headers=dict(response.headers)
    )
    
    return proxied_response
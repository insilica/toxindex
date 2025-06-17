from webserver.model.user import User
from webserver.forms.registration_form import RegistrationForm, ForgotPasswordForm, UpdatePasswordForm
from webserver.controller import sendgrid
from webserver import datastore as ds

from flask import url_for

import flask, flask_login, secrets, datetime, logging
import pprint


login_page = flask.Blueprint('login', __name__, template_folder='templates')

# EMAIL VERIFICATION ====================================================================================
def generate_validation_link(user):
  link_token = secrets.token_urlsafe(16)
  expiration = datetime.datetime.now() + datetime.timedelta(days=1)
  params = (user.user_id, link_token, expiration)
  ds.execute('INSERT INTO user_link (user_id,link_token,expiration) values (%s,%s,%s)',params)
  return link_token

def validate(user: User):
  link = generate_validation_link(user)
  route = url_for('verification', token=link, _external=True)
  msg = (
    f"<p>Dear {user.name if user.name else 'User'},</p>"
    "<p>Thank you for registering with Toxindex.com! To complete your registration and activate your account, please click on the link below:</p>"
    f"<p><a href='{route}'><b>Activate My Account</b></a></p>"
    "<p>If you did not request this verification, please ignore this email.</p>"
    "<p>Best regards,<br>The Toxindex Team</p>"
  )
  sendgrid.send(user.email, "Activate Your Toxindex.com Account", msg)



@login_page.route('/verify_message', methods=['GET'])
def verify_message():
    return flask.render_template('verify_message.html', email=flask.request.args.get('email'))

def verification(token):
  user_id = ds.find('SELECT user_id from user_link WHERE link_token=(%s)', (token,))['user_id']
  ds.execute('UPDATE users SET email_verified = TRUE WHERE user_id = (%s)', (user_id,))
  if User.get(user_id) is not None:
    flask_login.login_user(User.get(user_id))
  return flask.redirect('/')

# USER PASSWORD RESET ===================================================================================
def send_password_reset(email):
  link = generate_validation_link(User.get_user(email))
  route = url_for("reset_password", token=link, _external=True) # Added token parameter and _external=True
  msg = f"Reset your password at this <a href='{route}'>link</a>"
  sendgrid.send(email,"Toxindex.com Password Reset",msg)

def forgot_password():
  form = ForgotPasswordForm()
  if flask.request.method == "GET": 
    return flask.render_template('forgot_password.html', form=form, sent_email=False)
  
  email = flask.request.form['email']
  send_password_reset(email)
  
  return flask.render_template('forgot_password.html', form=form, sent_email=True) 

def reset_password(token):
  form = UpdatePasswordForm()

  if flask.request.method == "GET": 
    return flask.render_template('reset_password.html', form=form, token=token)
  
  passw = flask.request.form['password']
  user_id = ds.find('SELECT user_id from user_link WHERE link_token=(%s)',(token,))['user_id']
  ds.execute('UPDATE users set password = (%s) WHERE user_id = (%s)',(passw,user_id,))
  
  if User.get(user_id) is not None:
    flask_login.login_user(User.get(user_id))

  return flask.redirect('/')

# USER Registration =====================================================================================
def register():
    is_api = flask.request.path.startswith('/api/')
    if is_api:
        form = RegistrationForm(meta={'csrf': False})
    else:
        form = RegistrationForm()

    logging.debug(f"[register] Request method: {flask.request.method}")
    logging.debug(f"[register] Raw form data: {flask.request.form}")
    logging.debug(f"[register] Raw data: {flask.request.data}")
    logging.debug(f"[register] Form populated data: email={form.email.data}, password={form.password.data}, password_confirmation={form.password_confirmation.data}")

    if form.validate_on_submit(): # Validates form data
        logging.debug(f"[register] Form validated successfully for email: {form.email.data}")
        if User.user_exists(form.email.data):
            flask.flash(f"{form.email.data} already exists, please log in")
            logging.debug(f"[register] User already exists: {form.email.data}")
            return flask.redirect('register')

        user = User.create_user(form.email.data, form.password.data)

        if user is not None: # Checks if the user object is valid
            logging.debug(f"[register] User created: {user.email}")
            validate(user)
            return flask.redirect(url_for('verify_message', email=form.email.data))
        else:
          logging.debug('[register] user is none!')
          flask.flash("An error occurred while creating the user. Please try again.")
    else:
        logging.debug('[register] form did not validate')
        logging.debug(f"[register] form.errors: {pprint.pformat(form.errors)}")
        for field, errors in form.errors.items():
            for error in errors:
                logging.debug(f"[register] Field '{field}' error: {error}")
        flask.flash("An error occurred while creating the user. Please try again.")
        
    return flask.redirect("/login")

# USER LOGIN/OUT ==================================================================================
def login():

  if flask.request.method == "GET": return flask.render_template('login.html')

  email = flask.request.form['email']
  passw = flask.request.form['password']
  user  = User.get_user(email)

  if not User.user_exists(email) or not user.validate_password(passw):
    return "bad login"

  flask_login.login_user(user)
  return flask.redirect("/")

def logout():
  flask_login.logout_user()
  flask.session.clear()  # Explicitly clear the session
  return flask.jsonify({"success": True})
from webserver.model.user import User
from webserver.controller import sendgrid
from webserver import datastore as ds
from webserver.csrf import csrf
from flask import url_for, request, jsonify, Blueprint
import flask, flask_login, secrets, datetime, os
import random

auth_bp = Blueprint('auth', __name__, url_prefix='/api/auth')
FRONTEND_URL = os.environ.get('FRONTEND_URL', 'http://localhost:5173')

ERROR_MESSAGES = [
    "Oops! That combo didn't work. Try again?",
    "No luck! Double-check your email and password.",
    "That's not quite right. Give it another shot!",
    "Hmm, those credentials didn't match. Want to try again?",
    "Access denied! But don't give up!",
    "Almost! Check your email and password and try again.",
    "Whoops! That's not the right combo. Try once more?",
    "The login gods say no. Try again?",
    "Not quite! Maybe a typo snuck in?",
    "That didn't work, but you've got this!",
    "Passwords are tricky! Want to try again?",
    "Nope, not this time. But you're getting closer!",
    "The login fairies are not convinced. Try again!",
    "That wasn't it, but persistence pays off!",
    "Keep going! Even the best forget sometimes.",
    "Still locked out! Maybe a coffee break?",
    "Almost cracked it! One more try?",
    "The universe says: not quite. Try again!",
    "You shall not pass... yet! Try again.",
    "Mistakes happen! Give it another go."
]

def generate_validation_link(user):
  link_token = secrets.token_urlsafe(16)
  expiration = datetime.datetime.now() + datetime.timedelta(days=1)
  params = (user.user_id, link_token, expiration)
  ds.execute('INSERT INTO user_link (user_id,link_token,expiration) values (%s,%s,%s)',params)
  return link_token

def validate(user: User):
  link = generate_validation_link(user)
  # Use the React SPA route for verification
  route = f"{FRONTEND_URL}/verify/{link}"
  msg = (
    f"<p>Dear {user.name if user.name else 'User'},</p>"
    "<p>Thank you for registering with Toxindex.com! To complete your registration and activate your account, please click on the link below:</p>"
    f"<p><a href='{route}'><b>Activate My Account</b></a></p>"
    "<p>If you did not request this verification, please ignore this email.</p>"
    "<p>Best regards,<br>The Toxindex Team</p>"
  )
  sendgrid.send(user.email, "Activate Your Toxindex.com Account", msg)

@csrf.exempt
@auth_bp.route('/verification/<token>', methods=['GET'])
def api_verification(token):
    try:
        user_id = ds.find('SELECT user_id from user_link WHERE link_token=(%s)', (token,))['user_id']
        ds.execute('UPDATE users SET email_verified = TRUE WHERE user_id = (%s)', (user_id,))
        user = User.get(user_id)
        if user is not None:
            # flask_login.login_user(user)
            return jsonify({'success': True, 'message': 'Email verified!'})
        return jsonify({'error': 'Invalid or expired token.'}), 400
    except Exception:
        return jsonify({'error': 'Invalid or expired token.'}), 400

@csrf.exempt
@auth_bp.route('/login', methods=['POST'])
def api_login():
    # Accept both JSON and form data
    if request.is_json:
        data = request.get_json()
        email = data.get('email')
        password = data.get('password')
    else:
        email = request.form.get('email')
        password = request.form.get('password')

    user = User.get_user(email)
    if not user or not user.validate_password(password):
        return jsonify({'error': random.choice(ERROR_MESSAGES)}), 401
    # Require email verification except for test account
    if email != 'test@test.com' and not getattr(user, 'email_verified', False):
        return jsonify({'error': 'Please verify your email before logging in.'}), 403

    flask_login.login_user(user)
    return jsonify({'success': True, 'user_id': user.user_id, 'email': user.email})

@csrf.exempt
@auth_bp.route('/logout', methods=['POST'])
def api_logout():
    flask_login.logout_user()
    flask.session.clear()  # Explicitly clear the session
    return jsonify({"success": True})

def send_password_reset(email):
  link = generate_validation_link(User.get_user(email))
  route = url_for("reset_password", token=link, _external=True) # Added token parameter and _external=True
  msg = f"Reset your password at this <a href='{route}'>link</a>"
  sendgrid.send(email,"Toxindex.com Password Reset",msg)

@auth_bp.route('/forgot_password', methods=['POST'])
def api_forgot_password():
    data = flask.request.get_json()
    email = data.get('email')
    if not email:
        return flask.jsonify({'success': False, 'error': 'Email is required.'}), 400
    if not User.user_exists(email):
        return flask.jsonify({'success': True})  # Don't reveal if user exists
    send_password_reset(email)
    return flask.jsonify({'success': True})

@auth_bp.route('/reset_password/<token>', methods=['POST'])
def api_reset_password(token):
    data = flask.request.get_json()
    password = data.get('password')
    if not password:
        return flask.jsonify({'success': False, 'error': 'Password is required.'}), 400
    user_id = ds.find('SELECT user_id from user_link WHERE link_token=(%s)', (token,))['user_id']
    ds.execute('UPDATE users set password = (%s) WHERE user_id = (%s)', (password, user_id,))
    if User.get(user_id) is not None:
        flask_login.login_user(User.get(user_id))
    return flask.jsonify({'success': True})

@csrf.exempt
@auth_bp.route('/register', methods=['POST'])
def api_register():
    # Accept both JSON and form data
    if request.is_json:
        data = request.get_json()
        email = data.get('email')
        password = data.get('password')
        password_confirmation = data.get('password_confirmation')
    else:
        email = request.form.get('email')
        password = request.form.get('password')
        password_confirmation = request.form.get('password_confirmation')

    # Validate input
    if not email or not password or not password_confirmation:
        return jsonify({'error': 'All fields are required.'}), 400
    if len(email) < 6 or '@' not in email:
        return jsonify({'error': 'Please enter a valid email address.'}), 400
    if password != password_confirmation:
        return jsonify({'error': 'Passwords do not match.'}), 400
    if User.user_exists(email):
        return jsonify({'error': 'Email already registered.'}), 409

    user = User.create_user(email, password)
    if user is None:
        return jsonify({'error': 'Failed to create user.'}), 500
    validate(user)
    return jsonify({'success': True, 'message': 'Check your email to verify your account.'})
from webserver.model.user import User
from flask_login import LoginManager

login_manager = LoginManager()

@login_manager.user_loader
def load_user(user_id):
    return User.get(user_id)

@login_manager.request_loader
def request_loader(request):
  email = request.form.get('email')
  if User.user_exists(email):
    user = User.get_user(email)
    print(f'got user {user} with email {email}')
    return User.get_user(email)
  else:
    user = User.get_user(email)
    print(f'got user {user} with email {email}')
    return None

def init(app):
  login_manager.init_app(app)
  login_manager.login_view = "/"
  return login_manager

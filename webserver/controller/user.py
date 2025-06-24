from flask import Blueprint, jsonify
from webserver.model.user import User
import flask_login
from webserver.csrf import csrf

user_bp = Blueprint('users', __name__, url_prefix='/api/users')

@csrf.exempt
@user_bp.route('/<user_id>', methods=['GET'])
@flask_login.login_required
def get_user_by_id(user_id):
    user = User.get(user_id)
    if not user:
        return jsonify({'error': 'User not found'}), 404
    return jsonify({
        'user_id': user.user_id,
        'email': user.email,
    }) 
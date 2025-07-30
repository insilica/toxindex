from flask import Blueprint, request, jsonify
import flask_login
from webserver.model import Message
from webserver.csrf import csrf
from webserver.util import is_valid_uuid
from webserver.model.chat_session import ChatSession 
import logging

chat_bp = Blueprint('chat_sessions', __name__, url_prefix='/api/chat_sessions')

@chat_bp.route('/<session_id>/messages', methods=['GET'])
@flask_login.login_required
def get_chat_session_messages(session_id):
    if not is_valid_uuid(session_id):
        return jsonify({'error': 'Invalid session ID'}), 400
    messages = Message.get_messages_by_session(session_id)
    return jsonify({'messages': [m.to_dict() for m in messages]})

@csrf.exempt
@chat_bp.route('/<session_id>', methods=['DELETE'])
@flask_login.login_required
def delete_chat_session(session_id):
    if not is_valid_uuid(session_id):
        return jsonify({'error': 'Invalid session ID'}), 400
    user_id = flask_login.current_user.user_id
    
    ChatSession.delete_session(session_id, user_id)
    return jsonify({'success': True})

@csrf.exempt
@chat_bp.route('/<session_id>', methods=['PATCH'])
@flask_login.login_required
def rename_chat_session(session_id):
    if not is_valid_uuid(session_id):
        return jsonify({'error': 'Invalid session ID'}), 400
    data = request.get_json(force=True)
    new_title = data.get('title')
    if not new_title:
        return jsonify({'success': False, 'error': 'Missing title'}), 400
    ChatSession.update_title(session_id, new_title)
    return jsonify({'success': True})

@chat_bp.route('', methods=['GET'])
@flask_login.login_required
def list_chat_sessions():
    user_id = flask_login.current_user.user_id
    sessions = ChatSession.get_sessions_by_user(user_id)
    return jsonify({'sessions': [s.to_dict() for s in sessions]})

@csrf.exempt
@chat_bp.route('', methods=['POST'])
@flask_login.login_required
def create_chat_session():
    user_id = flask_login.current_user.user_id
    data = request.get_json(force=True) or {}
    title = data.get('title') or 'New chat'
    environment_id = data.get('environment_id')
    session = ChatSession.create_session(environment_id, user_id, title)
    if not session:
        return jsonify({'error': 'Failed to create chat session'}), 500
    return jsonify(session.to_dict()) 
from flask import Blueprint, request, jsonify
import flask_login
from webserver.model import Message, Task
from webserver.csrf import csrf
from webserver.util import is_valid_uuid
from webserver.ai_service import generate_title
from webserver.model.chat_session import ChatSession as CSModel
from workflows.plain_openai_tasks import plain_openai_task
from workflows.probra import probra_task
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
@chat_bp.route('/<session_id>/message', methods=['POST'])
@flask_login.login_required
def send_message_in_session(session_id):
    if not is_valid_uuid(session_id):
        return jsonify({'error': 'Invalid session ID'}), 400
    
    user_id = flask_login.current_user.user_id
    data = request.get_json() or {}
    prompt = data.get('prompt')
    model = data.get('model')
    environment_id = data.get('environment_id')
    if not prompt:
        return jsonify({'error': 'Missing prompt'}), 400
    # Create a new task for this message
    title = generate_title(prompt)
    workflow_id = 1 if model == 'toxindex-rap' else 2 if model == 'toxindex-vanilla' else 1
    task = Task.create_task(
        title=title,
        user_id=user_id,
        workflow_id=workflow_id,
        environment_id=environment_id,
        session_id=session_id,
    )
    # Add user message to session
    Message.create_message(task.task_id, user_id, 'user', prompt, session_id=session_id)
    # If this is the first user message, update the session title
    messages = Message.get_messages_by_session(session_id)
    user_messages = [m for m in messages if m.role == 'user']
    if len(user_messages) == 1:
        CSModel.update_title(session_id, generate_title(prompt))
        logging.info(f"[DEBUG] Updated chat session {session_id} title to: {generate_title(prompt)}")

    if model == 'toxindex-rap':
        celery_task = probra_task.delay({
            'payload': prompt,
            'task_id': task.task_id,
            'user_id': str(user_id),
        })
    elif model == 'toxindex-vanilla':
        celery_task = plain_openai_task.delay({
            'payload': prompt,
            'task_id': task.task_id,
            'user_id': str(user_id),
        })
    else:
        celery_task = probra_task.delay({
            'payload': prompt,
            'task_id': task.task_id,
            'user_id': str(user_id),
        })
    Task.update_celery_task_id(task.task_id, celery_task.id)
    return jsonify({'task_id': task.task_id, 'celery_id': celery_task.id})

@csrf.exempt
@chat_bp.route('/<session_id>', methods=['DELETE'])
@flask_login.login_required
def delete_chat_session(session_id):
    if not is_valid_uuid(session_id):
        return jsonify({'error': 'Invalid session ID'}), 400
    user_id = flask_login.current_user.user_id
    
    CSModel.delete_session(session_id, user_id)
    return jsonify({'success': True})

@csrf.exempt
@chat_bp.route('/<session_id>', methods=['PATCH'])
@flask_login.login_required
def rename_chat_session(session_id):
    if not is_valid_uuid(session_id):
        return jsonify({'error': 'Invalid session ID'}), 400
    user_id = flask_login.current_user.user_id
    data = request.get_json(force=True)
    new_title = data.get('title')
    if not new_title:
        return jsonify({'success': False, 'error': 'Missing title'}), 400
    CSModel.update_title(session_id, new_title)
    return jsonify({'success': True})

@chat_bp.route('', methods=['GET'])
@flask_login.login_required
def list_chat_sessions():
    user_id = flask_login.current_user.user_id
    sessions = CSModel.get_sessions_by_user(user_id)
    return jsonify({'sessions': [s.to_dict() for s in sessions]})

@csrf.exempt
@chat_bp.route('', methods=['POST'])
@flask_login.login_required
def create_chat_session():
    user_id = flask_login.current_user.user_id
    data = request.get_json(force=True) or {}
    title = data.get('title') or 'New chat'
    environment_id = data.get('environment_id')
    session = CSModel.create_session(environment_id, user_id, title)
    if not session:
        return jsonify({'error': 'Failed to create chat session'}), 500
    return jsonify(session.to_dict()) 
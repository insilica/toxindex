from flask import Blueprint, request, jsonify, send_file, abort
import flask_login
from webserver.model import Environment, File, Task
import os, logging, mimetypes, base64
from webserver.csrf import csrf
from webserver.util import is_valid_uuid

env_bp = Blueprint('environments', __name__, url_prefix='/api/environments')

@env_bp.route('', methods=['GET'])
@flask_login.login_required
def list_environments():
    user_id = flask_login.current_user.user_id
    envs = Environment.get_environments_by_user(user_id)
    if not envs:
        Environment.create_environment("Base environment", user_id, description="Your starter workspace")
        envs = Environment.get_environments_by_user(user_id)
    return jsonify({
        "environments": [e.to_dict() for e in envs]
    })

@csrf.exempt
@env_bp.route('', methods=['POST'])
@flask_login.login_required
def create_environment():
    env_data = request.get_json()
    title = env_data.get("title", "New Environment")
    description = env_data.get("description")
    env = Environment.create_environment(
        title, flask_login.current_user.user_id, description
    )
    return jsonify({"environment_id": env.environment_id})

@env_bp.route('/<env_id>', methods=['GET'])
@flask_login.login_required
def get_environment(env_id):
    env = Environment.get_environment(env_id)
    if not env:
        return jsonify({"error": "Environment not found"}), 404
    tasks = Task.get_tasks_by_environment(env_id, flask_login.current_user.user_id)
    files = File.get_files_by_environment(env_id)
    return jsonify({
        "environment": env.to_dict(),
        "tasks": [t.to_dict() for t in tasks],
        "files": [f.to_dict() for f in files]
    })

@csrf.exempt
@env_bp.route('/<env_id>', methods=['DELETE'])
@flask_login.login_required
def delete_environment(env_id):
    if not is_valid_uuid(env_id):
        logging.error(f"Attempted to delete environment with invalid UUID: {env_id}")
        return jsonify({"success": False, "error": "Invalid environment ID"}), 400
    try:
        files = File.get_files_by_environment(env_id)
        for file in files:
            try:
                if file.filepath and os.path.exists(file.filepath):
                    os.remove(file.filepath)
                    logging.info(f"Deleted file from disk: {file.filepath}")
            except Exception as e:
                logging.warning(f"Failed to remove file from disk: {file.filepath}, error: {e}")
        Environment.delete_environment(env_id, flask_login.current_user.user_id)
        return jsonify({"success": True})
    except Exception as e:
        logging.error(f"Failed to delete environment {env_id}: {e}")
        return jsonify({"success": False, "error": str(e)}), 500

@env_bp.route('/<env_id>/files', methods=['GET'])
@flask_login.login_required
def list_files(env_id):
    if not is_valid_uuid(env_id):
        return jsonify({'error': 'Invalid environment ID'}), 400
    files = File.get_files_by_environment(env_id)
    return jsonify({'files': [f.to_dict() for f in files]})

@csrf.exempt
@env_bp.route('/<env_id>/files', methods=['POST'])
@flask_login.login_required
def upload_file(env_id):
    if not is_valid_uuid(env_id):
        return jsonify({'error': 'Invalid environment ID'}), 400
    if 'file' not in request.files:
        return jsonify({'error': 'No file part'}), 400
    file = request.files['file']
    if file.filename == '':
        return jsonify({'error': 'No selected file'}), 400
    from werkzeug.utils import secure_filename
    filename = secure_filename(file.filename)
    upload_dir = os.path.join(os.path.dirname(__file__), '..', '..', 'uploads')
    os.makedirs(upload_dir, exist_ok=True)
    file_path = os.path.join(upload_dir, filename)
    file.save(file_path)
    user_id = flask_login.current_user.user_id
    File.create_file(
        task_id=None,
        user_id=user_id,
        filename=filename,
        filepath=file_path,
        environment_id=env_id
    )
    return jsonify({'success': True, 'filename': filename})

@env_bp.route('/<env_id>/files/<file_id>', methods=['GET'])
@flask_login.login_required
def get_file_info(env_id, file_id):
    if not is_valid_uuid(env_id) or not is_valid_uuid(file_id):
        return jsonify({'error': 'Invalid environment or file ID'}), 400
    file = File.get_file(file_id)
    if not file or str(file.environment_id) != str(env_id):
        return jsonify({'error': 'File not found'}), 404
    return jsonify(file.to_dict())

@csrf.exempt
@env_bp.route('/<env_id>/files/<file_id>', methods=['DELETE'])
@flask_login.login_required
def delete_file(env_id, file_id):
    if not is_valid_uuid(env_id) or not is_valid_uuid(file_id):
        return jsonify({'error': 'Invalid environment or file ID'}), 400
    try:
        file = File.get_file(file_id)
        if not file or str(file.environment_id) != str(env_id):
            return jsonify({'success': False, 'error': 'File not found'}), 404
        try:
            if file.filepath and os.path.exists(file.filepath):
                os.remove(file.filepath)
        except Exception as e:
            logging.warning(f"Failed to remove file from disk: {e}")
        File.delete_file(file_id, flask_login.current_user.user_id)
        return jsonify({'success': True})
    except Exception as e:
        logging.error(f"Failed to delete file {file_id}: {e}")
        return jsonify({'success': False, 'error': str(e)}), 500

# @env_bp.route('/<env_id>/files/<file_id>/inspect', methods=['GET'])
# @flask_login.login_required
# def inspect_file(env_id, file_id):
#     if not is_valid_uuid(env_id) or not is_valid_uuid(file_id):
#         return jsonify({'error': 'Invalid environment or file ID'}), 400
#     file = File.get_file(file_id)
#     if not file or str(file.environment_id) != str(env_id):
#         return jsonify({'error': 'File not found'}), 404
#     if not file.filepath:
#         return jsonify({'error': 'File not found'}), 404
#     if not os.path.exists(file.filepath):
#         return jsonify({'error': 'File not found'}), 404
#     ext = os.path.splitext(file.filename)[1].lower()
#     mimetype, _ = mimetypes.guess_type(file.filename)
#     try:
#         if ext in ['.txt', '.csv', '.json', '.md', '.markdown']:
#             with open(file.filepath, 'r', encoding='utf-8') as f:
#                 content = f.read()
#             if ext == '.json':
#                 import json
#                 try:
#                     parsed = json.loads(content)
#                     return jsonify({'type': 'json', 'content': parsed, 'filename': file.filename, 'mimetype': mimetype})
#                 except Exception as e:
#                     return jsonify({'type': 'text', 'content': content, 'filename': file.filename, 'mimetype': mimetype, 'warning': f'Invalid JSON: {e}'})
#             elif ext in ['.md', '.markdown']:
#                 import markdown
#                 html = markdown.markdown(content, extensions=["extra", "tables"])
#                 return jsonify({'type': 'markdown', 'content': html, 'filename': file.filename, 'mimetype': mimetype})
#             elif ext == '.csv':
#                 import csv
#                 import io
#                 reader = csv.reader(io.StringIO(content))
#                 rows = list(reader)
#                 return jsonify({'type': 'csv', 'content': rows, 'filename': file.filename, 'mimetype': mimetype})
#             else:
#                 return jsonify({'type': 'text', 'content': content, 'filename': file.filename, 'mimetype': mimetype})
#         elif ext in ['.xlsx', '.xls']:
#             try:
#                 import pandas as pd
#                 df = pd.read_excel(file.filepath)
#                 preview = df.head(100).to_dict(orient='records')
#                 columns = list(df.columns)
#                 return jsonify({'type': 'xlsx', 'content': preview, 'columns': columns, 'filename': file.filename, 'mimetype': mimetype})
#             except Exception as e:
#                 return jsonify({'error': f'Failed to parse Excel: {e}', 'type': 'error', 'filename': file.filename, 'mimetype': mimetype}), 400
#         elif ext in ['.png', '.jpg', '.jpeg', '.gif']:
#             with open(file.filepath, 'rb') as f:
#                 data = f.read()
#             b64 = base64.b64encode(data).decode('utf-8')
#             data_url = f"data:{mimetype};base64,{b64}"
#             return jsonify({'type': 'image', 'content': data_url, 'filename': file.filename, 'mimetype': mimetype})
#         else:
#             return jsonify({'error': 'Preview not supported for this file type', 'type': 'unsupported', 'filename': file.filename, 'mimetype': mimetype}), 415
#     except Exception as e:
#         logging.error(f"[inspect_file] Exception for file_id={file_id}: {e}")
#         return jsonify({'error': str(e), 'type': 'error', 'filename': file.filename, 'mimetype': mimetype}), 500

@env_bp.route('/<env_id>/files/<file_id>/download', methods=['GET'])
@flask_login.login_required
def download_file(env_id, file_id):
    if not is_valid_uuid(env_id) or not is_valid_uuid(file_id):
        return jsonify({'error': 'Invalid environment or file ID'}), 400
    file = File.get_file(file_id)
    if not file or not file.filepath or not os.path.exists(file.filepath) or str(file.environment_id) != str(env_id):
        return abort(404)
    return send_file(file.filepath, as_attachment=True, download_name=file.filename)

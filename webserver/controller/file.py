from flask import Blueprint, jsonify, send_file, abort
import flask_login
from webserver.model.file import File
import os, mimetypes, base64, logging

import json
import pandas as pd
import csv
import io

file_bp = Blueprint('files', __name__, url_prefix='/api/files')

@file_bp.route('', methods=['GET'])
@flask_login.login_required
def get_user_files():
    user_id = flask_login.current_user.user_id
    files = File.get_files_by_user(user_id)
    return jsonify({'files': [f.to_dict() for f in files]})

@file_bp.route('/<file_id>/inspect', methods=['GET'])
@flask_login.login_required
def inspect_file(file_id):
    file = File.get_file(file_id)
    if not file:
        return jsonify({'error': 'File not found'}), 404
    # Optionally, check user ownership here if needed
    if not file.filepath:
        return jsonify({'error': 'File not found'}), 404
    if not os.path.exists(file.filepath):
        return jsonify({'error': 'File not found'}), 404
    ext = os.path.splitext(file.filename)[1].lower()
    mimetype, _ = mimetypes.guess_type(file.filename)

    # supported file types:
    # .txt, .csv, .json, .md, .markdown, .xlsx, .xls, .png, .jpg, .jpeg, .gif
    try:
        if ext in ['.txt', '.csv', '.json', '.md', '.markdown']:
            with open(file.filepath, 'r', encoding='utf-8') as f:
                content = f.read()
            if ext == '.json':
                try:
                    parsed = json.loads(content)
                    return jsonify({'type': 'json', 'content': parsed, 'filename': file.filename, 'mimetype': mimetype})
                except Exception as e:
                    return jsonify({'type': 'text', 'content': content, 'filename': file.filename, 'mimetype': mimetype, 'warning': f'Invalid JSON: {e}'})
            elif ext in ['.md', '.markdown']:
                # Return raw markdown, not HTML
                return jsonify({'type': 'markdown', 'content': content, 'filename': file.filename, 'mimetype': mimetype})
            elif ext == '.csv':
                reader = csv.reader(io.StringIO(content))
                rows = list(reader)
                return jsonify({'type': 'csv', 'content': rows, 'filename': file.filename, 'mimetype': mimetype})
            else:
                return jsonify({'type': 'text', 'content': content, 'filename': file.filename, 'mimetype': mimetype})
        elif ext in ['.xlsx', '.xls']:
            try:
                df = pd.read_excel(file.filepath)
                preview = df.head(100).to_dict(orient='records')
                columns = list(df.columns)
                return jsonify({'type': 'xlsx', 'content': preview, 'columns': columns, 'filename': file.filename, 'mimetype': mimetype})
            except Exception as e:
                return jsonify({'error': f'Failed to parse Excel: {e}', 'type': 'error', 'filename': file.filename, 'mimetype': mimetype}), 400
        elif ext in ['.png', '.jpg', '.jpeg', '.gif']:
            with open(file.filepath, 'rb') as f:
                data = f.read()
            b64 = base64.b64encode(data).decode('utf-8')
            data_url = f"data:{mimetype};base64,{b64}"
            return jsonify({'type': 'image', 'content': data_url, 'filename': file.filename, 'mimetype': mimetype})
        else:
            return jsonify({'error': 'Preview not supported for this file type', 'type': 'unsupported', 'filename': file.filename, 'mimetype': mimetype}), 415
    except Exception as e:
        logging.error(f"[inspect_file] Exception for file_id={file_id}: {e}")
        return jsonify({'error': str(e), 'type': 'error', 'filename': file.filename, 'mimetype': mimetype}), 500

@file_bp.route('/<file_id>/download', methods=['GET'])
@flask_login.login_required
def download_file(file_id):
    file = File.get_file(file_id)
    if not file or not file.filepath or not os.path.exists(file.filepath):
        return abort(404)
    return send_file(file.filepath, as_attachment=True, download_name=file.filename) 
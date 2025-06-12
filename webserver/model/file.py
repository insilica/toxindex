import os
import logging
import datetime
from typing import Optional

import webserver.datastore as ds
from webserver.storage import S3FileStorage

class File:

    def __init__(self, file_id, task_id, user_id, filename, filepath, s3_url, created_at, environment_id=None):
        self.file_id = file_id
        self.task_id = task_id
        self.user_id = user_id
        self.filename = filename
        self.filepath = filepath
        self.s3_url = s3_url
        self.created_at = created_at
        self.environment_id = environment_id

    def to_dict(self):
        return {
            'file_id': self.file_id,
            'task_id': self.task_id,
            'environment_id': self.environment_id,
            'user_id': self.user_id,
            'filename': self.filename,
            'filepath': self.filepath,
            's3_url': self.s3_url,
            'created_at': self.created_at.strftime('%Y-%m-%d %H:%M:%S') if self.created_at else None
        }

    @staticmethod
    def from_row(row):
        return File(
            file_id=row['file_id'],
            task_id=row.get('task_id'),
            user_id=row['user_id'],
            filename=row['filename'],
            filepath=row['filepath'],
            s3_url=row['s3_url'],
            created_at=row['created_at'],
            environment_id=row.get('environment_id'),
        )

    @staticmethod
    def create_file(task_id, user_id, filename, filepath, s3_url, environment_id=None):
        # Suppress S3 URL in logs
        logging.info(f"[File.create_file] Storing file for task_id={task_id}, environment_id={environment_id}, filename={filename}, filepath={filepath}, s3_url=SUPPRESSED")
        # Duplicate check: does a file with this task_id/environment_id and filename already exist?
        if task_id:
            existing = ds.find("SELECT 1 FROM files WHERE task_id = %s AND filename = %s", (task_id, filename))
        elif environment_id:
            existing = ds.find("SELECT 1 FROM files WHERE environment_id = %s AND filename = %s", (environment_id, filename))
        else:
            existing = None
        if existing:
            logging.warning(f"[File.create_file] Duplicate file for task_id={task_id}, environment_id={environment_id}, filename={filename} -- skipping insert.")
            return
        params = (task_id, environment_id, user_id, filename, filepath, s3_url)
        try:
            ds.execute(
                "INSERT INTO files (task_id, environment_id, user_id, filename, filepath, s3_url) VALUES (%s, %s, %s, %s, %s, %s)",
                params
            )
            logging.info(f"[File.create_file] Successfully inserted file for task_id={task_id}, environment_id={environment_id}, filename={filename}")
        except Exception as e:
            logging.error(f"[File.create_file] Failed to insert file for task_id={task_id}, environment_id={environment_id}, filename={filename}: {e}")

    @staticmethod
    def upload_and_create(
        task_id: int,
        user_id: str,
        local_path: str,
        storage: S3FileStorage,
        key: Optional[str] = None,
        content_type: Optional[str] = None,
    ) -> str:
        """Upload a file to S3 and create a database record."""
        filename = os.path.basename(local_path)
        s3_key = key or f"{user_id}/{filename}"
        storage_key = storage.upload_file(local_path, s3_key, content_type)
        download_url = storage.generate_download_url(storage_key)
        File.create_file(task_id, user_id, filename, storage_key, download_url)
        return download_url

    @staticmethod
    def get_files(task_id):
        import os
        try:
            task_id_int = int(task_id)
            logging.info(f"DB ENV (get_files): PGHOST={os.getenv('PGHOST')}, PGPORT={os.getenv('PGPORT')}, PGDATABASE={os.getenv('PGDATABASE')}, PGUSER={os.getenv('PGUSER')}, PGPASSWORD={os.getenv('PGPASSWORD')}")
            logging.info(f"Calling get_files with task_id={task_id_int} (type: {type(task_id_int)})")
            rows = ds.find_all(
                "SELECT * FROM files WHERE task_id = %s ORDER BY created_at ASC",
                (task_id_int,)
            )
            logging.info(f"Raw rows returned from DB for task_id={task_id_int}: {rows}")
            if not rows:
                logging.warning(f"No files found in DB for task_id={task_id_int}")
            else:
                logging.info(f"Found {len(rows)} files in DB for task_id={task_id_int}")
            return [File.from_row(row) for row in rows]
        except Exception as e:
            logging.error(f"Error retrieving files for task_id={task_id}: {e}")
            return []

    @staticmethod  
    def process_event(task, event_data):
        user_id = event_data.get("user_id")
        filename = event_data.get("filename")
        filepath = event_data.get("filepath")
        s3_url = event_data.get("s3_url")
        logging.info(f"[File.process_event] Processing event for task_id={task.task_id}, filename={filename}, filepath={filepath}, s3_url=SUPPRESSED")
        if user_id and filename and filepath and s3_url:
            File.create_file(task.task_id, user_id, filename, filepath, s3_url)
            logging.info(f"[File.process_event] Stored file {filename} for task_id={task.task_id}")
        else:
            missing = []
            if not user_id:
                missing.append("user_id")
            if not filename:
                missing.append("filename") 
            if not filepath:
                missing.append("filepath")
            if not s3_url:
                missing.append("s3_url")
            logging.warning(f"[File.process_event] Malformed task_file event received - missing fields: {', '.join(missing)}")
            logging.debug(f"[File.process_event] Received payload: {event_data}")

    # ------------------------------------------------------------------
    def read_text(self):
        """Return file contents if the file is available on disk."""
        try:
            with open(self.filepath, "r", encoding="utf-8") as f:
                return f.read()
        except Exception as e:
            logging.warning(f"Could not read file {self.filepath}: {e}")
            return ""

    def markdown_to_html(self):
        text = self.read_text()
        if not text:
            return ""
        import markdown
        return markdown.markdown(text, extensions=["extra", "tables"])

    @staticmethod
    def get_files_by_environment(environment_id):
        rows = ds.find_all(
            "SELECT * FROM files WHERE environment_id = %s ORDER BY created_at DESC",
            (environment_id,)
        )
        return [File.from_row(row) for row in rows]

    @staticmethod
    def get_file(file_id):
        row = ds.find("SELECT * FROM files WHERE file_id = %s", (file_id,))
        return File.from_row(row) if row else None

    @staticmethod
    def delete_file(file_id, user_id=None):
        params = (file_id,)
        query = "DELETE FROM files WHERE file_id = %s"
        if user_id is not None:
            query += " AND user_id = %s"
            params = (file_id, user_id)
        ds.execute(query, params)


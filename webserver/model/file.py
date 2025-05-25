import os
import logging
import datetime
from typing import Optional

import webserver.datastore as ds
from webserver.storage import S3FileStorage

class File:

    def __init__(self, file_id, task_id, user_id, filename, filepath, s3_url, created_at):
        self.file_id = file_id
        self.task_id = task_id
        self.user_id = user_id
        self.filename = filename
        self.filepath = filepath
        self.s3_url = s3_url
        self.created_at = created_at

    def to_dict(self):
        return {
            'file_id': self.file_id,
            'task_id': self.task_id,
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
            task_id=row['task_id'],
            user_id=row['user_id'],
            filename=row['filename'],
            filepath=row['filepath'],
            s3_url=row['s3_url'],
            created_at=row['created_at']
        )

    @staticmethod
    def create_file(task_id, user_id, filename, filepath, s3_url):
        logging.info(f"Storing file for task_id={task_id}, filename={filename}")
        params = (task_id, user_id, filename, filepath, s3_url)
        ds.execute(
            "INSERT INTO files (task_id, user_id, filename, filepath, s3_url) VALUES (%s, %s, %s, %s, %s)",
            params
        )

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
        rows = ds.find_all(
            "SELECT * FROM files WHERE task_id = %s ORDER BY created_at ASC",
            (task_id,)
        )
        return [File.from_row(row) for row in rows]

    @staticmethod
    def process_event(task, event_data):
        payload = event_data.get("data", {})
        user_id = payload.get("user_id")
        filename = payload.get("filename")
        filepath = payload.get("filepath")
        s3_url = payload.get("s3_url")

        if user_id and filename and filepath and s3_url:
            File.create_file(task.task_id, user_id, filename, filepath, s3_url)
            logging.info(f"Stored file {filename} for task_id={task.task_id}")
        else:
            logging.warning("Malformed task_file event received")

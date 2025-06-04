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
        logging.info(f"DB ENV (create_file): PGHOST={os.getenv('PGHOST')}, PGPORT={os.getenv('PGPORT')}, PGDATABASE={os.getenv('PGDATABASE')}, PGUSER={os.getenv('PGUSER')}, PGPASSWORD={os.getenv('PGPASSWORD')}")
        logging.info(f"Storing file for task_id={task_id}, filename={filename}")
        params = (task_id, user_id, filename, filepath, s3_url)
        try:
            ds.execute(
                "INSERT INTO files (task_id, user_id, filename, filepath, s3_url) VALUES (%s, %s, %s, %s, %s)",
                params
            )
            logging.info(f"Successfully inserted file for task_id={task_id}, filename={filename}")
        except Exception as e:
            logging.error(f"Failed to insert file for task_id={task_id}, filename={filename}: {e}")

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

        if user_id and filename and filepath and s3_url:
            File.create_file(task.task_id, user_id, filename, filepath, s3_url)
            logging.info(f"Stored file {filename} for task_id={task.task_id}")
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
            logging.warning(f"Malformed task_file event received - missing fields: {', '.join(missing)}")
            logging.debug(f"Received payload: {event_data}")

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
        """Very small markdown to HTML converter supporting headings, lists,
        bold and italic text."""
        text = self.read_text()
        if not text:
            return ""

        import html
        import re

        def repl_heading(match):
            hashes = match.group(1)
            content = match.group(2)
            level = len(hashes)
            return f"<h{level}>" + html.escape(content) + f"</h{level}>"

        # Headings
        text = re.sub(r"^(#{1,6})\s*(.+)$", repl_heading, text, flags=re.MULTILINE)

        # Bold and italic
        text = re.sub(r"\*\*(.+?)\*\*", r"<strong>\1</strong>", text)
        text = re.sub(r"\*(.+?)\*", r"<em>\1</em>", text)

        # Lists
        lines = text.split("\n")
        html_lines = []
        in_list = False
        for line in lines:
            if re.match(r"^[-*]\s+", line):
                if not in_list:
                    html_lines.append("<ul>")
                    in_list = True
                item = re.sub(r"^[-*]\s+", "", line)
                html_lines.append(f"<li>{html.escape(item)}</li>")
            else:
                if in_list:
                    html_lines.append("</ul>")
                    in_list = False
                if line.strip() == "":
                    html_lines.append("<br>")
                else:
                    html_lines.append(f"<p>{html.escape(line)}</p>")
        if in_list:
            html_lines.append("</ul>")

        return "\n".join(html_lines)


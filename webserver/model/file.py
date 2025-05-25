import webserver.datastore as ds
import logging
import datetime

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


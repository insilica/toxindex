import webserver.datastore as ds
import logging
import json
import datetime
from pydantic import BaseModel

class MessageSchema(BaseModel):
    role: str
    content: str

class Message():

    def __init__(self, message_id, task_id, user_id, role, content, created_at):
        self.message_id = message_id
        self.task_id = task_id
        self.user_id = user_id
        self.role = role
        self.content = content
        self.created_at = created_at

    def to_dict(self):
        return {
            'message_id': str(self.message_id),
            'task_id': self.task_id,
            'user_id': str(self.user_id),
            'role': self.role,
            'content': self.content,
            'created_at': self.created_at.strftime('%Y-%m-%d %H:%M:%S') if self.created_at else None
        }

    
    def to_json(self):
        return json.dumps(self.to_dict())
    
    @staticmethod
    def from_row(row):
        return Message(
            message_id=row['message_id'],
            task_id=row['task_id'],
            user_id=row['user_id'],
            role=row['role'],
            content=row['content'],
            created_at=row['created_at']
        )

    @staticmethod
    def create_message(task_id, user_id, role, content):
        # Duplicate check: does a message with this task_id, role, and content already exist?
        existing = ds.find("SELECT 1 FROM messages WHERE task_id = %s AND role = %s AND content = %s", (task_id, role, content))
        if existing:
            logging.warning(f"[Message.create_message] Duplicate message for task_id={task_id}, role={role} -- skipping insert.")
            return
        logging.info(f"[Message.create_message] Storing message for task_id={task_id}, user_id={user_id}, role={role}, content={content}")
        params = (task_id, user_id, role, content)
        ds.execute(
            "INSERT INTO messages (task_id, user_id, role, content) VALUES (%s, %s, %s, %s)",
            params
        )
        logging.info(f"[Message.create_message] Message stored for task_id={task_id}")

    @staticmethod
    def get_messages(task_id):
        rows = ds.find_all(
            "SELECT * FROM messages WHERE task_id = %s ORDER BY created_at ASC",
            (task_id,)
        )
        logging.info(f"Retrieved {len(rows)} messages for task_id={task_id}")
        return [Message.from_row(row) for row in rows]

    @staticmethod
    def process_event(task, event_data):
        logging.info(f"[Message.process_event] Processing event for task_id={task.task_id}, role={event_data.get('role')}, content={event_data.get('content')}")
        role = event_data.get("role", "assistant")
        content = event_data.get("content")
        if content and role:
            Message.create_message(task.task_id, None, role, content)
            logging.info(f"[Message.process_event] Stored message for task_id={task.task_id} from role={role}")
        else:
            logging.warning(f"[Message.process_event] Malformed task_message event received data: {event_data}")
            logging.warning(f"[Message.process_event] Role: {role}")
            logging.warning(f"[Message.process_event] Content: {content}")

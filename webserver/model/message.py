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
        logging.info(f"Storing message for task_id={task_id}, user_id={user_id}, role={role}")
        params = (task_id, user_id, role, content)
        ds.execute(
            "INSERT INTO messages (task_id, user_id, role, content) VALUES (%s, %s, %s, %s)",
            params
        )

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
        role = event_data.get("role", "assistant")
        content = event_data.get("content")

        if content and role:
            Message.create_message(task.task_id, None, role, content)
            logging.info(f"Stored message for task_id={task.task_id} from role={role}")
        else:
            logging.warning(f"Malformed task_message event received data: {event_data}")
            logging.warning(f"Role: {role}")
            logging.warning(f"Content: {content}")

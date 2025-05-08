import webserver.datastore as ds
import logging
import datetime

class Conversation:
    def __init__(self, conversation_id, title, user_id, created_at=None):
        self.conversation_id = conversation_id
        self.title = title
        self.user_id = user_id
        self.created_at = created_at

    def to_dict(self):
        return {
            'conversation_id': self.conversation_id,
            'title': self.title,
            'user_id': self.user_id,
            'created_at': self.created_at.strftime('%Y-%m-%d %H:%M:%S') if self.created_at else None
        }

    @staticmethod
    def from_row(row):
        return Conversation(row['conversation_id'], row['title'], row['user_id'], row['created_at'])

    @staticmethod
    def create_conversation(title, user_id):
        params = (title, user_id)
        ds.execute("INSERT INTO conversations (title, user_id) VALUES (%s, %s)", params)
        
        # Fetch and return the newly created conversation
        res = ds.find("SELECT * FROM conversations WHERE title = %s AND user_id = %s ORDER BY created_at DESC LIMIT 1", 
                     (title, user_id))
        return Conversation.from_row(res) if res else None

    @staticmethod
    def get_conversations_by_user(user_id):
        rows = ds.find_all("SELECT * FROM conversations WHERE user_id = %s ORDER BY created_at DESC", (user_id,))
        return [Conversation.from_row(row) for row in rows]

    @staticmethod
    def add_message(conversation_id, user_id, role, content):
        params = (conversation_id, user_id, role, content)
        ds.execute("INSERT INTO messages (conversation_id, user_id, role, content) VALUES (%s, %s, %s, %s)", params)
        return True

    @staticmethod
    def get_messages(conversation_id, user_id):
        # Check if the conversation belongs to the user
        conv = ds.find("SELECT * FROM conversations WHERE conversation_id = %s AND user_id = %s", 
                      (conversation_id, user_id))
        if not conv:
            return []
        
        rows = ds.find_all("SELECT * FROM messages WHERE conversation_id = %s ORDER BY created_at ASC", 
                          (conversation_id,))
        return [{"role": row['role'], "content": row['content'], "created_at": row['created_at'].strftime('%Y-%m-%d %H:%M:%S')} 
                for row in rows]
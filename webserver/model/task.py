import webserver.datastore as ds
import logging
import datetime

class Task:

    def __init__(self, task_id, title, user_id, workflow_id, celery_task_id=None, description=None, created_at=None):
        self.task_id = task_id
        self.title = title
        self.user_id = user_id
        self.workflow_id = workflow_id
        self.celery_task_id = celery_task_id
        self.description = description
        self.created_at = created_at

    def to_dict(self):
        return {
            'task_id': self.task_id,
            'title': self.title,
            'user_id': self.user_id,
            'workflow_id': self.workflow_id,
            'description': self.description,
            'created_at': self.created_at.strftime('%Y-%m-%d %H:%M:%S') if self.created_at else None
        }

    @staticmethod
    def from_row(row):
        return Task(
            task_id=row['task_id'],
            title=row['title'],
            user_id=row['user_id'],
            workflow_id=row['workflow_id'],
            celery_task_id=row['celery_task_id'],
            description=row.get('description'),
            created_at=row['created_at']
        )

    @staticmethod
    def create_task(title, user_id, workflow_id, celery_task_id=None, description=None):
        logging.info(f"Creating new task with title='{title}' for user_id={user_id}")
        params = (title, user_id, celery_task_id, workflow_id, description)
        ds.execute("INSERT INTO tasks (title, user_id, celery_task_id, workflow_id, description) VALUES (%s, %s, %s, %s, %s)", params)
        
        # Fetch and return the newly created task
        res = ds.find("SELECT * FROM tasks WHERE title = %s AND user_id = %s ORDER BY created_at DESC LIMIT 1",
                     (title, user_id))
        if res:
            logging.info(f"Successfully created task with id={res['task_id']}")
            return Task.from_row(res)
        else:
            logging.error(f"Failed to create task for user_id={user_id}")
            return None

    @staticmethod
    def get_tasks_by_user(user_id):
        rows = ds.find_all("SELECT * FROM tasks WHERE user_id = %s ORDER BY created_at DESC", (user_id,))
        return [Task.from_row(row) for row in rows]
    
    @staticmethod
    def get_tasks_by_celery_task_id(celery_task_id):
        res = ds.find("SELECT * FROM tasks WHERE celery_task_id = %s", (celery_task_id,))
        return Task.from_row(res) if res else None
    
    @staticmethod
    def get_task(task_id):
        res = ds.find("SELECT * FROM tasks WHERE task_id = %s", (task_id,))
        return Task.from_row(res) if res else None

    @staticmethod
    def add_message(task_id, user_id, role, content):
        params = (task_id, user_id, role, content)
        ds.execute("INSERT INTO messages (task_id, user_id, role, content) VALUES (%s, %s, %s, %s)", params)
        return True
    
    @staticmethod
    def update_celery_task_id(task_id, celery_task_id):
        params = (celery_task_id, task_id)
        ds.execute("UPDATE tasks SET celery_task_id = %s WHERE task_id = %s", params)
        return True
    
    @staticmethod
    def get_messages(task_id, user_id):
        # Check if the task belongs to the user
        task = ds.find("SELECT * FROM tasks WHERE task_id = %s AND user_id = %s",
                      (task_id, user_id))
        if not task:
            return []

        rows = ds.find_all("SELECT * FROM messages WHERE task_id = %s ORDER BY created_at ASC",
                          (task_id,))
        return [{"role": row['role'], "content": row['content'], "created_at": row['created_at'].strftime('%Y-%m-%d %H:%M:%S')}
                for row in rows]

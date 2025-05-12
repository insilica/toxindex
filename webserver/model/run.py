import webserver.datastore as ds
import logging
import datetime

class Run:
    def __init__(self, run_id, title, user_id, workflow_id, description=None, created_at=None):
        self.run_id = run_id
        self.title = title
        self.user_id = user_id
        self.workflow_id = workflow_id
        self.description = description
        self.created_at = created_at

    def to_dict(self):
        return {
            'run_id': self.run_id,
            'title': self.title,
            'user_id': self.user_id,
            'workflow_id': self.workflow_id,
            'description': self.description,
            'created_at': self.created_at.strftime('%Y-%m-%d %H:%M:%S') if self.created_at else None
        }

    @staticmethod
    def from_row(row):
        return Run(
            run_id=row['run_id'], 
            title=row['title'], 
            user_id=row['user_id'],
            workflow_id=row['workflow_id'],
            description=row.get('description'),
            created_at=row['created_at']
        )

    @staticmethod
    def create_run(title, user_id, workflow_id, description=None):
        logging.info(f"Creating new run with title='{title}' for user_id={user_id}")
        params = (title, user_id, workflow_id, description)
        ds.execute("INSERT INTO runs (title, user_id, workflow_id, description) VALUES (%s, %s, %s, %s)", params)
        
        # Fetch and return the newly created run
        res = ds.find("SELECT * FROM runs WHERE title = %s AND user_id = %s ORDER BY created_at DESC LIMIT 1", 
                     (title, user_id))
        if res:
            logging.info(f"Successfully created run with id={res['run_id']}")
            return Run.from_row(res)
        else:
            logging.error(f"Failed to create run for user_id={user_id}")
            return None

    @staticmethod
    def get_runs_by_user(user_id):
        rows = ds.find_all("SELECT * FROM runs WHERE user_id = %s ORDER BY created_at DESC", (user_id,))
        return [Run.from_row(row) for row in rows]
    
    @staticmethod
    def get_run(run_id, user_id):
        res = ds.find("SELECT * FROM runs WHERE run_id = %s AND user_id = %s", (run_id, user_id))
        return Run.from_row(res) if res else None

    @staticmethod
    def add_message(run_id, user_id, role, content):
        params = (run_id, user_id, role, content)
        ds.execute("INSERT INTO messages (run_id, user_id, role, content) VALUES (%s, %s, %s, %s)", params)
        return True

    @staticmethod
    def get_messages(run_id, user_id):
        # Check if the run belongs to the user
        run = ds.find("SELECT * FROM runs WHERE run_id = %s AND user_id = %s", 
                      (run_id, user_id))
        if not run:
            return []
        
        rows = ds.find_all("SELECT * FROM messages WHERE run_id = %s ORDER BY created_at ASC", 
                          (run_id,))
        return [{"role": row['role'], "content": row['content'], "created_at": row['created_at'].strftime('%Y-%m-%d %H:%M:%S')} 
                for row in rows]
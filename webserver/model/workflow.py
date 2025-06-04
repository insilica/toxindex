import webserver.datastore as ds
import logging
import datetime
import json

class Workflow:
    def __init__(self, workflow_id, title, user_id, description=None, initial_prompt=None, created_at=None):
        self.workflow_id = workflow_id
        self.title = title
        self.user_id = user_id
        self.description = description
        self.initial_prompt = initial_prompt
        self.created_at = created_at

    def to_dict(self):
        return {
            'workflow_id': self.workflow_id,
            'title': self.title,
            'user_id': self.user_id,
            'description': self.description,
            'initial_prompt': self.initial_prompt,
            'created_at': self.created_at.strftime('%Y-%m-%d %H:%M:%S') if self.created_at else None
        }

    @staticmethod
    def from_row(row):
        return Workflow(
            workflow_id=row['workflow_id'], 
            title=row['title'], 
            user_id=row['user_id'],
            description=row.get('description'),
            initial_prompt=row.get('initial_prompt'),
            created_at=row['created_at']
        )

    @staticmethod
    def create_workflow(title, user_id=None, description=None, initial_prompt=None):
        params = (title, user_id, description, initial_prompt)
        ds.execute("INSERT INTO workflows (title, user_id, description, initial_prompt) VALUES (%s, %s, %s, %s)", params)
        logging.info(f"created workflow {title} for user {user_id}")
        res = ds.find("SELECT * FROM workflows WHERE title = %s ORDER BY created_at DESC LIMIT 1", (title,))
        return Workflow.from_row(res) if res else None
    
    @staticmethod
    def update_workflow(workflow_id, user_id=None,title=None, description=None, initial_prompt=None):
        params = (title, description, initial_prompt, workflow_id)
        ds.execute("UPDATE workflows SET title = %s, description = %s, initial_prompt = %s WHERE workflow_id = %s", params)
        logging.info(f"updated workflow {workflow_id}")
        return Workflow.get_workflow(workflow_id)
    
    @staticmethod
    def get_workflow(workflow_id):
        res = ds.find("SELECT * FROM workflows WHERE workflow_id = %s", (workflow_id,))
        return Workflow.from_row(res) if res else None

    @staticmethod
    def get_workflows_by_user(user_id):
        rows = ds.find_all("SELECT * FROM workflows WHERE user_id = %s ORDER BY created_at DESC", (user_id,))
        return [Workflow.from_row(row) for row in rows]

    @staticmethod
    def add_message(workflow_id, user_id, role, content):
        params = (workflow_id, user_id, role, content)
        ds.execute("INSERT INTO messages (workflow_id, user_id, role, content) VALUES (%s, %s, %s, %s)", params)
        return True

    @staticmethod
    def get_messages(workflow_id, user_id):
        # Check if the workflow belongs to the user
        workflow = ds.find("SELECT * FROM workflows WHERE workflow_id = %s AND user_id = %s", 
                      (workflow_id, user_id))
        if not workflow:
            return []
        
        rows = ds.find_all("SELECT * FROM messages WHERE workflow_id = %s ORDER BY created_at ASC", 
                          (workflow_id,))
        return [{"role": row['role'], "content": row['content'], "created_at": row['created_at'].strftime('%Y-%m-%d %H:%M:%S')} 
                for row in rows]
    
    @staticmethod
    def load_default_workflows():
        with open('resources/default_workflows.json', 'r') as f:
            for w in json.load(f)['workflows']:
                if not ds.find("SELECT 1 FROM workflows WHERE workflow_id = %s", (w['workflow_id'],)):
                    del w['workflow_id']
                    Workflow.create_workflow(**w)
                else:
                    Workflow.update_workflow(**w)
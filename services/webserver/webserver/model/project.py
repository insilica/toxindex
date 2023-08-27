import webserver.datastore as ds
import logging

class Project:
    def __init__(self, project_id, name, description, creator_id, creation_date=None):
        self.project_id = project_id
        self.name = name
        self.description = description
        self.creator_id = creator_id  # User who created the project
        self.creation_date = creation_date  # Date the project was created

    def to_dict(self):
        return {
            'project_id': self.project_id,
            'name': self.name,
            'description': self.description,
            'creator_id': self.creator_id,
            'creation_date': self.creation_date.strftime('%Y-%m-%d %H:%M:%S') if self.creation_date else None
        }

    @staticmethod
    def from_row(row):
        return Project(row['project_id'], row['name'], row['description'], row['user_id'], row['creation_date'])

    @staticmethod
    def get(project_id):
        res = ds.find("SELECT * from projects where project_id = %s", (project_id,))
        return Project.from_row(res) if res else None

    @staticmethod
    def create_project(name, description, creator_id):
        params = (name, description, creator_id)
        ds.execute("INSERT INTO projects (name, description, user_id) values (%s, %s, %s)", params)

        # Fetch and return the newly created project
        # We're assuming that 'name' is unique, but ideally, you should use the RETURNING clause to get the ID directly.
        res = ds.find("SELECT * from projects where name = %s AND user_id = %s", (name, creator_id))
        return Project.from_row(res) if res else None

    @staticmethod
    def get_projects_by_creator(creator_id):
        rows = ds.find_all("SELECT * from projects where user_id = %s", (creator_id,))
        return [Project.from_row(row) for row in rows]

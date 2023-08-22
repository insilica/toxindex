import webserver.datastore as ds

import logging

class Project:
  def __init__(self, project_id, name, description, creator_id):
    self.project_id = project_id
    self.name = name
    self.description = description
    self.creator_id = creator_id # User who created the project
  
  @staticmethod
  def from_row(row):
    return Project(row['project_id'], row['name'], row['description'], row['creator_id'])

  @staticmethod
  def get(project_id):
    res = ds.find("SELECT * from projects where project_id = (%s)",(project_id,))
    return Project.from_row(res) if res is not None else None

  @staticmethod
  def create_project(name, description, creator_id):
    params = (name, description, creator_id)
    ds.execute("INSERT INTO projects (name,description,creator_id) values (%s,%s,%s)", params)
    # Return the created project, you may want to retrieve its ID to do so

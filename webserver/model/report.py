import webserver.datastore as ds

import logging

class Report:
  def __init__(self, report_id, s3_reference, user_id, title, description, status):
    self.report_id = report_id
    self.s3_reference = s3_reference
    self.user_id = user_id
    self.title = title
    self.description = description
    self.status = status # You may want to map this to an actual enum or choices
  
  @staticmethod
  def from_row(row):
    return Report(row['report_id'], row['s3_reference'], row['user_id'], row['title'], row['description'], row['status'])

  @staticmethod
  def get(report_id):
    res = ds.find("SELECT * from reports where report_id = (%s)",(report_id,))
    return Report.from_row(res) if res is not None else None

  @staticmethod
  def create_report(s3_reference, user_id, title, description, status='PENDING'):
    params = (s3_reference, user_id, title, description, status)
    ds.execute("INSERT INTO reports (s3_reference,user_id,title,description,status) values (%s,%s,%s,%s,%s)", params)
    # Return the created report, you may want to retrieve its ID to do so

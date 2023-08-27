import boto3
import hashlib
import os
import report.datastore as ds
import logging

class Report:

    def __init__(self, report_id, s3_reference, user_id, title, description, created_at, updated_at, status):
        self.report_id = report_id
        self.s3_reference = s3_reference
        self.user_id = user_id
        self.title = title
        self.description = description
        self.created_at = created_at
        self.updated_at = updated_at
        self.status = status
        
    def to_dict(self):
        return {
            'report_id': self.report_id,
            's3_reference': self.s3_reference,
            'user_id': self.user_id,
            'title': self.title,
            'description': self.description,
            'created_at': self.created_at.strftime('%Y-%m-%d %H:%M:%S') if self.created_at else None,  # Convert datetime to string
            'updated_at': self.updated_at.strftime('%Y-%m-%d %H:%M:%S') if self.updated_at else None,  # Convert datetime to string
            'status': self.status
        }

    @staticmethod
    def from_row(row):
        logging.info(f"Creating report from row: {row}")
        return Report(row['report_id'], row['s3_reference'], row['user_id'], row['title'], row['description'], row['created_at'], row['updated_at'], row['status'])

    @staticmethod
    def get(report_id):
        res = ds.find("SELECT * from reports where report_id = %s", (report_id,))
        return Report.from_row(res) if res else None

    @staticmethod
    def upload_to_s3(file_path):
        filetype = os.path.splitext(file_path)[1]
        filehash = hashlib.md5(open(file_path, 'rb').read()).hexdigest()
        s3_filename = f"report/{filehash}{filetype}"
        
        s3_client = boto3.client('s3')
        s3_client.upload_file(file_path, "toxindex", s3_filename)
        
        return f"s3://toxindex/{s3_filename}"

    @staticmethod
    def create_report(project_id, file_path, user_id, title, description):
        
        s3_reference = Report.upload_to_s3(file_path)

        params = (s3_reference, user_id, title, description)
        ds.execute("INSERT INTO reports (s3_reference, user_id, title, description) values (%s, %s, %s, %s)", params)
        
        # Fetch and return the newly created report
        res = ds.find("SELECT * from reports where s3_reference = %s AND user_id = %s", (s3_reference, user_id))
        report = Report.from_row(res) if res else None
        ds.execute("INSERT INTO project_reports (project_id, report_id) values (%s, %s)", (project_id, report.report_id))
        
        return report 

    @staticmethod
    def associate_with_project(report_id, project_id):
        ds.execute("INSERT INTO project_reports (project_id, report_id) values (%s, %s)", (project_id, report_id))

    @staticmethod
    def get_reports_by_project(project_id):
        # TODO need to actually get by project_id, but routing doesn't work for that right now.
        query = "SELECT * FROM reports r INNER JOIN project_reports pr ON pr.report_id = r.report_id"
        # query = f"{query} WHERE pr.project_id = %s"
        
        rows = ds.find_all(query, (project_id,))
        return [Report.from_row(row) for row in rows]
    
    # Other methods, like update, delete, etc., can also be implemented as needed.

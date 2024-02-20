import boto3
import hashlib
import os
import report.datastore as ds
import logging
import types
import report.s3store as s3store

class Report:
    
    STATUS = types.SimpleNamespace(PENDING='PENDING', GENERATED='GENERATED', FAILED='FAILED')

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
    def create_report(project_id, user_id, title, description):
        con = ds.get_connection()
        try:
            # Start a transaction
            con.autocommit = False
            cur = con.cursor()
            
            # Insert report
            params = (user_id, title, description)
            query = """
                INSERT INTO reports (user_id, title, description, status, s3_reference)
                VALUES (%s, %s, %s, 'PENDING', 'none')
                RETURNING report_id
            """
            cur.execute(query, params)
            report_id = cur.fetchone()[0]
            
            # Associate report with project
            query = "INSERT INTO project_reports (project_id, report_id) VALUES (%s, %s)"
            cur.execute(query, (project_id, report_id))
            
            # Commit the transaction if both operations are successful
            con.commit()
            
            # Assuming you want to return the report_id of the created report
            return report_id
        except Exception as e:
            # Roll back the transaction in case of an error
            con.rollback()
            logging.error("Failed to create report and associate it with project", exc_info=e)
            raise  # Re-raise the exception for the caller to handle
        finally:
            # Clean up
            if con:
                con.autocommit = True  # Reset autocommit to its default state
                con.close()

    @staticmethod
    def get_reports_by_project(project_id):
        query = "SELECT * FROM reports r INNER JOIN project_reports pr ON pr.report_id = r.report_id AND pr.project_id = %s"
        
        rows = ds.find_all(query, (project_id,))
        return [Report.from_row(row) for row in rows]
    
    @staticmethod
    def update_status(report_id, status):
        
        if status not in Report.STATUS.__dict__.values():
            raise ValueError(f"Invalid status: {status}")
        
        query = "UPDATE reports SET status = %s WHERE report_id = %s"
        ds.execute(query, (status, report_id))
            
    @staticmethod
    def update_s3_reference(report_id, s3_reference):
        query = "UPDATE reports SET s3_reference = %s WHERE report_id = %s"
        ds.execute(query, (s3_reference, report_id))
    
    @staticmethod
    def create_report_object(report_id, s3_reference, mimetype):
        query = "INSERT INTO report_object (report_id, s3_reference, mimetype) VALUES (%s, %s, %s)"
        ds.execute(query, (report_id, s3_reference, mimetype))
    
    @staticmethod
    def get_report_object(report_id, mimetype):
        query = "SELECT s3_reference FROM report_object WHERE report_id = %s AND mimetype = %s"
        res = ds.find(query, (report_id, mimetype))
        return res['s3_reference'] if res else None
        
    @staticmethod                            
    def delete(report_id):
        
        report = Report.get(report_id)
        con = ds.get_connection()
        s3 = s3store.S3Store()
        
        try:
            # Start a transaction
            con.autocommit = False
            cur = con.cursor()
            
            query = "DELETE FROM project_reports WHERE report_id = %s"
            cur.execute(query, (report_id,))
            
            query = "DELETE FROM reports WHERE report_id = %s"
            cur.execute(query, (report_id,))
            
            s3.delete_object(report.s3_reference)
            
            con.commit()
        except Exception as e:
            con.rollback()
            logging.error("Failed to delete report", exc_info=e)
            raise
        finally:
            con.autocommit = True
            con.close()

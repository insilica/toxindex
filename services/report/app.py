import boto3, requests
import os, re, flask, logging, pdfkit
from flask import Flask, request, jsonify, render_template_string
from werkzeug.utils import secure_filename
from report.model.reports.reprotox import ReprotoxReport
from report.model.report import Report
from report import datastore as ds
from report import s3store 
from flask import url_for


app = Flask(__name__, template_folder="templates")
s3 = s3store.S3Store()
logging.basicConfig(level=logging.DEBUG) # Use INFO for less verbose output

@app.route('/')
def index():
    project_id = request.args.get('project_id')
    reports = [r.to_dict() for r in Report.get_reports_by_project(project_id)]
    logging.info(f"reports: {reports}")
    logging.info(f"project_id: {project_id}")
    return flask.render_template("reports.html", reports=reports, project_id=project_id)

@app.route('/new_report_form')
def new_report_form():
    project_id = request.args.get('project_id')
    return flask.render_template("new_report.html", project_id=project_id)

@app.route('/generate_report', methods=['POST'])
def generate_report():
    project_id = request.json['project_id']
    inchi = request.json['inchi']
    report_title = request.json['title']

    # Create a new 'pending' report
    logging.info(f"Creating report for {inchi} in project {project_id}")
    report_id = Report.create_report(project_id, 1, report_title, "This is a test report")
    
    # Generate the report using the ReprotoxReport class
    ReprotoxReport.generate(report_id, report_title, inchi)
    logging.info(f"Generated report for {report_title}")

    return "Report generated successfully!", 200

@app.route('/delete_report', methods=['POST'])
def delete_report():
    _ = request.json['project_id']
    report_id = request.json['report_id']
    Report.delete(report_id)
    
    # Create a new 'pending' report
    logging.info(f"Deleting report {report_id}")
    

    return "Report generated successfully!", 200

@app.route('/download_report', methods=['POST'])
def download_report():
    report_id = request.json['report_id']
    
    s3_url = s3.url_to_object(Report.get(report_id).s3_reference)
    presigned_url = s3.generate_presigned_url(s3_url)
    
    return {"url": presigned_url}, 200

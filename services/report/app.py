import boto3, requests
import os, re, flask, logging, pdfkit
from flask import Flask, request, jsonify, render_template_string
from werkzeug.utils import secure_filename
from report.model.reports.reprotox import ReprotoxReport
from report.model.report import Report
from report import datastore as ds
from report import s3store 
from flask import url_for
from threading import Thread
from concurrent.futures import ThreadPoolExecutor

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
    
    if ReprotoxReport.validate_form(inchi)[0] == False:
        return jsonify({"error": "Invalid InChI provided, please copy a valid inchi"}), 400
    
    # Create a new 'pending' report
    logging.info(f"Creating report for {inchi} in project {project_id}")
    report_id = Report.create_report(project_id, 1, report_title, "This is a test report")
    
    def generate_report_async():
        ReprotoxReport.generate(report_id, report_title, inchi)

    # Generate the report using the ReprotoxReport class
    thread = Thread(target=generate_report_async)
    thread.start()
    logging.info(f"Report generation started for {report_title}")
    
    return jsonify({"message": "Report generation initiated."}), 200

@app.route('/tmp')
def tmp():
    inchi = 'InChI=1S/C9H8O4/c1-6(10)13-8-5-3-2-4-7(8)9(11)12/h2-5H,1H3,(H,11,12)'
    predictions = ReprotoxReport._generate_predictions(inchi)
    report = ReprotoxReport('a test', inchi, predictions)
    group_preds = report.prediction_df.groupby('category').apply(lambda x: x.to_dict(orient='records')).to_dict()
    report_data = {
        "title": report.title,
        "generated_on": report.generated_on.strftime('%Y-%m-%d %H:%M'),
        "inchi": report.inchi,
        "smiles": report.smiles
    }
    return flask.render_template("reports/reprotox.html", report=report_data, grouped_data = group_preds)

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

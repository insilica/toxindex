import boto3, requests
import os, re, flask, logging, pdfkit
from flask import Flask, request, jsonify, render_template_string
from werkzeug.utils import secure_filename
from report.model.reports.reprotox import ReprotoxReport
from report.model.report import Report
from report import datastore as ds

app = Flask(__name__, template_folder="templates")
logging.basicConfig(level=logging.DEBUG) # Use INFO for less verbose output

@app.route('/p/<project_id>/report/')
def reports(project_id):
    reports = [r.to_dict() for r in Report.get_reports_by_project(project_id)]
    logging.info(f"reports: {reports}")
    return flask.render_template("reports.html", reports=reports, project_id=project_id)

@app.route('/p/<project_id>/report/new_report_form')
def new_report_form(project_id):
    return flask.render_template("new_report.html", project_id=project_id)

@app.route('/p/<project_id>/report/generate_report', methods=['POST'])
def generate_report(project_id):
    inchi = request.json['inchi']
    report_title = request.json['title']

    # Generate the report using the ReprotoxReport class
    report = ReprotoxReport.from_inchi(inchi)

    html_content = report.generate_html_content("reports/reprotox.html")
    pdf_content = report.generate_pdf_content(html_content)

    html_filename = f"report.html"
    pdf_filename = f"report.pdf"
    
    report.save_report(html_filename, pdf_filename, html_content, pdf_content)
    report = Report.create_report(project_id, pdf_filename, 1, report_title, "This is a test report")
    
    logging.info(f"Generated report for {inchi} and saved HTML to {html_filename} and PDF to {pdf_filename}")

    return "Report generated successfully!", 200

@app.route('/p/<project_id>/report/download_report/<path:path>')
def download_report(project_id,path):
    return flask.send_from_directory(os.getcwd(), path, as_attachment=True, download_name="report.pdf")
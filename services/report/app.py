import boto3, requests
import os, re, flask, logging, pdfkit
from flask import Flask, request, jsonify, render_template_string
from werkzeug.utils import secure_filename
from .model.report.reprotox import ReprotoxReport

app = Flask(__name__, template_folder="templates")
logging.basicConfig(level=logging.DEBUG) # Use INFO for less verbose output

reports_data = [
        {'name': 'Report 1', 'time': '12/12/2022', 'path': 'report.pdf'},
        {'name': 'Report 2', 'time': '10/10/2022', 'path': 'report.pdf'}
    ]

@app.route('/')
def reports():
    # You can fetch your reports data from your database
    # Here's a sample placeholder data
    
    return flask.render_template("reports.html", reports=reports_data)


@app.route('/generate_report', methods=['POST'])
def generate_report():
    inchi = request.form.get('inchi')
    report_name = request.form.get('report_name')

    # Generate the report using the ReprotoxReport class
    report = ReprotoxReport.from_inchi(inchi)

    html_content = report.generate_html_content("reports/reprotox.html")
    pdf_content = report.generate_pdf_content(html_content)

    html_filename = f"report.html"
    pdf_filename = f"report.pdf"
    
    report.save_report(html_filename, pdf_filename, html_content, pdf_content)

    # Store the details in the reports_data object
    reports_data.append({
        'name': report_name,
        'time': '12/12/2022',
        'path_html': html_filename,
        'path_pdf': pdf_filename
    })

    logging.info(f"Generated report for {inchi} and saved HTML to {html_filename} and PDF to {pdf_filename}")

    return "Report generated successfully!", 200

@app.route('/download_report/<path:path>')
def download_report(path):
    return flask.send_from_directory(os.getcwd(), path, 
                                     as_attachment=True, 
                                     download_name="report.pdf")
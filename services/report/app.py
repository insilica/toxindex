import boto3
import os, re, flask
from flask import Flask, request, jsonify, render_template_string
from werkzeug.utils import secure_filename
from reportlab.pdfgen import canvas

app = Flask(__name__, template_folder="templates")

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

    # Placeholder code to generate a PDF
    pdf_filename = "report.pdf"
    pdf_path = os.path.join(os.getcwd(), pdf_filename)
    c = canvas.Canvas(pdf_path)
    c.drawString(100, 750, f"InChI: {inchi}")
    c.drawString(100, 730, f"Report: {report_name}")
    c.save()

    reports_data.append({'name': report_name, 'time': '12/12/2022', 
                         'path': pdf_path})
    return jsonify({
        "status": "success",
        "message": "Report generated successfully",
        "path": pdf_path
    })

@app.route('/download_report/<path:path>')
def download_report(path):
    return flask.send_from_directory(os.getcwd(), path, 
                                     as_attachment=True, 
                                     download_name="report.pdf")
import boto3, requests
import os, re, flask, logging, pdfkit
from flask import Flask, request, jsonify, render_template_string
from werkzeug.utils import secure_filename

app = Flask(__name__, template_folder="templates")
logging.basicConfig(level=logging.DEBUG) # Use INFO for less verbose output

@app.route('/p/<project_id>/search/')
def reports(project_id):
    return flask.render_template("index.html")

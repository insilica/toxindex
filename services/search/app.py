import boto3, requests
import os, re, flask, logging, pdfkit
from flask import Flask, request, jsonify, render_template_string
from werkzeug.utils import secure_filename

app = Flask(__name__, template_folder="templates")
logging.basicConfig(level=logging.DEBUG) # Use INFO for less verbose output

@app.route('/')
def index():
    return flask.render_template("index.html")

@app.route('/search', methods=['GET'])
def search():
    # Assuming 'smiles' is sent as a query parameter in the URL
    smiles = request.args.get('smiles')
    
    # Now you can use 'smiles' variable to perform your search logic

    return flask.render_template("search.html", smiles=smiles)
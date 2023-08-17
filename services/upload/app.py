import boto3
import os
import re
import psycopg
from flask import Flask, request, jsonify
from werkzeug.utils import secure_filename

app = Flask(__name__)

# Define the InChI pattern for parsing
INCHI_PATTERN = r'InChI=\d+[a-zA-Z]+.*'

def parse_inchi_from_file(file_path):
    inchis = []
    with open(file_path, 'r') as file:
        for line in file.readlines():
            line = line.strip()
            match = re.match(INCHI_PATTERN, line)
            if match:
                inchis.append(line)
    return inchis

def upload_to_s3(file, filename):
    s3 = boto3.client('s3')
    s3.upload_file(file.filename, 'toxindex', filename)


@app.route('/upload', methods=['POST'])
def upload_file():
    if 'file' not in request.files:
        return jsonify({'error': 'No file part'}), 400

    file = request.files['file']
    if file.filename == '':
        return jsonify({'error': 'No selected file'}), 400

    if file:
        filename = secure_filename(file.filename)
        file_path = os.path.join("/tmp", filename)
        file.save(file_path)
        
        inchis = parse_inchi_from_file(file_path)

        # Connect to the database
        conn = psycopg.connect(
            host="your_host",
            database="your_database",
            user="your_user",
            password="your_password"
        )
        
        # Create a cursor to execute SQL commands
        cur = conn.cursor()
        
        # Insert into source table
        cur.execute("INSERT INTO source (file_reference) VALUES (%s) RETURNING id", (filename,))
        source_id = cur.fetchone()[0]
        
        # Insert the parsed InChIs into the source_chemicals table
        for inchi in inchis:
            cur.execute("INSERT INTO source_chemicals (source_id, inchi) VALUES (%s, %s)", (source_id, inchi,))
        
        # Commit the transaction
        conn.commit()
        
        # Close the connection
        cur.close()
        conn.close()

        # Upload file to S3 bucket
        upload_to_s3(file, filename)

        return jsonify({'message': 'File uploaded and processed successfully'}), 200

    return jsonify({'error': 'Unexpected error'}), 500

if __name__ == "__main__":
    app.run(debug=True)

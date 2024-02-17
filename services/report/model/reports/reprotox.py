import requests, pdfkit, flask, json, tempfile, uuid, boto3, dotenv, os, logging
from rdkit import Chem
from report.model.report import Report
from report import s3store
import requests
import sqlite3
import numpy as np

def inchi2smi(inchi):
    mol = Chem.MolFromInchi(inchi)
    smiles = Chem.MolToSmiles(mol)
    return smiles

class CVAEModel:
    
    @staticmethod
    def predict(inchi, property_token):
        url = "http://localhost:6515/predict"
        params = {'inchi': inchi, 'property_token': property_token}
        response = requests.get(url, params=params)
        return response.json()  # Assuming the response is JSON-formatted

import pandas as pd
cvaesql = sqlite3.connect('report/data/cvae.sqlite')
# inchi = 'InChI=1S/C9H8O4/c1-6(10)13-8-5-3-2-4-7(8)9(11)12/h2-5H,1H3,(H,11,12)'
def repro_report(inchi):
    
    # Build the property table for inchi =========================================================
    query = f"""
    SELECT * FROM property p 
    INNER JOIN property_category pc ON p.property_id = pc.property_id
    INNER JOIN category c ON pc.category_id = c.category_id
    """    
    property = pd.read_sql_query(query, cvaesql)
    query = f"SELECT property_token, value FROM activity where inchi = '{inchi}' AND value='positive'"
    activity = pd.read_sql_query(query, cvaesql)
    property = property.merge(activity, on='property_token', how='left')
    property['value'] = property['value'].fillna('unknown')
    
    # Predict properties for reprotox =========================================================
    repro = property[(property['category'] == "reproductive toxicity") & (property['strength'] > 7)]['property_token']
    devtox = property[(property['category'] == "developmental toxicity") & (property['strength'] > 7)]['property_token']
    endo = property[(property['category'] == "endocrine disruption") & (property['strength'] > 8)]['property_token']
    
    property_tokens = np.unique(list(repro) + list(devtox) + list(endo))
    
    ## generate predictions
    pdf = pd.DataFrame({'property_token': property_tokens, 'prediction': np.nan})
    for i, property_token in enumerate(pdf['property_token']):
        prediction = CVAEModel.predict(inchi, property_token)['positive_prediction']
        pdf.at[i, 'prediction'] = prediction
    
    # put it together
    output = pdf.merge(property, on='property_token')
    return output[['title','prediction','category','prediction','value']]

    
    
class ReprotoxReport:

    def __init__(self, title, inchi):
        self.title = title
        self.inchi = inchi
        self.smiles = inchi2smi(inchi)   
        
    @staticmethod
    def from_json(json):
        return ReprotoxReport(json['title'], json['inchi'])
    
    @classmethod
    def generate(cls, report_id, title, inchi):
        """Generates a pdf report, an html report and a json report"""
        
        s3 = s3store.S3Store()
        
        report = ReprotoxReport(title, inchi)
        html_content = flask.render_template("reports/reprotox.html", report=report)

        config = pdfkit.configuration(wkhtmltopdf="/usr/bin/wkhtmltopdf")
        
        with tempfile.NamedTemporaryFile(suffix='.html') as html_file, \
             tempfile.NamedTemporaryFile(suffix='.pdf') as pdf_file:

            html_file.write(html_content.encode())
            html_file.flush() 

            pdfkit.from_file(html_file.name, pdf_file.name, configuration=config)

            object_url = s3.upload_file(pdf_file.name, f'{uuid.uuid4()}.pdf')
            Report.update_s3_reference(report_id, object_url)
            Report.update_status(report_id, Report.STATUS.GENERATED)
        
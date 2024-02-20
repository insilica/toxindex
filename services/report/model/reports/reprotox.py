import requests, pdfkit, flask, json, tempfile, uuid, boto3, dotenv, os, logging
from rdkit import Chem
from report.model.report import Report
from report import s3store
import requests
import sqlite3
import numpy as np
import pandas as pd
from jinja2 import Environment, FileSystemLoader, select_autoescape

def inchi2smi(inchi):
    mol = Chem.MolFromInchi(inchi)
    smiles = Chem.MolToSmiles(mol)
    return smiles

class CVAEModel:
    
    @staticmethod
    def predict(inchi, property_token):
        url = "https://api.insilica.co/service/run/chemsim/predict"
        params = {'inchi': inchi, 'property_token': property_token}
        response = requests.get(url, params=params)
        return response.json()


class ReprotoxReport:

    jinja = Environment(
        loader=FileSystemLoader('report/templates/reports'),
        autoescape=select_autoescape(['html', 'xml'])
    )

    def __init__(self, title, inchi, prediction_df):
        self.title = title
        self.inchi = inchi
        self.smiles = inchi2smi(inchi)
        self.prediction_df = prediction_df
        self.generated_on = pd.Timestamp.now(tz='UTC')
        
    @staticmethod
    def validate_form(inchi):
        try:
            mol = Chem.MolFromInchi(inchi)
            if mol is None:
                return (False, f"Invalid InChI: {inchi}")
            return (True, "Valid InChI")
        except Exception as e:
            return (False, f"Error processing InChI: {str(e)}")
    
    # inchi = 'InChI=1S/C9H8O4/c1-6(10)13-8-5-3-2-4-7(8)9(11)12/h2-5H,1H3,(H,11,12)'
    @staticmethod
    def _generate_predictions(inchi):
        cvaesql = sqlite3.connect('report/data/cvae.sqlite')
        # Build the property table for inchi =========================================================
        query = f"""
        SELECT * FROM property p 
        INNER JOIN property_category pc ON p.property_id = pc.property_id
        INNER JOIN category c ON pc.category_id = c.category_id
        WHERE strength > 8
        """    
        property = pd.read_sql_query(query, cvaesql)
        query = f"SELECT property_token, value FROM activity where inchi = '{inchi}' AND value='positive'"
        activity = pd.read_sql_query(query, cvaesql)
        property = property.merge(activity, on='property_token', how='left')
        property['value'] = property['value'].fillna('unknown')
        
        # Predict properties for reprotox =========================================================
        repro = property[(property['category'] == "reproductive toxicity")]['property_token']
        devtox = property[(property['category'] == "developmental toxicity")]['property_token']
        endo = property[(property['category'] == "endocrine disruption") & (property['strength'] > 9)]['property_token']
        
        property_tokens = np.unique(list(repro) + list(devtox) + list(endo))
        
        ## generate predictions
        pdf = pd.DataFrame({'property_token': property_tokens, 'prediction': np.nan})
        for i, property_token in enumerate(pdf['property_token']):
            logging.info(f"Predicting property {property_token}")
            prediction = CVAEModel.predict(inchi, property_token)['positive_prediction']
            pdf.at[i, 'prediction'] = prediction
        
        # put it together
        output = pdf.merge(property, on='property_token')
        cvaesql.close()
        return output[['property_token','title','prediction','category','prediction','value']]
    
    @classmethod
    def generate(cls, report_id, title, inchi):
        """Generates a pdf report, an html report and a json report"""
        
        s3 = s3store.S3Store()
        
        predictions = cls._generate_predictions(inchi)
        report = ReprotoxReport(title, inchi, predictions)
        group_preds = report.prediction_df.groupby('category').apply(lambda x: x.to_dict(orient='records')).to_dict()
        report_data = {
            "title": report.title,
            "generated_on": report.generated_on.strftime('%Y-%m-%d %H:%M'),
            "inchi": report.inchi,
            "smiles": report.smiles
        }
        template = cls.jinja.get_template('reprotox.html')
        html_content = template.render(report=report_data, grouped_data = group_preds)

        config = pdfkit.configuration(wkhtmltopdf="/usr/bin/wkhtmltopdf")
        
        with tempfile.NamedTemporaryFile(suffix='.html') as html_file, \
             tempfile.NamedTemporaryFile(suffix='.pdf') as pdf_file:

            html_file.write(html_content.encode())
            html_file.flush() 

            pdfkit.from_file(html_file.name, pdf_file.name, configuration=config)

            object_url = s3.upload_file(pdf_file.name, f'{uuid.uuid4()}.pdf')
            Report.update_s3_reference(report_id, object_url)
            
            html_object_url = s3.upload_file(html_file.name, f'{uuid.uuid4()}.html')
            Report.create_report_object(report_id, html_object_url, 'text/html')
            
            Report.update_status(report_id, Report.STATUS.GENERATED)
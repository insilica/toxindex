import requests, pdfkit, flask, json, tempfile, uuid, boto3, dotenv, os, logging
from rdkit import Chem
from report.model.report import Report
from report import s3store

def inchi2smi(inchi):
    mol = Chem.MolFromInchi(inchi)
    smiles = Chem.MolToSmiles(mol)
    return smiles

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
        
import requests, pdfkit, flask, json
from rdkit import Chem

def inchi2smi(inchi):
    mol = Chem.MolFromInchi(inchi)
    smiles = Chem.MolToSmiles(mol)
    return smiles

class ReprotoxReport:

    def __init__(self, inchi, prediction, analogue_info):
        self.inchi = inchi
        self.smiles = inchi2smi(inchi)
        self.prediction = prediction
        self.analogue_info = analogue_info

    @classmethod
    def from_inchi(cls, inchi):

        # TODO need a model service
        # TODO need to update cvae with reprotox labels
        url = "https://api.insilica.co/service/run/chemsim/predict"
        params = { "inchi": inchi, "label": "25", "k": 5 }
        data = json.loads(requests.get(url,params).text)
        
        # Accessing specific keys in the parsed data
        weights = data['weights']
        analogue_smiles = data['analogue_smiles']
        analogue_values = data['analogue_values']
        prediction = data['prediction']
        
        analogue_info = []

        for weight, smile, value in zip(weights, analogue_smiles, analogue_values):
            analogue = { 'weight': weight, 'smile': smile, 'value': value}
            analogue_info.append(analogue)

        return cls(inchi, prediction, analogue_info)

    def generate_html_content(self, template_path):
        # Render HTML content using the given template path and this object's data
        return flask.render_template(template_path, report=self)

    def generate_pdf_content(self, html_content, config_path='/usr/bin/wkhtmltopdf'):
        # Convert HTML content to a PDF using pdfkit
        config = pdfkit.configuration(wkhtmltopdf=config_path)
        return pdfkit.from_string(html_content, False, configuration=config)

    def save_report(self, html_filename, pdf_filename, html_content, pdf_content):
        # Save HTML content
        with open(html_filename, 'w') as html_file:
            html_file.write(html_content)

        # Save PDF content
        with open(pdf_filename, 'wb') as pdf_file:
            pdf_file.write(pdf_content)


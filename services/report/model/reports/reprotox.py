import requests, pdfkit, flask, json

class ReprotoxReport:
    def __init__(self, prediction, analogue_info):
        self.prediction = prediction
        self.analogue_info = analogue_info

    @classmethod
    def from_inchi(cls, inchi):
        # API call logic
        # api_url = f"http://localhost:6500/predict?inchi={inchi}&label=25&k=2"
        # response = requests.get(api_url)
        # api_data = response.json()

        response = """{"weights": [1.0000001192092896, 0.6444193720817566, 0.6379547119140625, 0.6308428049087524, 0.6285170912742615], "analogue_smiles": ["CC(=O)Oc1ccccc1C(=O)O", "COC(=O)c1cccc(C(=O)OC)c1", "O=C(OCc1ccccc1)c1ccccc1", "O=C(OCc1ccccc1)c1ccccc1O", "CC(C=O)Cc1ccc(C(C)C)cc1"], "analogue_values": [0.0, 0.0, 0.0, 0.0, 0.0], "prediction": 0.0}"""
        api_data = json.loads(response)

        analogue_info = []

        for weight, smile, value in zip(api_data['weights'], api_data['analogue_smiles'], api_data['analogue_values']):
            analogue = {
                'weight': weight,
                'smile': smile,
                'value': value
            }
            analogue_info.append(analogue)

        return cls(api_data['prediction'], analogue_info)

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


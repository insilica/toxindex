from celery import Celery
import re

# Task: write a simple Celery task (take in a query from GCS, extract chemical names, check if the names are SMILES strings, return the SMILES strings)
# This is a simple example to demonstrate the structure of a Celery task.
app = Celery('tasks', broker='pyamqp://guest@localhost//')

def is_smiles(chemical_name):
    # Simple SMILES validation: contains only allowed characters
    smiles_pattern = r'^[A-Za-z0-9@+\-\[\]\(\)=#$%\\/\.]+$'
    return bool(re.match(smiles_pattern, chemical_name))



@app.task
def extract_smiles_from_query(query):
    # Example: assume query is a dict with 'chemical_name' key
    chemical_name = query.get('chemical_name', '')
    if is_smiles(chemical_name):
        return {'smiles': chemical_name}
    else:
        return {'smiles': None}
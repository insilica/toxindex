import redis
import json
import os
import uuid
import logging
import tempfile

# Libraries for chemical data handling
import pubchempy as pcp
from rdkit import Chem

from pathlib import Path

# For data manipulation
import pandas as pd

# Celery for task management
from workflows.celery_worker import celery

# Webserver models and storage
from webserver.model.message import MessageSchema
from webserver.model.task import Task
from webserver.model.file import File
from webserver.storage import GCSFileStorage

# Logging
logger = logging.getLogger(__name__)

# Utility functions
from workflows.utils import emit_status, download_gcs_file_to_temp, upload_local_file_to_gcs
## if this does not work, add the functions directly

# Import the Sygma tool
from sygma_predictor import SygmaMetabolitePredictor, predict_metabolites_dataframe, predict_single, get_metabolites

# Define the task for metabolite SYGMA analysis (main function)
@celery.task(bind=True)
def metabolite_sygma_task():
    """Task to perform metabolite SYGMA analysis.

    General Steps:
    1. Download input files from GCS. Format?? 
    Should be a chemical name or SMILES string for Sygma.
    2. If input is a chemical name, convert to SMILES. Try NIH tool or PubChemPy.
    3. Call the sygma_predictor tool with the SMILES string.
    4. Save the results to a temporary file.
    5. Upload the results file to GCS.
    6. Emit task status updates.
    7. Return the GCS path of the results file.

    Notes: 
    * test with hardcoded SMILES or chemical names first.
    * is workflow linear? ask Yifan/Zaki since they made it
        * how to interpret the query -> tool-calling
    * probra task has code snippet to get message payload
    * check cache for similar queries (if applicable)

    Inputs: 

    Outputs:

    """
    try:
        logger.info(f"[Metabolite SyGMa GCS] Running metabolite_sygma_task for task_id={payload.get('task_id')}, user_id={payload.get('user_id')}, payload={payload}")
        # set up the 
        r = redis.Redis(
            host=os.environ.get("REDIS_HOST", "localhost"),
            port=int(os.environ.get("REDIS_PORT", "6379"))
        )
        task_id = payload.get("task_id")
        user_id = payload.get("user_id")

        # Allow users to upload a file with chemical names or SMILES strings
        file_id = payload.get("payload")
        if not all([task_id, user_id, file_id]):
            raise ValueError(f"Missing required fields. task_id={task_id}, user_id={user_id}, file_id={file_id}")
    pass

# Define helper functions for the task (Add unit tests for each helper function)
def convert_chemical_name_to_smiles(chemical_name: str) -> str:
    """Convert a chemical name to SMILES using PubChemPy or NIH tool."""
    import pubchempy as pcp
    compound = pcp.get_compounds(chemical_name, 'name')
    if compound:
        return compound[0].canonical_smiles
    else:
        # Add user input validation and error handling
        raise ValueError(f"Could not find SMILES for chemical name: {chemical_name}")
    ### Reminder: Add unit tests for this function

def check_if_SMILES(chemical_string: str) -> bool:
    """Check if a given string is a valid SMILES."""
    try:
        if Chem.MolFromSmiles(chemical_string):
            return True
    except Exception as e:
        logger.error(f"String is not a valid SMILES: {chemical_string}, error: {e}")
        return False

def sygma_get_metabolites(parent_smiles: str, phase1_cycles=1, phase2_cycles=1) -> pd.DataFrame:
    """Get metabolites using the SygmaMetabolitePredictor."""
    predictor = SygmaMetabolitePredictor()
    metabolites = predictor.get_metabolites(parent_smiles, phase1_cycles, phase2_cycles)
    return pd.DataFrame(metabolites)

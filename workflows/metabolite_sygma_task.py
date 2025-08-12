# nix develop failed because of dependency conflicts with other tools -> move to separate Docker containers + Orchestrate
# Implement this tool as a stand alone Docker container

import re
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
def metabolite_sygma_task(self, payload):
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

    More notes:
    * SyGMa generates metabolites in isolation, but in practice, the context of metabolic pathways matters.
    This task should run in conjunction with the pathway analysis task to provide a more comprehensive view of the metabolites.
    Should customize the reaction rules based on the chemical and its context.

    Inputs: 
    - payload: Dictionary containing the task ID, user ID, query, and path to the input file.

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

        user_query = payload.get("payload", "Is Gentamicin nephrotoxic?")

        if not all([task_id, user_id, user_query]):
            raise ValueError(f"Missing required fields. task_id={task_id}, user_id={user_id}, user_query={user_query}")
        # Allow users to upload a file with chemical names or SMILES strings
        # file_id = payload.get("payload")

        chemicals = get_chemical_from_query(user_query)

        for chemical in chemicals:
            is_smiles = check_if_SMILES(chemical)
            if not is_smiles:
                # If not a valid SMILES, convert chemical name to SMILES
                chemical = convert_chemical_name_to_smiles(chemical)
                logger.info(f"Converted chemical name to SMILES: {chemical}")
            else:
                logger.info(f"Using provided SMILES: {chemical}")
            # Define reaction rules based on the chemical
            reaction_rules = define_reaction_rules(chemical)
            logger.info(f"Defined reaction rules for {chemical}: {reaction_rules}")
            # Get metabolites using the SygmaMetabolitePredictor
            metabolites_df = sygma_get_metabolites(chemical, reaction_rules=reaction_rules)
            logger.info(f"Retrieved metabolites for {chemical}: {metabolites_df.shape[0]} metabolites found.")
            # Save the metabolites DataFrame to a temporary file
            with tempfile.NamedTemporaryFile(delete=False, suffix=".csv") as temp_file:
                metabolites_df.to_csv(temp_file.name, index=False)
                temp_file_path = temp_file.name
                logger.info(f"Saved metabolites to temporary file: {temp_file_path}")
            # Upload the temporary file to GCS
            gcs_path = f"metabolites/{user_id}/{task_id}/{uuid.uuid4()}.csv"
            upload_local_file_to_gcs(temp_file_path, gcs_path)
            logger.info(f"Uploaded metabolites file to GCS: {gcs_path}")
            # Emit task status update
            emit_status(
                task_id=task_id,
                user_id=user_id,
                status="completed",
                message=f"Metabolites for {chemical} saved to {gcs_path}",
                gcs_path=gcs_path
            )

        logger.info(f"[Metabolite SyGMa GCS] Completed metabolite_sygma_task for task_id={task_id}, user_id={user_id}")
    except Exception as e:
        logger.error(f"[Metabolite SyGMa GCS] Error in metabolite_sygma_task: {e}")
        emit_status(
            task_id=task_id,
            user_id=user_id,
            status="failed",
            message=str(e)
        )
    
        # Handle any exceptions that occur during the task
        # This could include logging the error, sending a notification, etc.
        # For now, we just log the error and emit a failure status
        logger.error(f"Error in metabolite_sygma_task: {e}")

# Define helper functions for the task (Add unit tests for each helper function)
def parse_user_query(user_query: str) -> dict:
    """Parse the user query to determine steps for metabolite prediction.
    This function can be extended to handle more complex queries in the future."""
    # This is a placeholder for more complex parsing logic
    # For now, we assume the user query is a single chemical name or SMILES string
    return {"chemical": user_query.strip()}

def get_chemical_from_query(user_query: str) -> str:
    """
        Extract chemical names and/or SMILES strings from the user query using a simple LLM agent approach.
        For now, use regex to find SMILES-like strings and fallback to extracting chemical names.
    """
    # Regex for SMILES (very basic, can be improved)
    smiles_pattern = r'([A-Za-z0-9@+\-\[\]\(\)\\\/%=#$]+)'
    # Try to find SMILES in the query
    smiles_matches = re.findall(smiles_pattern, user_query)
    # Filter out matches that are likely not chemical names (heuristic: length > 5)
    smiles_candidates = [s for s in smiles_matches if len(s) > 5]
    if smiles_candidates:
        return smiles_candidates[0]
    # Otherwise, fallback to extracting chemical name (assume the whole query is the name)
    return {"chemical": user_query.strip()}

def define_reaction_rules(chemical: str) -> List[str]:
    """Define reaction rules based on the chemical name or SMILES string.
    Use an LLM to generate reaction rules that are appropriate for the chemical."""
    # For now, we return a dummy set of rules.
    # In a real implementation, this would involve calling an LLM or a database of chemical reactions.
    # This is a placeholder for more complex logic.
    rules = []
    if not chemical:
        raise ValueError("Chemical name or SMILES string cannot be empty.")
    try:
        logger.info(f"Defining reaction rules for chemical: {chemical}")
    # Example of a simple rule definition
    # In practice, this would be more complex and based on chemical properties
        rules = ["Oxidation", "Reduction", "Hydrolysis", "Hydroxylation", "Dealkylation", "Deamination"]

    except Exception as e:
        logger.error(f"Error defining reaction rules for chemical {chemical}: {e}")
        raise ValueError(f"Failed to define reaction rules for chemical: {chemical}")
    return rules


def convert_chemical_name_to_smiles(chemical_name: str) -> str:
    """Convert a chemical name to SMILES using PubChemPy or NCI tool.
    Kyu suggested using an LLM agent."""
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
    """Get metabolites using the SygmaMetabolitePredictor.
    
    Notes: Add in `reaction_rules: List[str]` parameter to customize the reaction rules 
    based on the chemical and its context.
    """
    predictor = SygmaMetabolitePredictor()
    metabolites = predictor.get_metabolites(parent_smiles, phase1_cycles, phase2_cycles)

    if not metabolites:
        raise ValueError(f"No metabolites found for SMILES: {parent_smiles}")
        logger.error(f"No metabolites found for SMILES: {parent_smiles}")
    with open("metabolites.json", "w") as f:
        json.dump(metabolites, f, indent=4)
    logger.info(f"Found {len(metabolites)} metabolites for SMILES: {parent_smiles}")
    
    return pd.DataFrame(metabolites)


# Just for testing purposes, remove later
if __name__ == "__main__":
    test_payload = {
        "task_id": "dev-test",
        "user_id": "local-user",
        "payload": "acetaminophen"
    }
    metabolite_sygma_task.apply(args=(test_payload,))

""" celery_template_simple.py
========================================
This module defines a basic Celery task that demonstrates the
structure of a workflow in the application. 

The example workflow reads a text query from a file stored in Google Cloud Storage (GCS), 
tokenizes the contents to identify potential chemical names, filters those
tokens for valid SMILES strings, and returns the list of SMILES.

This code herein can serve as scaffolding for more complex
workflows. It demonstrates how to interact with the shared ``celery``
application, emit status updates via Redis, and download input
data using the helper functions in the ``workflows.utils`` module.

"""

from __future__ import annotations

import logging
import re
import tempfile
from pathlib import Path
from typing import Dict, List

from celery import Celery  # maybe not necessary to import here, but for clarity
from workflows.celery_worker import celery
from workflows.utils import emit_status, download_gcs_file_to_temp, upload_local_file_to_gcs


def is_smiles(token: str) -> bool:
    """
    Performs a lightweight heuristic check to determine whether or not a given token looks like a SMILES string.

    Args:
        token (str): The string to check. This token should be a chemical name or SMILES string, provided by the user.
    
    Outputs:
        bool: True if the token appears to be a valid SMILES string, False otherwise.

    Exceptions:
        TypeError: If the input is not a string.
        ValueError: If the input string is empty.
    """
    if not isinstance(token, str):
        raise TypeError("Input must be a string.")
    if not token:
        raise ValueError("Input string cannot be empty.")
    
    # Allowed characters in typical SMILES strings.  This includes
    # uppercase and lowercase letters, digits, and common symbols such
    # as parentheses, brackets, bond indicators and stereochemistry
    # markers.  See https://en.wikipedia.org/wiki/Simplified_molecular_input_line_entry_system
    if not re.match(r'^[A-Za-z0-9@+\-\[\]\(\)=#$%/\\\.]+$', token):
        return False
    
    # Require at least one alphabetic character so that plain numbers
    # like ring closure digits are not treated as SMILES on their own.
    return any(ch.isalpha() for ch in token)

@celery.task(bind=True, name="extract_smiles_task")
def extract_smiles_task(self, gcs_path: str, payload: Dict[str, str]) -> Dict[str, List[str]]:
    """
    Celery task to extract SMILES strings from a query file, stored in Google Cloud Storage (GCS), containing chemical names or SMILES strings.
    This task is designed to be simple and illustrative of how to structure a Celery workflow.

    The task workflow is as follows:

    1. Emit a status update indicating that the task has started.
    2. Download the query text from Google Cloud Storage into a
       temporary directory using the ``download_gcs_file_to_temp`` helper.
    3. Tokenize the text on whitespace to find potential chemical names.
    4. Filter these tokens using a heuristic to identify SMILES strings.
    5. Emit a completion status and return the list of SMILES strings.

    If any exception occurs during processing, the task will emit a
    failure status and re‑raise the exception so that Celery marks the
    task as failed.

    Args:
        gcs_path (str): The Google Cloud Storage (GCS) path to the file containing chemical names or SMILES strings.
        payload (Dict[str, str]): Additional payload data, if needed.
    
    Returns:
        Dict[str, List[str]]: A dictionary with a list of valid SMILES strings extracted from the file.
    """
    logger = logging.getLogger(__name__)
    
    task_id = payload.get('task_id')
    logger.info(f"Starting extract_smiles_task for task {task_id}")

    emit_status(task_id, "fetching file from GCS")

    try:
        # Download the GCS file to a temporary location
        with tempfile.TemporaryDirectory() as tmpdir:
            temp_dir = Path(tmpdir)
            local_file_path = download_gcs_file_to_temp(gcs_path, temp_dir)
            logger.info(
                "Downloaded query file from %s to %s",
                gcs_path,
                local_file_path
            )
        # Read the contents of the query file
        # Assume UTF-8 encoding for simplicity. May need to adjust based on actual file encoding.
            with open(local_file_path, 'r', encoding='utf-8') as file:
                query_content = file.read()
        
        # Split the query content into tokens based on whitespace
        # For more sophisticated chemical entity recognition you could employ
        # natural language processing techniques.
        tokens = re.findall(r"\S+", query_content)
        logger.debug(f"Extracted {len(tokens)} tokens from the query file.")

        # Filter tokens to find valid SMILES strings
        smiles_list = [token for token in tokens if is_smiles(token)]
        logger.info(f"Found {len(smiles_list)} valid SMILES strings in the query file.")

        emit_status(task_id, "processing complete")

        return {"smiles": smiles_list}
    
    except Exception as exc:
        # Log the error and emit a failure status.  Re‑raising the
        # exception allows Celery to handle retries if configured.
        logger.exception("Error processing extract_smiles_task: %s", exc)
        emit_status(task_id, status="FAILED")
        raise
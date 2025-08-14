# nix develop failed because of dependency conflicts with other tools -> move to separate Docker containers + Orchestrate
# Implement this tool as a stand alone Docker container

import re
import redis
import json
import os
import uuid
import logging
import tempfile

import hashlib

# Libraries for chemical data handling
from rdkit import Chem

from pathlib import Path

# For data manipulation
import pandas as pd
from typing import List, Dict, Any, Optional, Tuple, Iterable, Union

# Celery for task management
from celery_worker_sygma import celery

# Webserver models and storage
from webserver.model.message import MessageSchema
from webserver.model.task import Task
from webserver.model.file import File
from webserver.storage import GCSFileStorage
from webserver.cache_manager import cache_manager
from webserver.ai_service import convert_pydantic_to_markdown

# Utility functions
from workflows.utils import download_gcs_file_to_temp, upload_local_file_to_gcs
## if this does not work, add the functions directly

# Import the Sygma tool
from sygma_predictor import SygmaMetabolitePredictor

# Set up the logger
logging.getLogger().info("metabolite_sygma_task.py module loaded")
logger = logging.getLogger(__name__)

# <----- Utility Functions ----->
def get_redis_connection():
    """Get Redis connection with consistent configuration"""
    return redis.Redis(
        host=os.environ.get("REDIS_HOST", "localhost"),
        port=int(os.environ.get("REDIS_PORT", "6379"))
    )


def publish_to_celery_updates(event_type, task_id, data):
    """Publish event to celery_updates channel for database processing"""
    r = get_redis_connection()
    event = {
        "type": event_type,
        "task_id": task_id,
        "data": data,
    }
    r.publish("celery_updates", json.dumps(event, default=str))
    logger.info(f"Published {event_type} to celery_updates for task {task_id}")


def publish_to_socketio(event_name, room, data):
    """Publish event to Socket.IO Redis channel for real-time updates"""
    r = get_redis_connection()
    socketio_event = {
        "method": "emit",
        "event": event_name,
        "room": room,
        "data": data
    }
    r.publish("socketio", json.dumps(socketio_event, default=str))
    logger.info(f"Published {event_name} to Socket.IO for room {room}")


def emit_status(task_id, status):
    """Emit task status update to both database and real-time channels"""
    logger.info(f"[emit_status] {task_id} -> {status}")
    
    # Update database directly
    Task.set_status(task_id, status)
    task = Task.get_task(task_id)
    
    # Publish to celery_updates for any additional database processing
    publish_to_celery_updates("task_status_update", task.task_id, task.to_dict())
    
    # Publish to Socket.IO for real-time updates
    publish_to_socketio("task_status_update", f"task_{task_id}", task.to_dict())


def emit_task_message(task_id, message_data):
    """Emit task message to both database and real-time channels"""
    # Publish to celery_updates for database processing
    publish_to_celery_updates("task_message", task_id, message_data)
    
    # Publish to Socket.IO for real-time updates
    task = Task.get_task(task_id)
    if task and getattr(task, 'session_id', None):
        # Emit to chat session room
        publish_to_socketio("new_message", f"chat_session_{task.session_id}", message_data)
    
    # Emit to task room
    publish_to_socketio("task_message", f"task_{task_id}", {
        "type": "task_message",
        "data": message_data,
        "task_id": str(task_id),  # Convert UUID to string for JSON serialization
    })


def emit_task_file(task_id, file_data):
    """Emit task file to both database and real-time channels"""
    # Publish to celery_updates for database processing
    publish_to_celery_updates("task_file", task_id, file_data)
    
    # Publish to Socket.IO for real-time updates
    publish_to_socketio("task_file", f"task_{task_id}", file_data)


# <----- Helper functions: Query parsing & validation ----->
def get_chemical_from_query(user_query: str) -> List[str]:
    """
    Extract candidate chemicals from a free-text query.

    Returns a list of strings:
      - If valid SMILES are found in the text, returns the unique set of SMILES.
      - Otherwise, returns [user_query.strip()] (treat as a single name).

    Notes:
      - This is intentionally conservative; we only keep tokens that RDKit can parse.
    """
    if not user_query or not isinstance(user_query, str):
        return []

    text = user_query.strip()
    # Split on whitespace/commas/semicolons; keep “chemical-ish” tokens.
    rough_tokens = re.findall(r"[A-Za-z0-9@+\-\[\]\(\)\\/%=#$\.]+", text)

    valid_smiles: List[str] = []
    seen = set()
    for tok in rough_tokens:
        if _is_smiles(tok) and tok not in seen:
            seen.add(tok)
            valid_smiles.append(tok)

    if valid_smiles:
        logger.info(f"Detected {len(valid_smiles)} SMILES token(s) in query.")
        return valid_smiles

    # Fallback: treat the whole thing as a single chemical name
    return [text]

def _is_smiles(s: str) -> bool:
    """True if s parses as a valid SMILES with RDKit."""
    if not s or not isinstance(s, str):
        return False
    try:
        return Chem.MolFromSmiles(s) is not None
    except Exception:
        return False

def check_if_SMILES(chemical_string: str) -> bool:
    """Alias kept for backward compatibility with your code."""
    return _is_smiles(chemical_string)

def convert_chemical_name_to_smiles(chemical_name: str) -> Optional[str]:
    """
    Resolve a chemical name to a SMILES string via PubChem.
    Returns None if no hit; callers should handle None.
    """
    try:
        import pubchempy as pcp
        hits = pcp.get_compounds(chemical_name, "name")
        if hits:
            smi = hits[0].canonical_smiles
            logger.info(f"Resolved name '{chemical_name}' → SMILES '{smi}'")
            return smi
        logger.warning(f"No PubChem hit for '{chemical_name}'")
    except Exception as e:
        logger.warning(f"PubChem lookup failed for '{chemical_name}': {e}")
    return None

# <----- Helper functions: Reaction rule stub (not used now, kept for future) ----->
def define_reaction_rules(chemical: str) -> List[str]:
    """
    Placeholder for custom rule selection. SyGMa’s public API doesn’t take rules directly;
    you can wire this in later if you extend the predictor.
    """
    if not chemical:
        raise ValueError("Chemical name or SMILES string cannot be empty.")
    logger.info(f"Defining reaction rules for chemical: {chemical}")
    return ["Oxidation", "Reduction", "Hydrolysis", "Hydroxylation", "Dealkylation", "Deamination"]

# <----- SyGMa integration and normalization ----->
def _normalize_sygma_output(result: Any) -> List[Dict[str, Any]]:
    """
    Normalize SyGMa outputs to a list of {'smiles': str, 'probability': float|None}.

    Accepts:
      - list[dict] that already contain 'smiles'/'probability'
      - list[tuple|list] like (smiles, probability, *rest)
    """
    out: List[Dict[str, Any]] = []
    if not isinstance(result, list):
        return out

    for item in result:
        if isinstance(item, dict):
            smi = item.get("smiles")
            prob = item.get("probability")
        elif isinstance(item, (list, tuple)) and len(item) >= 2:
            smi, prob = item[0], item[1]
        else:
            continue

        if isinstance(smi, str):
            try:
                prob_f = float(prob) if prob is not None else None
            except Exception:
                prob_f = None
            out.append({"smiles": smi, "probability": prob_f})
    return out


def sygma_get_metabolites(
    parent_smiles: str,
    phase1_cycles: int = 1,
    phase2_cycles: int = 1,
    as_dataframe: bool = False,
) -> Union[List[Dict[str, Any]], pd.DataFrame]:
    """
    Run SyGMa and return either a normalized list of dicts or a DataFrame.

    Raises ValueError if no metabolites are found.
    """
    predictor = SygmaMetabolitePredictor()
    raw = predictor.get_metabolites(parent_smiles, phase1_cycles, phase2_cycles)
    mets = _normalize_sygma_output(raw)

    if not mets:
        logger.error(f"No metabolites found for SMILES: {parent_smiles}")
        raise ValueError(f"No metabolites found for SMILES: {parent_smiles}")

    logger.info(f"Found {len(mets)} metabolites for SMILES: {parent_smiles}")

    if as_dataframe:
        return pd.DataFrame(mets, columns=["smiles", "probability"])
    return mets

# <----- Main Function ----->
# Define the task for metabolite SYGMA analysis (main function)
@celery.task(bind=True, queue="sygma")
def metabolite_sygma_task(self, payload: Dict[str, Any]):
    """
    Celery task to run SyGMa on a name or SMILES, cache results, emit progress,
    send a markdown summary, and upload JSON/CSV to GCS—mirrors probra_task style.
    """
    logger.info("=== TASK STARTED: metabolite_sygma_task ===")
    logger.info(f"Task ID: {getattr(self.request, 'id', None)}")
    logger.info(f"Payload: {payload}")

    task_id = payload.get("task_id")
    user_id = payload.get("user_id")
    query = payload.get("payload")  # either a chemical name or a SMILES

    try:
        if not all([task_id, user_id, query]):
            raise ValueError(f"Missing required fields: task_id={task_id}, user_id={user_id}, payload={query}")

        # match probra_task structure
        emit_status(task_id, "starting")
        logger.info(f"Processing SyGMa task {task_id} for user {user_id}. Query: {query}")

        # Redis
        if get_redis_connection:
            r = get_redis_connection()
        else:
            import redis as _redis
            r = _redis.Redis(
                host=os.environ.get("REDIS_HOST", "localhost"),
                port=int(os.environ.get("REDIS_PORT", "6379")),
                decode_responses=True,
            )

        emit_status(task_id, "resolving input")

        # Resolve to SMILES:
        if _is_smiles(query):
            smiles = query
            logger.info(f"Input detected as SMILES: {smiles}")
        else:
            logger.info(f"Input detected as name; resolving '{query}' to SMILES via PubChem")
            smiles = _name_to_smiles(query)
            if not smiles:
                raise ValueError(f"Could not resolve '{query}' to a SMILES string.")

        # Cache keys (like probra_task)
        cache_key_base = f"sygma_cache:{hashlib.sha256(f'{query}|{smiles}'.encode()).hexdigest()}"
        json_cache_key = f"{cache_key_base}:json"
        csv_cache_key = f"{cache_key_base}:csv"
        md_cache_key = f"{cache_key_base}:markdown"

        emit_status(task_id, "checking cache")
        cached_json = r.get(json_cache_key)
        cached_csv = r.get(csv_cache_key)
        cached_md = r.get(md_cache_key)

        if cached_json and cached_md and cached_csv:
            logger.info("Using cached SyGMa results")
            metabolites = json.loads(cached_json)
            markdown = cached_md
            csv_data = cached_csv
            emit_status(task_id, "using cache")
        else:
            # Run SyGMa
            emit_status(task_id, "running sygma")
            logger.info("Instantiating SygmaMetabolitePredictor")
            predictor = SygmaMetabolitePredictor()

            logger.info(f"Predicting metabolites for SMILES: {smiles}")
            raw = predictor.get_metabolites(smiles, phase1_cycles=1, phase2_cycles=1)
            metabolites = _normalize_metabolites(raw)

            emit_status(task_id, "formatting results")
            logger.info(f"Normalized {len(metabolites)} metabolites")

            # Make CSV
            import pandas as pd
            df = pd.DataFrame(metabolites, columns=["smiles", "probability"])
            csv_data = df.to_csv(index=False)

            # Markdown summary (like probra_task makes markdown output)
            markdown = _to_markdown_summary(query=query, smiles=smiles, rows=metabolites)

            # Cache for 24h
            emit_status(task_id, "caching results")
            r.set(json_cache_key, json.dumps(metabolites, ensure_ascii=False), ex=60 * 60 * 24)
            r.set(csv_cache_key, csv_data, ex=60 * 60 * 24)
            r.set(md_cache_key, markdown, ex=60 * 60 * 24)

        # Emit message to user (markdown)
        emit_status(task_id, "sending message")
        message = MessageSchema(role="assistant", content=markdown)
        emit_task_message(task_id, message.model_dump())

        # Create temp files and upload to GCS
        emit_status(task_id, "uploading files to GCS")
        logger.info("Preparing temp files for upload")

        json_filename = f"sygma_result_{uuid.uuid4().hex}.json"
        csv_filename = f"sygma_result_{uuid.uuid4().hex}.csv"
        md_filename = f"sygma_result_{uuid.uuid4().hex}.md"

        with tempfile.NamedTemporaryFile(mode="w", suffix=".json", delete=False, encoding="utf-8") as tj:
            json_path = tj.name
            tj.write(json.dumps({
                "query": query,
                "resolved_smiles": smiles,
                "metabolites": metabolites
            }, indent=2, ensure_ascii=False))

        with tempfile.NamedTemporaryFile(mode="w", suffix=".csv", delete=False, encoding="utf-8") as tc:
            csv_path = tc.name
            tc.write(csv_data)

        with tempfile.NamedTemporaryFile(mode="w", suffix=".md", delete=False, encoding="utf-8") as tm:
            md_path = tm.name
            tm.write(markdown)

        try:
            gcs = GCSFileStorage()
            base = f"tasks/{task_id}"

            json_gcs_path = f"{base}/{json_filename}"
            csv_gcs_path = f"{base}/{csv_filename}"
            md_gcs_path = f"{base}/{md_filename}"

            gcs.upload_file(json_path, json_gcs_path, content_type="application/json")
            gcs.upload_file(csv_path, csv_gcs_path, content_type="text/csv")
            gcs.upload_file(md_path, md_gcs_path, content_type="text/markdown")

            # optional: warm file cache for UI (like probra_task)
            if cache_manager:
                cache_manager.cache_file_content(json_gcs_path, open(json_path, "r", encoding="utf-8").read())
                cache_manager.cache_file_content(csv_gcs_path, open(csv_path, "r", encoding="utf-8").read())
                cache_manager.cache_file_content(md_gcs_path, open(md_path, "r", encoding="utf-8").read())

            emit_status(task_id, "files uploaded")

            # Emit file events
            emit_task_file(task_id, {
                "user_id": user_id,
                "filename": json_filename,
                "filepath": json_gcs_path,
                "file_type": "json",
                "content_type": "application/json",
            })
            emit_task_file(task_id, {
                "user_id": user_id,
                "filename": csv_filename,
                "filepath": csv_gcs_path,
                "file_type": "csv",
                "content_type": "text/csv",
            })
            emit_task_file(task_id, {
                "user_id": user_id,
                "filename": md_filename,
                "filepath": md_gcs_path,
                "file_type": "markdown",
                "content_type": "text/markdown",
            })

        finally:
            # cleanup temps
            for p in (locals().get("json_path"), locals().get("csv_path"), locals().get("md_path")):
                if p:
                    try:
                        os.unlink(p)
                    except Exception as e:
                        logger.warning(f"Failed to remove temp file {p}: {e}")

        # Finish
        logger.info(f"SyGMa task {task_id} completed successfully.")
        finished_at = Task.mark_finished(task_id)
        emit_status(task_id, "done")
        return {"done": True, "finished_at": finished_at}

    except Exception as e:
        logger.error(f"Error in metabolite_sygma_task: {e}", exc_info=True)
        try:
            emit_status(task_id, "error")
        except Exception:
            pass
        raise

# @celery.task(bind=True)
# def metabolite_sygma_task(self, payload):
#     """Task to perform metabolite SYGMA analysis.

#     General Steps:
#     1. Download input files from GCS. Format?? 
#     Should be a chemical name or SMILES string for Sygma.
#     2. If input is a chemical name, convert to SMILES. Try NIH tool or PubChemPy.
#     3. Call the sygma_predictor tool with the SMILES string.
#     4. Save the results to a temporary file.
#     5. Upload the results file to GCS.
#     6. Emit task status updates.
#     7. Return the GCS path of the results file.

#     Notes: 
#     * test with hardcoded SMILES or chemical names first.
#     * is workflow linear? ask Yifan/Zaki since they made it
#         * how to interpret the query -> tool-calling
#     * probra task has code snippet to get message payload
#     * check cache for similar queries (if applicable)

#     More notes:
#     * SyGMa generates metabolites in isolation, but in practice, the context of metabolic pathways matters.
#     This task should run in conjunction with the pathway analysis task to provide a more comprehensive view of the metabolites.
#     Should customize the reaction rules based on the chemical and its context.

#     Inputs: 
#     - payload: Dictionary containing the task ID, user ID, query, and path to the input file.

#     Outputs:

#     """
#     try:
#         logger.info(f"[Metabolite SyGMa GCS] Running metabolite_sygma_task for task_id={payload.get('task_id')}, user_id={payload.get('user_id')}, payload={payload}")
#         # set up the 
#         r = redis.Redis(
#             host=os.environ.get("REDIS_HOST", "localhost"),
#             port=int(os.environ.get("REDIS_PORT", "6379"))
#         )
#         task_id = payload.get("task_id")
#         user_id = payload.get("user_id")

#         user_query = payload.get("payload", "Is Gentamicin nephrotoxic?")

#         if not all([task_id, user_id, user_query]):
#             raise ValueError(f"Missing required fields. task_id={task_id}, user_id={user_id}, user_query={user_query}")
#         # Allow users to upload a file with chemical names or SMILES strings
#         # file_id = payload.get("payload")

#         chemicals = get_chemical_from_query(user_query)

#         for chemical in chemicals:
#             is_smiles = check_if_SMILES(chemical)
#             if not is_smiles:
#                 # If not a valid SMILES, convert chemical name to SMILES
#                 chemical = convert_chemical_name_to_smiles(chemical)
#                 logger.info(f"Converted chemical name to SMILES: {chemical}")
#             else:
#                 logger.info(f"Using provided SMILES: {chemical}")
#             # Define reaction rules based on the chemical
#             reaction_rules = define_reaction_rules(chemical)
#             logger.info(f"Defined reaction rules for {chemical}: {reaction_rules}")
#             # Get metabolites using the SygmaMetabolitePredictor
#             metabolites_df = sygma_get_metabolites(chemical, reaction_rules=reaction_rules)
#             logger.info(f"Retrieved metabolites for {chemical}: {metabolites_df.shape[0]} metabolites found.")
#             # Save the metabolites DataFrame to a temporary file
#             with tempfile.NamedTemporaryFile(delete=False, suffix=".csv") as temp_file:
#                 metabolites_df.to_csv(temp_file.name, index=False)
#                 temp_file_path = temp_file.name
#                 logger.info(f"Saved metabolites to temporary file: {temp_file_path}")
#             # Upload the temporary file to GCS
#             gcs_path = f"metabolites/{user_id}/{task_id}/{uuid.uuid4()}.csv"
#             upload_local_file_to_gcs(temp_file_path, gcs_path)
#             logger.info(f"Uploaded metabolites file to GCS: {gcs_path}")
#             # Emit task status update
#             emit_status(
#                 task_id=task_id,
#                 user_id=user_id,
#                 status="completed",
#                 message=f"Metabolites for {chemical} saved to {gcs_path}",
#                 gcs_path=gcs_path
#             )

#         logger.info(f"[Metabolite SyGMa GCS] Completed metabolite_sygma_task for task_id={task_id}, user_id={user_id}")
#     except Exception as e:
#         logger.error(f"[Metabolite SyGMa GCS] Error in metabolite_sygma_task: {e}")
#         emit_status(
#             task_id=task_id,
#             user_id=user_id,
#             status="failed",
#             message=str(e)
#         )
    
#         # Handle any exceptions that occur during the task
#         # This could include logging the error, sending a notification, etc.
#         # For now, we just log the error and emit a failure status
#         logger.error(f"Error in metabolite_sygma_task: {e}")

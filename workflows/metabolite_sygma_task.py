# # nix develop failed because of dependency conflicts with other tools -> move to separate Docker containers + Orchestrate
# # Implement this tool as a stand alone Docker container

import re
import json
import os
import uuid
import logging
import tempfile
import hashlib

# Libraries for chemical data handling
from rdkit import Chem

# For data manipulation
import pandas as pd
from typing import List, Dict, Any, Optional, Union

# Celery for task management
from workflows.celery_app import celery

# Webserver models and storage
from webserver.model.message import MessageSchema
from webserver.model.task import Task
from webserver.storage import GCSFileStorage
from webserver.cache_manager import cache_manager

# Utility functions
from workflows.utils import emit_status, get_redis_connection, publish_to_celery_updates, publish_to_socketio

# Import the Sygma tool (installed from submodule in Dockerfile.sygma)
from sygma_predictor import SygmaMetabolitePredictor

# --- Gating imports ---
from workflows.sygma_docker.gating.runner import load_rule_config, gate_reactions

# Set up the logger
logging.getLogger().info("metabolite_sygma_task.py module loaded")
logger = logging.getLogger(__name__)


# --- Gating config locations (Option A paths) ---
RULES_YAML = "workflows/sygma_docker/resources/sygma_gating/rules.yaml"
CONTEXTS_YAML = "workflows/sygma_docker/resources/sygma_gating/contexts.yaml"

# Default biological context if not supplied by client payload
DEFAULT_GATING_CONTEXT = {"species": "human", "tissue": "liver"}


# <----- Utility Functions ----->
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


# Aliases for missing functions
def _name_to_smiles(chemical_name: str) -> Optional[str]:
    """Alias for convert_chemical_name_to_smiles"""
    return convert_chemical_name_to_smiles(chemical_name)


def _normalize_metabolites(raw_result: Any) -> List[Dict[str, Any]]:
    """Alias for _normalize_sygma_output"""
    return _normalize_sygma_output(raw_result)


def _to_markdown_summary(query: str, smiles: str, rows: List[Dict[str, Any]]) -> str:
    """Generate a markdown summary of the metabolite analysis results"""
    markdown = f"""# Metabolite Analysis Results

## Input
- **Query**: {query}
- **SMILES**: `{smiles}`

## Results
Found {len(rows)} metabolites:

"""
    for i, metabolite in enumerate(rows[:len(rows)], 1):  # Show first 10
        prob = metabolite.get('probability', 'N/A')
        markdown += f"{i}. **SMILES**: `{metabolite.get('smiles', 'N/A')}`\n"
        if prob != 'N/A' and prob is not None:
            try:
                markdown += f"   - **Probability**: {float(prob):.3f}\n"
            except Exception:
                markdown += f"   - **Probability**: {prob}\n"
        markdown += "\n"

    return markdown

# <----- Gating helpers ----->
def _gating_context_from_payload(payload: Dict[str, Any]) -> Dict[str, str]:
    """
    Extract gating context from payload if present; fallback to defaults.
    Accepts fields like payload['context'] = {'species': 'human', 'tissue': 'liver'}.
    """
    ctx = {}
    try:
        user_ctx = (payload.get("context") or {}) if isinstance(payload.get("context"), dict) else {}
        ctx["species"] = (user_ctx.get("species") or DEFAULT_GATING_CONTEXT["species"]).strip().lower()
        ctx["tissue"] = (user_ctx.get("tissue") or DEFAULT_GATING_CONTEXT["tissue"]).strip().lower()
    except Exception:
        ctx = dict(DEFAULT_GATING_CONTEXT)
    return ctx


def _gating_to_markdown(gating_report: Any) -> str:
    """
    Make a short, UI-friendly markdown block summarizing gating results.
    """
    try:
        d = gating_report.model_dump() if hasattr(gating_report, "model_dump") else (
            gating_report.to_dict() if hasattr(gating_report, "to_dict") else gating_report
        )
        ctx = d.get("context", {})
        rows = d.get("decisions", []) or []
        applicable = [r for r in rows if r.get("applicable")]
        relevant = [r for r in rows if r.get("relevant")]
        both = [r for r in rows if r.get("applicable") and r.get("relevant")]

        lines = []
        lines.append("## Pre-analysis Gating (SMARTS + Context)")
        lines.append(f"- **Context resolved**: `{ctx.get('resolved_context_id', 'unknown')}`")
        lines.append(f"- **Rules**: total {d.get('n_rules', len(rows))}, "
                     f"applicable {len(applicable)}, relevant {len(relevant)}, "
                     f"applicable ∧ relevant {len(both)}")
        if both:
            lines.append("\n**Likely-active rule hits:**")
            for r in both[:8]:  # keep it brief in chat
                rid = r.get("rule_id"); nm = r.get("rule_name"); tg = r.get("tags") or []
                lines.append(f"- `{rid}` — {nm} (tags: {', '.join(tg) if tg else '—'})")
        return "\n".join(lines) + "\n"
    except Exception:
        return "## Pre-analysis Gating\n- (Could not render gating summary)\n"


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
                prob_f = float(prob) if (prob is not None and prob != "") else None
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
    query = payload.get("user_query")  # either a chemical name or a SMILES

    try:
        if not all([task_id, user_id, query]):
            raise ValueError(f"Missing required fields: task_id={task_id}, user_id={user_id}, payload={query}")

        # match probra_task structure
        emit_status(task_id, "starting")
        logger.info(f"Processing SyGMa task {task_id} for user {user_id}. Query: {query}")

        # Redis
        r = get_redis_connection()

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
        gating_cache_key = f"{cache_key_base}:gating_json"
        gating_md_cache_key = f"{cache_key_base}:gating_md"

        emit_status(task_id, "checking cache")
        cached_json = r.get(json_cache_key)
        cached_csv = r.get(csv_cache_key)
        cached_md = r.get(md_cache_key)
        cached_gating = r.get(gating_cache_key)
        cached_gating_md = r.get(gating_md_cache_key)

        # --- Always compute gating unless both JSON + MD are cached ---
        emit_status(task_id, "gating")
        try:
            gating_cfg = load_rule_config(RULES_YAML, CONTEXTS_YAML)
            gating_ctx = _gating_context_from_payload(payload)
            if cached_gating and cached_gating_md:
                logger.info("Using cached gating report")
                gating_report_json = json.loads(cached_gating)
                gating_md = cached_gating_md.decode("utf-8") if isinstance(cached_gating_md, (bytes, bytearray)) else cached_gating_md
            else:
                logger.info(f"Running gating on SMILES with context {gating_ctx}")
                gating_report_obj = gate_reactions(
                    parent_smiles=smiles, rules=gating_cfg, context=gating_ctx
                )
                # serialize
                if hasattr(gating_report_obj, "model_dump"):
                    gating_report_json = gating_report_obj.model_dump()
                elif hasattr(gating_report_obj, "to_dict"):
                    gating_report_json = gating_report_obj.to_dict()
                else:
                    gating_report_json = gating_report_obj  # already dict-like

                gating_md = _gating_to_markdown(gating_report_obj)
                # cache for 24h
                r.set(gating_cache_key, json.dumps(gating_report_json, ensure_ascii=False), ex=60 * 60 * 24)
                r.set(gating_md_cache_key, gating_md, ex=60 * 60 * 24)
        except Exception as ge:
            logger.warning(f"Gating step failed or unavailable: {ge}")
            gating_report_json = None
            gating_md = "## Pre-analysis Gating\n- (unavailable)\n"

        if cached_json and cached_md and cached_csv:
            logger.info("Using cached SyGMa results")
            metabolites = json.loads(cached_json)
            markdown = cached_md.decode("utf-8") if isinstance(cached_md, (bytes, bytearray)) else cached_md
            csv_data = cached_csv.decode("utf-8") if isinstance(cached_csv, (bytes, bytearray)) else cached_csv
            # augment the markdown with gating section if not present already
            markdown = f"{gating_md}\n{markdown}"
        else:
            # Run SyGMa
            emit_status(task_id, "running sygma")
            logger.info("Instantiating SygmaMetabolitePredictor")
            predictor = SygmaMetabolitePredictor()

            # If the predictor supports a rule filter, pass IDs of (applicable ∧ relevant) rules.
            allowed_rule_ids = []
            try:
                if gating_report_json and isinstance(gating_report_json, dict):
                    for d in gating_report_json.get("decisions", []):
                        if d.get("applicable") and d.get("relevant"):
                            rid = d.get("rule_id")
                            if rid:
                                allowed_rule_ids.append(rid)
            except Exception:
                allowed_rule_ids = []

            logger.info(f"Predicting metabolites for SMILES: {smiles}")
            try:
                # optimistic path: submodule may accept a filter kwarg
                raw = predictor.get_metabolites(
                    smiles, phase1_cycles=1, phase2_cycles=1, allowed_rule_ids=allowed_rule_ids
                )
            except TypeError:
                # fallback: older API with no filter support
                raw = predictor.get_metabolites(smiles, phase1_cycles=1, phase2_cycles=1)

            metabolites = _normalize_metabolites(raw)

            emit_status(task_id, "formatting results")
            logger.info(f"Normalized {len(metabolites)} metabolites")

            # Make CSV
            df = pd.DataFrame(metabolites, columns=["smiles", "probability"])
            csv_data = df.to_csv(index=False)

            # Markdown summary (prepend gating section)
            results_md = _to_markdown_summary(query=query, smiles=smiles, rows=metabolites)
            markdown = f"{gating_md}\n{results_md}"

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
        md_filename  = f"sygma_result_{uuid.uuid4().hex}.md"
        gating_filename = f"sygma_gating_{uuid.uuid4().hex}.json"   # NEW

        with tempfile.NamedTemporaryFile(mode="w", suffix=".json", delete=False, encoding="utf-8") as tj:
            json_path = tj.name
            tj.write(json.dumps({
                "query": query,
                "resolved_smiles": smiles,
                "gating": gating_report_json,                    # embed gating in main JSON
                "metabolites": metabolites
            }, indent=2, ensure_ascii=False))

        # Gating JSON (separate downloadable artifact)
        gating_path = None
        try:
            if gating_report_json is not None:
                with tempfile.NamedTemporaryFile(mode="w", suffix=".json", delete=False, encoding="utf-8") as tg:
                    gating_path = tg.name
                    tg.write(json.dumps(gating_report_json, indent=2, ensure_ascii=False))
        except Exception as e:
            logger.warning(f"Failed to write gating JSON temp file: {e}")

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
            csv_gcs_path  = f"{base}/{csv_filename}"
            md_gcs_path   = f"{base}/{md_filename}"
            gating_gcs_path = f"{base}/{gating_filename}" if gating_path else None

            gcs.upload_file(json_path, json_gcs_path, content_type="application/json")
            gcs.upload_file(csv_path,  csv_gcs_path,  content_type="text/csv")
            gcs.upload_file(md_path,   md_gcs_path,   content_type="text/markdown")

            if gating_path and gating_gcs_path:
                gcs.upload_file(gating_path, gating_gcs_path, content_type="application/json")

            # optional: warm file cache for UI (like probra_task)
            if cache_manager:
                cache_manager.cache_file_content(json_gcs_path, open(json_path, "r", encoding="utf-8").read())
                cache_manager.cache_file_content(csv_gcs_path,  open(csv_path, "r", encoding="utf-8").read())
                cache_manager.cache_file_content(md_gcs_path,   open(md_path, "r", encoding="utf-8").read())
                if gating_path and gating_gcs_path:
                    cache_manager.cache_file_content(gating_gcs_path, open(gating_path, "r", encoding="utf-8").read())

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
            if gating_path and gating_gcs_path:
                emit_task_file(task_id, {
                    "user_id": user_id,
                    "filename": gating_filename,
                    "filepath": gating_gcs_path,
                    "file_type": "json",
                    "content_type": "application/json",
                })

        finally:
            # cleanup temps
            for p in (locals().get("json_path"), locals().get("csv_path"), locals().get("md_path"), locals().get("gating_path")):
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

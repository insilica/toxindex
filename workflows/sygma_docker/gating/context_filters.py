from typing import Dict, List, Set
import yaml
from .schemas import ContextDef, Rule


def load_contexts_yaml(path: str) -> List[ContextDef]:
    with open(path, "r", encoding="utf-8") as f:
        data = yaml.safe_load(f) or {}
    items = data.get("contexts", [])
    ctxs = []
    for it in items:
        try:
            ctxs.append(ContextDef(**it))
        except Exception as e:
            raise ValueError(f"Invalid context entry in {path}: {it}\nError: {e}")
    return ctxs


def _context_key(species: str | None, tissue: str | None) -> str:
    s = (species or "").strip().lower() or "any"
    t = (tissue or "").strip().lower() or "any"
    return f"{s}:{t}"


def allowed_and_denied_tags_for_context(
    user_ctx: Dict[str, str], contexts: List[ContextDef]
) -> tuple[Set[str], Set[str], str]:
    """
    Resolve the user's {species, tissue} to the most specific ContextDef.
    Fallback priority: (species,tissue) → (species,any) → (any,tissue) → default/any:any → empty.
    Returns: (allowed_tags, deny_tags, resolved_context_id)
    """
    species = (user_ctx.get("species") or "").strip().lower() or None
    tissue = (user_ctx.get("tissue") or "").strip().lower() or None

    # Build lookup map
    by_key: Dict[str, ContextDef] = {}
    default_id = None
    for c in contexts:
        key = _context_key(c.species, c.tissue)
        by_key[key] = c
        if c.id.lower() in {"default", "any:any"}:
            default_id = c.id

    candidates = [
        _context_key(species, tissue),
        _context_key(species, None),
        _context_key(None, tissue),
        "any:any",
    ]
    for k in candidates:
        if k in by_key:
            ctx = by_key[k]
            return set(ctx.allowed_tags or []), set(ctx.deny_tags or []), ctx.id

    # If nothing matched
    return set(), set(), default_id or "unresolved"


def relevant_for_context(rule: Rule, allowed: set[str], denied: set[str]) -> tuple[bool, List[str]]:
    """
    Decide relevance:
      - If rule has no tags → assume relevant (no context restriction)
      - If any rule tag is denied → not relevant
      - If any rule tag is allowed → relevant
      - Else → not relevant (no overlap with allowed)
    """
    reasons: List[str] = []
    rtags = set(t.lower() for t in (rule.tags or []))
    if not rtags:
        reasons.append("No tags on rule → default relevant.")
        return True, reasons

    d = set(t.lower() for t in denied)
    if rtags & d:
        reasons.append(f"Denied by context: overlap with deny_tags {sorted(rtags & d)}.")
        return False, reasons

    a = set(t.lower() for t in allowed)
    if rtags & a:
        reasons.append(f"Context-allowed tags overlap {sorted(rtags & a)}.")
        return True, reasons

    reasons.append("No overlap with context allowed tags.")
    return False, reasons

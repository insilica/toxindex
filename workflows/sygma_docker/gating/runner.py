from typing import Dict, Optional
from rdkit import Chem
from .schemas import Rule, ContextDef, GatingConfig, GatingReport, RuleDecision
from .smarts_rules import load_rules_yaml, match_rule
from .context_filters import load_contexts_yaml, allowed_and_denied_tags_for_context, relevant_for_context


def load_rule_config(rules_yaml_path: str, contexts_yaml_path: str) -> GatingConfig:
    rules = load_rules_yaml(rules_yaml_path)
    contexts = load_contexts_yaml(contexts_yaml_path)
    return GatingConfig(rules=rules, contexts=contexts)


def gate_reactions(
    parent_smiles: str,
    rules: GatingConfig | None = None,
    *,
    rules_yaml_path: Optional[str] = None,
    contexts_yaml_path: Optional[str] = None,
    context: Optional[Dict[str, str]] = None,
) -> GatingReport:
    """
    Apply SMARTS-based applicability and simple context relevance gating.

    Supply either a loaded `rules: GatingConfig` OR the YAML paths.
    `context` may include {"species": "...", "tissue": "..."}; anything else is carried through.
    """
    if rules is None:
        if not (rules_yaml_path and contexts_yaml_path):
            raise ValueError("Provide either `rules` or both YAML paths.")
        rules = load_rule_config(rules_yaml_path, contexts_yaml_path)

    context = context or {}
    parent_mol = Chem.MolFromSmiles(parent_smiles)
    if parent_mol is None:
        raise ValueError(f"Invalid parent SMILES: {parent_smiles}")

    allowed_tags, deny_tags, resolved_ctx_id = allowed_and_denied_tags_for_context(context, rules.contexts)
    # Store what we resolved for transparency
    ctx_out = dict(context)
    ctx_out["resolved_context_id"] = resolved_ctx_id
    ctx_out["allowed_tags"] = sorted(list(allowed_tags))
    ctx_out["deny_tags"] = sorted(list(deny_tags))

    decisions: list[RuleDecision] = []
    n_applicable = 0
    n_relevant = 0

    for rule in rules.rules:
        applicable, matches, reason_app = match_rule(parent_mol, rule)
        if applicable:
            n_applicable += 1

        relevant, reasons_ctx = relevant_for_context(rule, allowed_tags, deny_tags)
        if relevant:
            n_relevant += 1

        decisions.append(
            RuleDecision(
                rule_id=rule.id,
                rule_name=rule.name,
                smarts=rule.reactant_smarts,
                tags=rule.tags or [],
                applicable=applicable,
                relevant=relevant,
                matched_atoms=[list(m) for m in matches],
                reasons=[reason_app, *reasons_ctx],
            )
        )

    report = GatingReport(
        parent_smiles=parent_smiles,
        context=ctx_out,
        n_rules=len(rules.rules),
        n_applicable=n_applicable,
        n_relevant=n_relevant,
        decisions=decisions,
    )
    return report

from typing import List, Tuple
from rdkit import Chem
import yaml
from .schemas import Rule


def load_rules_yaml(path: str) -> List[Rule]:
    with open(path, "r", encoding="utf-8") as f:
        data = yaml.safe_load(f) or {}
    items = data.get("rules", [])
    rules = []
    for it in items:
        try:
            rules.append(Rule(**it))
        except Exception as e:
            raise ValueError(f"Invalid rule entry in {path}: {it}\nError: {e}")
    return rules


def compile_smarts(smarts: str) -> Chem.Mol:
    """Compile SMARTS to an RDKit Mol, raising clear errors if invalid."""
    mol = Chem.MolFromSmarts(smarts)
    if mol is None:
        raise ValueError(f"Invalid SMARTS: {smarts}")
    return mol


def match_rule(parent_mol: Chem.Mol, rule: Rule) -> Tuple[bool, List[Tuple[int, ...]], str]:
    """
    Return (applicable, matches, reason).
    - applicable: True if SMARTS found at least once
    - matches: tuple of atom index tuples
    - reason: short message
    """
    try:
        pattern = compile_smarts(rule.reactant_smarts)
    except ValueError as e:
        return (False, [], f"SMARTS compile error: {e}")

    hits = parent_mol.GetSubstructMatches(pattern, uniquify=True)
    if len(hits) > 0:
        return (True, list(hits), f"Matched {len(hits)} occurrence(s).")
    return (False, [], "No substructure match.")

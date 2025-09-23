import pytest
from pathlib import Path

# Skip gracefully if RDKit isn't available in the current env
pytest.importorskip("rdkit")

from workflows.sygma_docker.gating.runner import load_rule_config, gate_reactions


# Resolve resources relative to the current file location
CWD = Path(__file__).resolve().parent
RULES_YAML = CWD / "sygma_docker" / "resources" / "sygma_gating" / "rules.yaml"
CONTEXTS_YAML = CWD / "sygma_docker" / "resources" / "sygma_gating" / "contexts.yaml"


@pytest.fixture(scope="session")
def cfg():
    assert RULES_YAML.exists(), f"Missing rules file: {RULES_YAML}"
    assert CONTEXTS_YAML.exists(), f"Missing contexts file: {CONTEXTS_YAML}"
    return load_rule_config(str(RULES_YAML), str(CONTEXTS_YAML))


def _decisions_by_id(report):
    return {d.rule_id: d for d in report.decisions}


def test_config_loads(cfg):
    # Basic sanity on loaded config
    assert len(cfg.rules) >= 5  # from the starter rules.yaml
    rule_ids = {r.id for r in cfg.rules}
    # Check a few expected IDs from the starter file
    for rid in {"R_AROM_OH", "R_AMIDE_HYD", "R_O_DEALK", "R_GLUC", "R_SULF"}:
        assert rid in rule_ids, f"Expected rule {rid} not found"


@pytest.mark.parametrize(
    "smiles,label",
    [
        ("Oc1ccccc1", "phenol"),                 # phenolic OH on aromatic ring
        ("COc1ccccc1", "anisole"),               # methoxy on aromatic ring (O-dealkylation)
        ("CC(=O)NC1=CC=CC=C1", "acetanilide"),   # anilide (amide + ring)
    ],
)
def test_gating_basic_applicability(cfg, smiles, label):
    ctx = {"species": "human", "tissue": "liver"}
    report = gate_reactions(parent_smiles=smiles, rules=cfg, context=ctx)

    assert report.parent_smiles == smiles
    assert report.n_rules == len(cfg.rules)
    assert "resolved_context_id" in report.context
    assert "human:liver" == report.context["resolved_context_id"]

    decisions = _decisions_by_id(report)

    # Convenience accessors
    def app(rule_id): return decisions[rule_id].applicable
    def rel(rule_id): return decisions[rule_id].relevant
    def hits(rule_id): return decisions[rule_id].matched_atoms

    # Relevance expectations for human:liver (from contexts.yaml):
    #   - CYP/UGT/SULT/phaseI/phaseII → relevant
    #   - hydrolase (amidase) → not relevant (no overlap with allowed_tags)
    assert rel("R_AROM_OH")
    assert rel("R_O_DEALK")
    assert rel("R_GLUC")
    assert rel("R_SULF")
    assert not rel("R_AMIDE_HYD")

    if label == "phenol":
        # Phenol: aromatic ring + phenolic OH
        assert app("R_AROM_OH")
        assert app("R_GLUC")
        assert app("R_SULF")
        assert not app("R_O_DEALK")
        assert not app("R_AMIDE_HYD")

        assert len(hits("R_AROM_OH")) >= 1
        assert len(hits("R_GLUC")) >= 1
        assert len(hits("R_SULF")) >= 1

    elif label == "anisole":
        # Anisole: aromatic ring + O-CH3 (O-dealkylation)
        assert app("R_AROM_OH")
        assert app("R_O_DEALK")
        assert not app("R_GLUC")      # no free OH hydrogen on the oxygen
        assert not app("R_SULF")      # no phenolic OH
        assert not app("R_AMIDE_HYD") # no amide

        assert len(hits("R_AROM_OH")) >= 1
        assert len(hits("R_O_DEALK")) >= 1

    elif label == "acetanilide":
        # Acetanilide: aromatic ring + amide
        assert app("R_AROM_OH")
        assert app("R_AMIDE_HYD")
        assert not app("R_O_DEALK")
        assert not app("R_GLUC")
        assert not app("R_SULF")

        assert len(hits("R_AROM_OH")) >= 1
        assert len(hits("R_AMIDE_HYD")) >= 1

    # Count invariants
    calc_applicable = sum(1 for d in report.decisions if d.applicable)
    assert calc_applicable == report.n_applicable

    calc_relevant = sum(1 for d in report.decisions if d.relevant)
    assert calc_relevant == report.n_relevant


def test_invalid_smiles_raises(cfg):
    with pytest.raises(ValueError):
        gate_reactions(parent_smiles="not_a_smiles", rules=cfg, context={"species": "human"})


# Run with: pytest -q workflows/test_gating_minimal.py
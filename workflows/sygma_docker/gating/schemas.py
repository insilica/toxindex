"""Pydantic models for type checking"""

from typing import List, Optional, Dict, Any
from pydantic import BaseModel, Field

class _CompatBaseModel(BaseModel):
    """Pydantic v1/v2 compatibility helper."""
    def to_dict(self) -> Dict[str, Any]:
        if hasattr(self, "model_dump"):
            return self.model_dump()  # pydantic v2
        return self.dict()  # pydantic v1


class Rule(_CompatBaseModel):
    id: str
    name: str
    reactant_smarts: str = Field(..., description="SMARTS that must be present in the parent molecule")
    tags: List[str] = Field(default_factory=list)
    description: Optional[str] = None  # free text, optional


class ContextDef(_CompatBaseModel):
    id: str
    title: Optional[str] = None
    species: Optional[str] = None
    tissue: Optional[str] = None
    allowed_tags: List[str] = Field(default_factory=list)
    deny_tags: List[str] = Field(default_factory=list)


class RuleDecision(_CompatBaseModel):
    rule_id: str
    rule_name: str
    smarts: str
    tags: List[str] = Field(default_factory=list)

    applicable: bool
    relevant: bool
    matched_atoms: List[List[int]] = Field(default_factory=list)  # each match = list of atom indices
    reasons: List[str] = Field(default_factory=list)  # human-readable notes


class GatingConfig(_CompatBaseModel):
    rules: List[Rule]
    contexts: List[ContextDef]


class GatingReport(_CompatBaseModel):
    parent_smiles: str
    context: Dict[str, Any] = Field(default_factory=dict)

    n_rules: int
    n_applicable: int
    n_relevant: int

    decisions: List[RuleDecision] = Field(default_factory=list)

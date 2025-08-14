import os
import sys
import json

import redis
import uuid
import logging
import tempfile
from datetime import datetime

from rdkit import Chem
from sygma_predictor import SygmaMetabolitePredictor

# Below are helper functions for the metabolite-sygma task workflow

def is_valid_smiles(smiles: str) -> bool:
    return Chem.MolFromSmiles(smiles) is not None

def run_sygma_on_smiles(smiles: str):
    predictor = SygmaMetabolitePredictor()
    metabolites = predictor.get_metabolites(smiles, phase1_cycles=1, phase2_cycles=1)
    return metabolites
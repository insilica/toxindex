import sys
import json
# import pubchempy as pcp
from rdkit import Chem
from sygma_predictor import SygmaMetabolitePredictor

# def convert_chemical_name_to_smiles(chemical_name: str) -> str:
#     compounds = pcp.get_compounds(chemical_name, 'name')
#     if not compounds:
#         raise ValueError(f"Could not find SMILES for: {chemical_name}")
#     return compounds[0].canonical_smiles

def is_valid_smiles(smiles: str) -> bool:
    return Chem.MolFromSmiles(smiles) is not None

def run_sygma_on_smiles(smiles: str):
    predictor = SygmaMetabolitePredictor()
    metabolites = predictor.get_metabolites(smiles, phase1_cycles=1, phase2_cycles=1)
    return metabolites

def main():
    if len(sys.argv) < 2:
        print("Usage: python sygma_scratch.py <SMILES_string>")
        sys.exit(1)

    # chemical_name = sys.argv[1]
    # print(f"Converting '{chemical_name}' to SMILES...")
    # smiles = convert_chemical_name_to_smiles(chemical_name)
    # print(f"SMILES: {smiles}")
    smiles = sys.argv[1]
    print(f"Running SyGMa on SMILES: {smiles}")

    if not is_valid_smiles(smiles):
        print(f"Invalid SMILES string: {smiles}")
        sys.exit(1)

    print("Running SyGMa metabolite prediction...")
    metabolites = run_sygma_on_smiles(smiles)

    print(f"\nFound {len(metabolites)} metabolites:")
    for idx, (smi, prob) in enumerate(metabolites):
        print(f"{idx + 1}. SMILES: {smi} | Probability: {prob:.4f}")

    output = [{"smiles": smi, "probability": prob} for smi, prob in metabolites]
    with open("metabolites.json", "w") as f:
        json.dump(output, f, indent=2)
    print("\nSaved results to metabolites.json")

if __name__ == "__main__":
    main()
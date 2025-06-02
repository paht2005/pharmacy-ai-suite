# Drug Similarity Finder 
# Phat Nguyen Cong - https://github.com/paht2005


import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs

# A basic drug list with SMILES (chemical representations)
drugs = pd.DataFrame({
    "name": ["Aspirin", "Ibuprofen", "Paracetamol", "Naproxen"],
    "smiles": [
        "CC(=O)OC1=CC=CC=C1C(=O)O",
        "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O",
        "CC(=O)NC1=CC=C(O)C=C1",
        "CC(C)CC1=CC=C(C=C1)C(C)CO"
    ]
})

def smiles_to_fp(smiles_str):
    # convert SMILES to fingerprint (vector)
    m = Chem.MolFromSmiles(smiles_str)
    if not m:
        raise Exception(f"Invalid SMILES: {smiles_str}")
    return AllChem.GetMorganFingerprintAsBitVect(m, 2, 1024)

def search_similar(input_smiles, k=3):
    """
    Try to find the most similar drugs by Tanimoto score.
    Just returns a list of top matches.
    """
    try:
        input_fp = smiles_to_fp(input_smiles)
    except Exception as err:
        return [f"Error: {err}"]

    results = []
    for _, drug in drugs.iterrows():
        try:
            fp = smiles_to_fp(drug['smiles'])
            sim = DataStructs.TanimotoSimilarity(input_fp, fp)
            results.append((drug['name'], round(sim, 3)))
        except:
            continue

    results.sort(key=lambda x: x[1], reverse=True)
    return results[:k]

# For quick test:
# output = search_similar("CC(=O)OC1=CC=CC=C1C(=O)O")
# print(output)

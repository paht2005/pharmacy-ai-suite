# ===  Drug Similarity Search===
# Phat Nguyen Cong - Github: https://github.com/paht2005

from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs
import pandas as pd

drug_df = pd.DataFrame({
    "drug": ["Aspirin", "Ibuprofen", "Paracetamol", "Naproxen"],
    "smiles": [
        "CC(=O)OC1=CC=CC=C1C(=O)O",
        "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O",
        "CC(=O)NC1=CC=C(O)C=C1",
        "CC(C)CC1=CC=C(C=C1)C(C)CO"
    ]
})
def get_fp(smiles):
    mol = Chem.MolFromSmiles(smiles)
    return AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=1024)

def search_similar_drugs(query_smiles, top_k=3):
    query_fp = get_fp(query_smiles)
    results = []
    for _, row in drug_df.iterrows():
        sim = DataStructs.TanimotoSimilarity(query_fp, get_fp(row['smiles']))
        results.append((row['drug'], sim))
    results.sort(key=lambda x: x[1], reverse=True)
    return results[:top_k]



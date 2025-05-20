# === Warhouse Management and Medical Consultant ===
# Phat Nguyen Cong - Github: https://github.com/paht2005

# Dependencies
from transformers import pipeline
import pandas as pd

# Mock Knowledge Base: Symptoms -> Disease
symptom_condition_map = {
    "fever cough fatigue": "Flu or COVID-19",
    "chest pain shortness of breath": "Possible heart issue or pneumonia",
    "headache nausea light sensitivity": "Migraine",
    "abdominal pain diarrhea": "Gastroenteritis",
    "joint pain stiffness swelling": "Arthritis"
}

# Mock Drug Similarity Data (SMILES)
from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs

drugs = [
    {"name": "Aspirin", "smiles": "CC(=O)OC1=CC=CC=C1C(=O)O"},
    {"name": "Ibuprofen", "smiles": "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O"},
    {"name": "Paracetamol", "smiles": "CC(=O)NC1=CC=C(O)C=C1"},
    {"name": "Naproxen", "smiles": "CC(C)CC1=CC=C(C=C1)C(C)CO"}
]

# Convert SMILES to fingerprints
def get_fingerprint(smiles):
    mol = Chem.MolFromSmiles(smiles)
    return AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=1024)

def find_similar_drugs(query_smiles, top_k=3):
    query_fp = get_fingerprint(query_smiles)
    scores = []
    for drug in drugs:
        target_fp = get_fingerprint(drug["smiles"])
        similarity = DataStructs.TanimotoSimilarity(query_fp, target_fp)
        scores.append((drug["name"], similarity))
    scores.sort(key=lambda x: x[1], reverse=True)
    return scores[:top_k]

# Basic Rule-based diagnosis
def diagnose(symptoms):
    for key, disease in symptom_condition_map.items():
        if all(word in symptoms.lower() for word in key.split()):
            return disease
    return "❌ Unable to identify the illness. Please consult a doctor."

# LLM-based Triage Assistant
qa_pipeline = pipeline("text-generation", model="tiiuae/falcon-7b-instruct")

def analyze_with_llm(symptoms):
    prompt = f"""
    You are a medical assistant. Analyze the following symptoms and suggest the possible illness.
    Symptoms: {symptoms}
    Suggestions:
    """
    return qa_pipeline(prompt, max_new_tokens=100)[0]['generated_text']


# Example Test

symptoms = "I have a fever, cough, and fatigue."
print("🩺 Rule-based:", diagnose(symptoms))
print("🤖 LLM:", analyze_with_llm(symptoms))

# Test Drug Similarity Finder
print("🔍 Medicine similar to Paracetamol:", find_similar_drugs("CC(=O)NC1=CC=C(O)C=C1"))

# Enhanced Medical Assistant & Drug Finder 
# Phat Nguyen Cong - https://github.com/paht2005



import logging
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem, DataStructs
from transformers import pipeline



logging.basicConfig(level=logging.INFO)

# Simulated mapping of symptoms to conditions (could be expanded)
SYMPTOM_CONDITION_MAP = {
    "fever cough fatigue": "Likely Flu or COVID-19",
    "chest pain shortness of breath": "Cardiac issue or pneumonia suspected",
    "headache nausea light sensitivity": "Symptoms align with migraine",
    "abdominal pain diarrhea": "Possibly gastroenteritis",
    "joint pain stiffness swelling": "Consistent with arthritis"
}

# Local in-memory drug DB
DRUG_DATABASE = [
    {"name": "Aspirin", "smiles": "CC(=O)OC1=CC=CC=C1C(=O)O"},
    {"name": "Ibuprofen", "smiles": "CC(C)CC1=CC=C(C=C1)C(C)C(=O)O"},
    {"name": "Paracetamol", "smiles": "CC(=O)NC1=CC=C(O)C=C1"},
    {"name": "Naproxen", "smiles": "CC(C)CC1=CC=C(C=C1)C(C)CO"},
]

class DrugSimilarityEngine:
    def __init__(self, drug_db):
        self.drug_db = drug_db

    def compute_fingerprint(self, smiles):
        mol = Chem.MolFromSmiles(smiles)
        if not mol:
            raise ValueError(f"Invalid SMILES: {smiles}")
        return AllChem.GetMorganFingerprintAsBitVect(mol, 2, 1024)

    def find_similar(self, query_smiles, top_k=3):
        try:
            query_fp = self.compute_fingerprint(query_smiles)
        except Exception as e:
            logging.warning(f"Fingerprint computation failed: {e}")
            return [str(e)]

        scored = []
        for drug in self.drug_db:
            try:
                fp = self.compute_fingerprint(drug["smiles"])
                sim = DataStructs.TanimotoSimilarity(query_fp, fp)
                scored.append((drug["name"], sim))
            except Exception as e:
                logging.error(f"Comparison failed for {drug['name']}: {e}")

        return sorted(scored, key=lambda x: x[1], reverse=True)[:top_k]


def rule_based_diagnose(input_text):
    """Diagnose based on rule-matching known symptom patterns"""
    text = input_text.strip().lower()
    for phrase, diagnosis in SYMPTOM_CONDITION_MAP.items():
        keywords = phrase.split()
        if all(k in text for k in keywords):
            return diagnosis
    return "[!] Unable to determine condition. Clinical consultation advised."


# Try loading a transformer model
try:
    llm_pipeline = pipeline("text-generation", model="tiiuae/falcon-7b-instruct")
except Exception as err:
    llm_pipeline = None
    logging.error(f"LLM pipeline load failed: {err}")


def analyze_symptoms_with_llm(symptom_description):
    if not llm_pipeline:
        return "[Error] LLM model is not available."

    base_prompt = (
        "Act as a medical expert. Given the patient's symptoms, provide possible diagnoses "
        "and recommendations. Be concise yet thorough.\n\nSymptoms: " + symptom_description + "\n\nFindings:"
    )
    try:
        gen = llm_pipeline(base_prompt, max_new_tokens=120)
        return gen[0].get("generated_text", "[No output returned]").strip()
    except Exception as e:
        return f"[LLM error: {e}]"


if __name__ == "__main__":
    sample_symptoms = "Experiencing chest pain and shortness of breath."
    print("[Rule-Based Diagnosis]", rule_based_diagnose(sample_symptoms))

    llm_response = analyze_symptoms_with_llm(sample_symptoms)
    print("[LLM Suggestion]", llm_response)

    drug_finder = DrugSimilarityEngine(DRUG_DATABASE)
    result = drug_finder.find_similar("CC(=O)NC1=CC=C(O)C=C1")
    print("[Drug Similarity Search]", result)

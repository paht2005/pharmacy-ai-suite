# === Symptom-based Drug Suggestion ===
# Phat Nguyen Cong - Github: https://github.com/paht2005

# Rule-based symptom to condition mapping
medical_knowledge = {
    "fever cough fatigue": "Flu",
    "sore throat headache": "Viral Infection",
    "chest pain shortness of breath": "Pneumonia or Heart Issue",
    "nausea vomiting diarrhea": "Gastroenteritis",
    "joint pain swelling": "Arthritis"
}

# Symptom -> Suggested Drug (simplified)
condition_to_drug = {
    "Flu": "Paracetamol",
    "Viral Infection": "Ibuprofen",
    "Pneumonia or Heart Issue": "Consult doctor",
    "Gastroenteritis": "ORS + Antiemetic",
    "Arthritis": "Naproxen"
}

def suggest_condition(symptoms: str) -> str:
    for key, condition in medical_knowledge.items():
        if all(word in symptoms.lower() for word in key.split()):
            return condition
    return "Undiagnosed condition"

def suggest_drug(symptoms: str) -> str:
    condition = suggest_condition(symptoms)
    return condition_to_drug.get(condition, "No suggested medication available.")
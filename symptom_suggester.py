"""
Symptom-Based Drug Recommendation System
Phat Nguyen Cong - https://github.com/paht2005

Description: Maps reported symptoms to likely conditions and suggests appropriate treatments.
"""

from typing import List, Optional
import re

# Database of symptoms grouped by condition
CONDITION_DATABASE = {
    "Flu": ["fever", "cough", "fatigue", "body ache"],
    "Viral Infection": ["sore throat", "headache", "mild fever"],
    "Pneumonia": ["chest pain", "shortness of breath", "productive cough"],
    "Heart Issue": ["chest pain", "shortness of breath", "dizziness"],
    "Gastroenteritis": ["nausea", "vomiting", "diarrhea", "abdominal cramps"],
    "Arthritis": ["joint pain", "swelling", "stiffness"]
}

TREATMENT_GUIDE = {
    "Flu": ["Paracetamol", "Rest", "Fluids"],
    "Viral Infection": ["Ibuprofen", "Rest"],
    "Pneumonia": ["Antibiotics", "Medical supervision"],
    "Heart Issue": ["Consult a cardiologist immediately"],
    "Gastroenteritis": ["ORS", "Antiemetic", "Hydration"],
    "Arthritis": ["Naproxen", "Physical therapy"]
}

def normalize_symptoms(input_text: str) -> List[str]:
    # Extract keywords and lowercases them
    return list(set(re.findall(r'\b[a-z\s]+\b', input_text.lower())))

def match_condition(symptom_list: List[str]) -> Optional[str]:
    matched_scores = {}
    for condition, symptoms in CONDITION_DATABASE.items():
        match_count = sum(1 for s in symptoms if any(s in user_symptom for user_symptom in symptom_list))
        if match_count > 0:
            matched_scores[condition] = match_count / len(symptoms)
    if matched_scores:
        return max(matched_scores, key=matched_scores.get)
    return None

def suggest_treatment(condition: Optional[str]) -> str:
    if condition is None:
        return "Unable to determine condition. Please consult a healthcare provider."
    return ", ".join(TREATMENT_GUIDE.get(condition, ["No known treatment available."]))

def diagnose_and_suggest(symptom_input: str) -> str:
    symptoms = normalize_symptoms(symptom_input)
    condition = match_condition(symptoms)
    treatment = suggest_treatment(condition)
    return f"Likely condition: {condition or 'Unknown'}\nRecommended treatment: {treatment}"

# Example usage:
if __name__ == "__main__":
    user_input = input("Describe your symptoms: ")
    print((user_input))

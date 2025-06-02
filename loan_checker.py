# Phat Nguyen Cong - https://github.com/paht2005

import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.pipeline import Pipeline
from sklearn.compose import ColumnTransformer
from sklearn.preprocessing import OneHotEncoder, StandardScaler
from sklearn.model_selection import train_test_split
from sklearn.exceptions import NotFittedError
import logging

# Set up basic logging
logging.basicConfig(level=logging.INFO)

# Sample dataset
data = pd.DataFrame({
    "age": [25, 45, 35, 29, 62],
    "income": [50000, 120000, 70000, 45000, 90000],
    "loan_amount": [10000, 20000, 15000, 8000, 25000],
    "credit_score": [700, 750, 680, 620, 800],
    "employment": ["employed", "self-employed", "employed", "unemployed", "employed"],
    "approved": [1, 1, 1, 0, 1]
})

# Split features and target
X = data.drop("approved", axis=1)
y = data["approved"]

# Define transformers
numeric_features = ["age", "income", "loan_amount", "credit_score"]
categorical_features = ["employment"]

preprocessor = ColumnTransformer(
    transformers=[
        ("num", StandardScaler(), numeric_features),
        ("cat", OneHotEncoder(drop="first"), categorical_features)
    ]
)

# Define model pipeline
model_pipeline = Pipeline(steps=[
    ("preprocessor", preprocessor),
    ("classifier", RandomForestClassifier(random_state=42))
])

# Train the pipeline
model_pipeline.fit(X, y)

def check_loan_eligibility(applicant: dict) -> str:
    """
    Predict loan approval for a single applicant.

    Parameters:
        applicant (dict): Dictionary with keys matching training data columns.

    Returns:
        str: Approval decision.
    """
    try:
        applicant_df = pd.DataFrame([applicant])
        pred = model_pipeline.predict(applicant_df)[0]
        return "✅ Approved" if pred == 1 else "❌ Disapproved"
    except NotFittedError:
        logging.error("Model pipeline has not been fitted.")
        return "Model is not trained."
    except Exception as e:
        logging.exception("Unexpected error during prediction.")
        return f"Prediction failed: {str(e)}"

# Example usage
# applicant_info = {
#     "age": 32,
#     "income": 60000,
#     "loan_amount": 10000,
#     "credit_score": 690,
#     "employment": "employed"
# }
# print(check_loan_eligibility(applicant_info))

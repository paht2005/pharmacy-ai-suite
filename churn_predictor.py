# Customer Churn Prediction Model
# Phat Nguyen Cong - https://github.com/paht2005

import pandas as pd
import numpy as np
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import OneHotEncoder
from sklearn.compose import ColumnTransformer
from sklearn.pipeline import Pipeline
from sklearn.metrics import classification_report

# Simulated dataset
data = pd.DataFrame({
    "tenure": [1, 12, 24, 5, 36],
    "monthly_charges": [70, 50, 65, 80, 60],
    "total_charges": [70, 600, 1500, 400, 2160],
    "contract": ["month-to-month", "yearly", "yearly", "month-to-month", "yearly"],
    "tech_support": ["no", "yes", "yes", "no", "yes"],
    "churn": [1, 0, 0, 1, 0]
})

# Define feature columns
categorical_features = ['contract', 'tech_support']
numeric_features = ['tenure', 'monthly_charges', 'total_charges']
target_column = 'churn'

# Split data
X = data.drop(columns=target_column)
y = data[target_column]
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

# Pipeline for preprocessing and classification
preprocessor = ColumnTransformer(
    transformers=[
        ('cat', OneHotEncoder(handle_unknown='ignore'), categorical_features)
    ],
    remainder='passthrough'  # Leave numeric columns as is
)

model = Pipeline(steps=[
    ('preprocessor', preprocessor),
    ('classifier', RandomForestClassifier(n_estimators=100, random_state=42))
])

# Train model
model.fit(X_train, y_train)

# Evaluate
y_pred = model.predict(X_test)
print("Model Evaluation:\n", classification_report(y_test, y_pred))

# Prediction function
def predict_churn(customer_info: dict) -> str:
    """
    Predicts if a customer is likely to churn based on their data.

    Args:
        customer_info (dict): Customer features.

    Returns:
        str: Human-readable churn prediction.
    """
    try:
        input_df = pd.DataFrame([customer_info])
        prediction = model.predict(input_df)[0]
        return "Customer is likely to churn." if prediction == 1 else "Customer is likely to stay."
    except Exception as e:
        return f"Prediction error: {e}"

# === Loan Approval Predictor ===
# Phat Nguyen Cong - Github: https://github.com/paht2005

from sklearn.ensemble import RandomForestClassifier
import pandas as pd

loan_data = pd.DataFrame({
    "age": [25, 45, 35, 29, 62],
    "income": [50000, 120000, 70000, 45000, 90000],
    "loan_amount": [10000, 20000, 15000, 8000, 25000],
    "credit_score": [700, 750, 680, 620, 800],
    "employment": ["employed", "self-employed", "employed", "unemployed", "employed"],
    "approved": [1, 1, 1, 0, 1]
})

loan_data = pd.get_dummies(loan_data, columns=["employment"], drop_first=True)
X_loan = loan_data.drop("approved", axis=1)
y_loan = loan_data["approved"]
model_loan = RandomForestClassifier().fit(X_loan, y_loan)


def check_loan_eligibility(user_data: dict) -> str:
    df = pd.DataFrame([user_data])
    df = pd.get_dummies(df)
    df = df.reindex(columns=X_loan.columns, fill_value=0)
    pred = model_loan.predict(df)[0]
    return "✅ Approved" if pred == 1 else "❌ Disapproved"


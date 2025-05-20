# ===  Churn Predictor ===
# Phat Nguyen Cong 
# Github: https://github.com/paht2005


import pandas as pd
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import LabelEncoder

# Training example
churn_data = pd.DataFrame({
    "tenure": [1, 12, 24, 5, 36],
    "monthly_charges": [70, 50, 65, 80, 60],
    "total_charges": [70, 600, 1500, 400, 2160],
    "contract": ["month-to-month", "yearly", "yearly", "month-to-month", "yearly"],
    "tech_support": ["no", "yes", "yes", "no", "yes"],
    "churn": [1, 0, 0, 1, 0]
})

le_contract = LabelEncoder()
le_support = LabelEncoder()

churn_data['contract'] = le_contract.fit_transform(churn_data['contract'])
churn_data['tech_support'] = le_support.fit_transform(churn_data['tech_support'])

X = churn_data.drop("churn", axis=1)
y = churn_data["churn"]
model_churn = RandomForestClassifier().fit(X, y)

def predict_churn(customer: dict) -> str:
    df = pd.DataFrame([customer])
    df['contract'] = le_contract.transform(df['contract'])
    df['tech_support'] = le_support.transform(df['tech_support'])
    prediction = model_churn.predict(df)[0]
    return "⚠️ Leaving Risk" if prediction == 1 else "✅ Stable"


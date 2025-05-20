# === CLI App ===

# Phat Nguyen Cong - https://github.com/paht2005

from symptom_suggester import suggest_condition, suggest_drug
from drug_similarity import search_similar_drugs
from conversation import get_doctor_reply, summarize_conversation
from management_medical import diagnose, analyze_with_llm, find_similar_drugs
from review_analyzer import analyze_review
from pos_backend import create_invoice
from churn_predictor import predict_churn
from loan_checker import check_loan_eligibility
from sales_forecast import load_sales_data, forecast_sales

import os
os.environ["USE_TF"] = "0"
def menu():
    print("\n=== Pharmacy AI CLI ===")
    print("1. Symptom Checker")
    print("2. Drug Similarity")
    print("3. Conversational Doctor")
    print("4. Medical Consultant")
    print("5. Review Analyzer")
    print("6. POS Invoice")
    print("7. Churn Prediction")
    print("8. Loan Approval")
    print("9. Sales Forecast (Print to File)")
    print("0. Exit")
    return input("Choose module (0–9): ")

def main():
    while True:
        choice = menu()

        if choice == "1":
            symptoms = input("Enter symptoms: ")
            print("Condition:", suggest_condition(symptoms))
            print("Suggested Drug:", suggest_drug(symptoms))

        elif choice == "2":
            smiles = input("Enter SMILES: ")
            for drug, score in search_similar_drugs(smiles):
                print(f"{drug}: {score:.2f}")

        elif choice == "3":
            chat = ""
            while True:
                user = input("You: ")
                if user.strip().lower() == "exit": break
                chat += f"\nPatient: {user}"
                reply = get_doctor_reply(chat)
                print("Doctor:", reply)
                chat += f"\nDoctor: {reply}"
            print("--- Summary ---")
            print(summarize_conversation(chat))

        elif choice == "4":
            symptoms = input("Enter symptoms: ")
            print("Rule-based:", diagnose(symptoms))
            print("LLM-based:", analyze_with_llm(symptoms))

        elif choice == "5":
            text = input("Paste review: ")
            result = analyze_review(text)
            print(f"Sentiment: {result['sentiment']} ({result['confidence']:.2f})")
            print("Aspects:", result['aspects'])
            print("Summary:", result['summary'])

        elif choice == "6":
            print("Enter quantities:")
            cart = {
                "001": int(input("Paracetamol: ") or 0),
                "002": int(input("Ibuprofen: ") or 0),
                "003": int(input("ORS Pack: ") or 0)
            }
            lines, total = create_invoice(cart)
            for line in lines:
                print(f"{line['product']} x{line['qty']} = ${line['subtotal']:.2f}")
            print("Total:", total)

        elif choice == "7":
            print("Enter customer data:")
            data = {
                "tenure": int(input("Tenure (months): ")),
                "monthly_charges": float(input("Monthly Charges: ")),
                "total_charges": float(input("Total Charges: ")),
                "contract": input("Contract (month-to-month/yearly): "),
                "tech_support": input("Tech Support (yes/no): ")
            }
            print(predict_churn(data))

        elif choice == "8":
            print("Enter loan info:")
            info = {
                "age": int(input("Age: ")),
                "income": float(input("Income: ")),
                "loan_amount": float(input("Loan amount: ")),
                "credit_score": int(input("Credit score: ")),
                "employment": input("Employment (employed/self-employed/unemployed): ")
            }
            print(check_loan_eligibility(info))

        elif choice == "9":
            print("Generating forecast (saved to forecast.csv)...")
            df = load_sales_data()
            model, forecast = forecast_sales(df)
            forecast[['ds', 'yhat']].to_csv("forecast.csv", index=False)
            print("Done.")

        elif choice == "0":
            print("Goodbye!")
            break

        else:
            print("Invalid choice.")

if __name__ == "__main__":
    main()
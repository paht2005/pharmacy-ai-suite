# CLI Pharmacy App - Interactive Assistant
# Phat Nguyen Cong - https://github.com/paht2005

import os
import pandas as pd

from symptom_suggester import diagnose_and_suggest
from drug_similarity import search_similar
from conversation import respond_to_conversation, get_chat_summary
from management_medical import rule_based_diagnose, analyze_symptoms_with_llm
from review_analyzer import analyze_review
from pos_backend import create_invoice
from churn_predictor import predict_churn
from loan_checker import check_loan_eligibility
from sales_forecast import generate_synthetic_sales_data, build_and_train_model, extend_forecast

os.environ["USE_TF"] = "0"  # Disable TensorFlow-based operations by default


class ModuleRunner:
    def __init__(self):
        self.running = True

    def display_menu(self):
        print("\n--- Welcome to the Pharmacy AI Assistant ---")
        print("Select a feature below:")
        options = [
            "Symptom Checker",
            "Drug Similarity Search",
            "Conversational Doctor",
            "Medical Consultant (Rule & LLM)",
            "Review Sentiment Analyzer",
            "Point-of-Sale (POS) Invoice Generator",
            "Customer Churn Prediction",
            "Loan Eligibility Checker",
            "Sales Forecasting (CSV Output)",
            "Exit"
        ]
        for idx, option in enumerate(options):
            print(f"{idx}. {option}")
        return input("Enter option (0-9): ").strip()

    def run(self):
        while self.running:
            try:
                choice = self.display_menu()
                method = getattr(self, f"run_{choice}", None)
                if method:
                    method()
                else:
                    print("Invalid selection. Please choose a valid option.")
            except Exception as e:
                print(f"[Error] {str(e)}\nPlease try again.")

    def run_0(self):
        print(" Exiting application. Stay healthy!")
        self.running = False

    def run_1(self):
        symptoms = input("Please describe the symptoms you're experiencing:\n> ")
        print("\nDiagnosis & Suggestion:")
        print(diagnose_and_suggest(symptoms))

    def run_2(self):
        smiles = input("Input SMILES string for drug structure:\n> ")
        results = search_similar(smiles)
        print("\nSimilarity Scores:")
        for name, score in results:
            print(f"  - {name}: {score:.2f}")

    def run_3(self):
        print("Chat with the AI Doctor (type 'exit' to quit):")
        chat_log = ""
        while True:
            user_input = input("You: ")
            if user_input.lower() == "exit":
                break
            chat_log += f"\nPatient: {user_input}"
            reply = respond_to_conversation(chat_log)
            print("Doctor:", reply)
            chat_log += f"\nDoctor: {reply}"
        print("\n--- Chat Summary ---")
        print(get_chat_summary(chat_log))

    def run_4(self):
        symptoms = input("Enter symptoms to consult:\n> ")
        print("\n[Rule-Based Diagnosis]")
        print(rule_based_diagnose(symptoms))
        print("\n[LLM-Based Analysis]")
        print(analyze_symptoms_with_llm(symptoms))

    def run_5(self):
        review = input("Paste a product or service review:\n> ")
        result = analyze_review(review)
        print("\n[Review Insights]")
        print(f"  - Sentiment: {result['sentiment']} ({result['confidence']:.2f})")
        print(f"  - Aspects: {', '.join(result['aspects'])}")
        print(f"  - Summary: {result['summary']}")

    def run_6(self):
        print("Enter quantities for the following products:")
        try:
            cart = {
                "001": int(input("Paracetamol: ") or 0),
                "002": int(input("Ibuprofen: ") or 0),
                "003": int(input("ORS Pack: ") or 0)
            }
        except ValueError:
            print("Invalid input. Please enter numbers.")
            return
        lines, total = create_invoice(cart)
        print("\n[Invoice]")
        for item in lines:
            print(f"  {item['product']} x{item['qty']} = ${item['subtotal']:.2f}")
        print(f"Total: ${total:.2f}")

    def run_7(self):
        print("Customer info for churn prediction:")
        try:
            data = {
                "tenure": int(input("Tenure (months): ")),
                "monthly_charges": float(input("Monthly Charges: ")),
                "total_charges": float(input("Total Charges: ")),
                "contract": input("Contract Type (month-to-month/yearly): ").strip(),
                "tech_support": input("Tech Support (yes/no): ").strip().lower()
            }
            print("\nPrediction Result:")
            print(predict_churn(data))
        except ValueError:
            print("Invalid numeric input.")

    def run_8(self):
        print("Provide loan applicant details:")
        try:
            info = {
                "age": int(input("Age: ")),
                "income": float(input("Monthly Income: ")),
                "loan_amount": float(input("Loan Amount: ")),
                "credit_score": int(input("Credit Score (300-850): ")),
                "employment": input("Employment Status: ").strip()
            }
            print("\nEligibility Result:")
            print(check_loan_eligibility(info))
        except ValueError:
            print("Please check your inputs and try again.")

    def run_9(self):
        print("Generating and saving sales forecast...")
        df = generate_synthetic_sales_data()
        model = build_and_train_model(df)
        forecast = extend_forecast(model, future_days=30)
        forecast[['ds', 'yhat']].to_csv("forecast.csv", index=False)
        print("Sales forecast saved to forecast.csv.")


if __name__ == "__main__":
    app = ModuleRunner()
    app.run()

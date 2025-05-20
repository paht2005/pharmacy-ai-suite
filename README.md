# 💊 Pharmacy AI-suite Toolkit (10-in-1 CLI & API System)

A modular AI-powered system for pharmacies, clinics, and medical retail chains. It includes **10 intelligent modules** for tasks like inventory management, medical consulting, customer prediction, drug similarity search, OCR-based pricing, and more — all unified into a **CLI interface** and **REST API service**.



---

## 📌 Table of Contents

1. [✨ Project Overview](#-project-overview)  
2. [🚀 Features](#-features)  
3. [🗂️ Project Structure](#-project-structure)
4. [🧰 Tech Stack](#-tech-stack)
5. [⚙️ Installation](#-installation)  
6. [✅ Example Output](#-example-output)
7. [🚀 Enhancements in CLI & API](#-enhancements-in-cli-&-api)
8. [🧭 Future Work](#-future-work)  
9. [📄 License](#-license)
10. [🤝 Contributing](#-contributing)
11. [📬 Contact](#-contact)

---

## ✨ Project Overview

### 1. This is a multi-module AI system designed to assist pharmacy operations and decision-making. It combines:
- AI-based medical consultation & triage
- Sales forecasting using Prophet
- OCR-powered barcode & price extraction
- Drug suggestion from symptoms
- Inventory checkout/invoice system
- NLP-based review analytics
- Machine learning-based loan/churn predictions
### 2. Run from:
- `main_cli_app.py` for interactive CLI
- `main_api_app.py` for RESTful API

---

## 🚀 Features

| Module                     | Description                                         |
|----------------------------|-----------------------------------------------------|
| ``Symptom Checker``        | Map symptoms to conditions & suggest drugs          |
| ``Drug Similarity``        | Tanimoto similarity via RDKit                       |
| ``Medical Chatbot``        | Falcon7B-based conversation agent                   |
| ``Rule & LLM Diagnosis``   | Rule-based + LLM triage assistant                   |
| ``Review Analyzer``        | Extract aspects, sentiment, summary                 |
| ``POS Invoice``            | Product checkout & subtotal calculator              |
| ``Churn Predictor``        |  Classify if customer will leave                    |
| ``Loan Approver``          | Predict loan approval based on profile              |
| ``Sales Forecast``         | Predict future sales via Prophet                    |
| ``Barcode/Price OCR``      |Detect & extract price from image                    |

---
## 🗂️ Project Structure
```
├── main_cli_app.py           # CLI app for testing all modules
├── main_api_app.py           # Flask REST API
├── requirements.txt          # Python dependencies
├── churn_predictor.py
├── conversation.py           # GPT2-based doctor agent
├── drug_similarity.py        # RDKit similarity search
├── loan_checker.py
├── management_medical.py     # Rule-based + LLM symptom analysis
├── pos_backend.py
├── review_analyzer.py        # spaCy + Transformers
├── sales_forecast.py         # Prophet forecasting
├── scan_and_price_check.py   # EasyOCR based
├── symptom_suggester.py      # Symptom rule-based lookup
├── style.css                 # Optional styling for UI
└── README.md
└── LICENSE

```
---

## 🧰 Tech Stack

| Purpose                  | Library                                 |
|------------------------  |-----------------------------------------|
| **LLM/Chat/Summary**     | ``Transformers``, ``Falcon-7B-Instruct``|
| **Medical Triage**       | ``Rule-based + LLM``                    |
| **Forecasting**          | ``Prophet``                             |
| **OCR Barcode Reader**   | ``EasyOCR``, ``OpenCV``                 |
| **NLP Review Analyzer**  | ``spaCy``, ``Transformers``             |
| **ML Prediction**        | ``Scikit-learn`` (RF)                   |
| **Drug Similarity**      | ``RDKit``                               |
| **Invoice Management**   | ``pandas``                              |
| **CLI Interface**        | ``Python I/O``                          |
| **API (optional)**       | ``Flask``                               |

---

## ⚙️ Installation

```bash
git clone https://github.com/paht2005/pharmacy-ai-suite.git
cd pharmacy-ai-suite
pip install -r requirements.txt

# CLI Mode:
python main_cli_app.py

# REST API Mode:
python main_api_app.py
```
---
## ✅ Example Output

### 1. Symptom Checker
```bash
Input: fever cough fatigue
Output: Condition: Flu | Suggested Drug: Paracetamol
```
### 2. Drug Similarity
```bash
Input: SMILES of Paracetamol
Output: Ibuprofen (0.72), Naproxen (0.68)...
```
### 3. Conversational Doctor
```bash
User: I have a headache
Doctor: Have you also experienced nausea or light sensitivity?
```
### 4. Rule vs LLM Diagnosis
```bash
Rule-based: Flu
LLM-based: Possible viral infection. Rest and hydration advised.
```
### 5. Review Analyzer
```bash
Sentiment: POSITIVE (0.96)
Aspects: ['staff', 'service']
Summary: The service was fast and staff very helpful.
```
### 6. POS Invoice
```bash
Paracetamol x2 = $3.00
ORS Pack x1 = $0.80
Total = $3.80
```
### 7. Churn Prediction
```bash
Input: {'tenure': 1, 'contract': 'month-to-month'}
Result: ⚠️ Leaving Risk
```
### 8. Loan Approval
```bash
Input: age=25, income=50000, credit_score=700
Result: ✅ Approved
```
### 9. Sales Forecast
```bash
Output file: forecast.csv with predicted 'yhat' values
```
### 10. Barcode & Price Scanner
```bash
Output: ['$5.99', '$12.49'] from image with shelf tags
```

--- 
## 🚀 Enhancements in CLI & AP
### CLI Mode (``main_cli_app.py``)
- Select options 1–10 from a simple numbered interface
- Inputs handled via keyboard
- Forecast saved as CSV

### API Mode (``main_api_app.py``)
- Exposes ``/predict_churn``, ``/diagnose``, ``/chat``, ``/forecast``, etc.
- Accepts JSON payloads
- Can be plugged into a frontend later (e.g., React/Flutter)
--- 
## 🧭 Future Work
- ✅ Add user authentication layer
- ✅ Add database backend (SQLite/PostgreSQL)
- 📷 Real-time camera price scanner
- 🌐 Host web UI (Streamlit or Flask Web)
- 🧪 Add automated tests (PyTest)
---
## 📄 License
This project is licensed under the MIT License. See the [LICENSE](./LICENSE) file for details.


---
## 🤝 Contributing
I welcome contributions to improve this project!
Feel free to:
- Submit pull requests
- Report bugs
- Suggest new features

--- 
## 📬 Contact
Contact for work: **Nguyễn Công Phát** – congphatnguyen.work@gmail.com

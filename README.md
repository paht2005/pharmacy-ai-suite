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
git clone https://github.com/paht2005/Healthcare-project.git
cd Healthcare-project
pip install -r requirements.txt

# CLI Mode:
python main_cli_app.py

# REST API Mode:
python main_api_app.py
```
---
## ✅ Example Output
- 📄 John_Doe_CV.pdf — Match Score: 87.6%
- 📄 Jane_Smith_CV.pdf — Match Score: 65.4%

--- 
## 🧭 Future Work
- Export results to CSV
- Use named-entity recognition (NER) for skill extraction
- Add job role suggestions via LLM
- Multi-language support (English + Vietnamese)

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
## 🧠 Acknowledgements

### 1. What is Ollama?
- [Ollama](https://ollama.ai/) is a tool for running lightweight open-source LLMs (like LLaMA 2, Mistral, Phi, etc.) **locally on your machine** with a simple CLI.
- In this project, Ollama is used to **summarize PDF resumes** into 3–5 bullet points using models like `llama2`.

### 2. Install Ollama (Optional for LLM Summary)
If you want to enable **AI-powered resume summarization**, install Ollama:
- Download Ollama: https://ollama.com/download
- Install a model (example: `llama2`):
```bash
ollama run llama2
```
- After installation, make sure the `ollama` command is available in your system PATH.

### 3. Troubleshooting Ollama Issues
If you see this error:
```bash
FileNotFoundError: [WinError 2] The system cannot find the file specified
```
Make sure that:
- Ollama is installed correctly
- The `ollama` CLI is in your system `PATH`
- You’ve downloaded a model (like `llama2`) with `ollama run llama2`


--- 
## 📬 Contact
Contact for work: **Nguyễn Công Phát** – congphatnguyen.work@gmail.com

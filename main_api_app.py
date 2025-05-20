# === API App ===
# FastAPI REST API 
# Phat Nguyen Cong 
# Github: https://github.com/paht2005

from fastapi import FastAPI
from pydantic import BaseModel
from symptom_suggester import suggest_condition, suggest_drug
from drug_similarity import search_similar_drugs
from conversation import get_doctor_reply, summarize_conversation
from management_medical import diagnose, analyze_with_llm
from review_analyzer import analyze_review
from pos_backend import create_invoice
from churn_predictor import predict_churn
from loan_checker import check_loan_eligibility

app = FastAPI(title="Pharmacy AI API", description="All-in-one REST API")

# ==== Request Models ====

class SymptomRequest(BaseModel):
    symptoms: str

class SMILESRequest(BaseModel):
    smiles: str

class ChatRequest(BaseModel):
    chat_log: str

class ReviewRequest(BaseModel):
    text: str

class CartRequest(BaseModel):
    cart: dict[str, int]

class ChurnRequest(BaseModel):
    tenure: int
    monthly_charges: float
    total_charges: float
    contract: str
    tech_support: str

class LoanRequest(BaseModel):
    age: int
    income: float
    loan_amount: float
    credit_score: int
    employment: str

# ==== Endpoints ====

@app.post("/symptom-checker")
def symptom_checker(req: SymptomRequest):
    return {
        "condition": suggest_condition(req.symptoms),
        "drug": suggest_drug(req.symptoms)
    }

@app.post("/drug-similarity")
def drug_similarity(req: SMILESRequest):
    return search_similar_drugs(req.smiles)

@app.post("/chat-doctor")
def chat_doctor(req: ChatRequest):
    reply = get_doctor_reply(req.chat_log)
    return {"reply": reply}

@app.post("/chat-summary")
def summarize_chat(req: ChatRequest):
    return {"summary": summarize_conversation(req.chat_log)}

@app.post("/medical-consult")
def medical_consult(req: SymptomRequest):
    return {
        "rule_based": diagnose(req.symptoms),
        "llm_based": analyze_with_llm(req.symptoms)
    }

@app.post("/review-analyze")
def review_analyze(req: ReviewRequest):
    return analyze_review(req.text)

@app.post("/pos-invoice")
def pos_invoice(req: CartRequest):
    lines, total = create_invoice(req.cart)
    return {"items": lines, "total": total}

@app.post("/churn")
def churn(req: ChurnRequest):
    return {"status": predict_churn(req.dict())}

@app.post("/loan")
def loan(req: LoanRequest):
    return {"approval": check_loan_eligibility(req.dict())}

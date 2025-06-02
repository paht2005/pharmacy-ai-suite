from fastapi import FastAPI, HTTPException
from pydantic import BaseModel, Field
from typing import Dict, Optional

# Local imports for specialized modules
from symptom_suggester import diagnose_and_suggest
from drug_similarity import search_similar
from conversation import respond_to_conversation, get_chat_summary
from management_medical import rule_based_diagnose, analyze_symptoms_with_llm
from review_analyzer import analyze_review
from pos_backend import create_invoice
from churn_predictor import predict_churn
from loan_checker import check_loan_eligibility

app = FastAPI(
    title="Pharmacy AI Gateway",
    description="Comprehensive API for pharmacy analytics and services",
    version="1.1.2"
)

# Models 
class SymptomRequest(BaseModel):
    symptoms: str = Field(..., description="Comma-separated list of symptoms")

    
class SMILESRequest(BaseModel):
    smiles: str = Field(...)

class ChatRequest(BaseModel):
    chat_log: str

class ReviewRequest(BaseModel):
    text: str
    user_id: Optional[str] = None  # Optional metadata

class CartRequest(BaseModel):
    cart: Dict[str, int]

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

# Endpoints 
@app.post("/symptoms/diagnose")
def diagnose(req: SymptomRequest):
    result = diagnose_and_suggest(req.symptoms)
    if not result:
        raise HTTPException(status_code=404, detail="No diagnosis found.")
    return {"diagnosis": result}

@app.post("/drugs/similarity")
def similar_drugs(req: SMILESRequest):
    return {"matches": search_similar(req.smiles)}

@app.post("/chat/doctor")
def doctor_chat(req: ChatRequest):
    response = respond_to_conversation(req.chat_log)
    return {"response": response}

@app.post("/chat/summary")
def summarize_chat(req: ChatRequest):
    summary = get_chat_summary(req.chat_log)
    return {"summary": summary}

@app.post("/consultation")
def consult(req: SymptomRequest):
    rule_result = rule_based_diagnose(req.symptoms)
    llm_result = analyze_symptoms_with_llm(req.symptoms)
    return {
        "rule_based": rule_result,
        "llm_analysis": llm_result
    }

@app.post("/reviews/analyze")
def review_insights(req: ReviewRequest):
    return {
        "analysis": analyze_review(req.text),
        "user_tag": req.user_id or "anonymous"
    }

@app.post("/pos/invoice")
def invoice(req: CartRequest):
    items, total = create_invoice(req.cart)
    if total <= 0:
        raise HTTPException(status_code=400, detail="Empty or invalid cart")
    return {"items": items, "total": total}

@app.post("/customer/churn")
def churn(req: ChurnRequest):
    risk = predict_churn(req.dict())
    return {"risk_score": risk}

@app.post("/loan/eligibility")
def loan_eligibility(req: LoanRequest):
    decision = check_loan_eligibility(req.dict())
    return {"eligible": decision}

# === Conversational Medical Assistant ===
# Phat Nguyen Cong - Github: https://github.com/paht2005
from transformers import pipeline

# Load conversation model (Falcon or GPT-based)
doctor_bot = pipeline("text-generation", model="tiiuae/falcon-7b-instruct")

def get_doctor_reply(chat_log: str) -> str:
    prompt = f"""
    You are a virtual doctor assistant.
    Here is the conversation so far:
    {chat_log}

    Please provide a follow-up question or advice:
    """
    reply = doctor_bot(prompt, max_new_tokens=150)[0]['generated_text']
    return reply.strip()

def summarize_conversation(chat_log: str) -> str:
    summarizer = pipeline("summarization", model="sshleifer/distilbart-cnn-12-6")
    summary = summarizer(chat_log, max_length=100, min_length=30, do_sample=False)[0]['summary_text']
    return summary


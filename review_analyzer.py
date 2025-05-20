# === Review Analyzer ===
# Phat Nguyen Cong - Github: https://github.com/paht2005

from transformers import pipeline
import spacy

# Load spaCy for aspect extraction
nlp = spacy.load("en_core_web_sm")

# Load sentiment model
sentiment_model = pipeline("sentiment-analysis")

# Summarizer
summarizer = pipeline("summarization", model="sshleifer/distilbart-cnn-12-6")

def analyze_review(text: str):
    sentiment = sentiment_model(text)[0]
    doc = nlp(text)
    aspects = set(chunk.text.lower() for chunk in doc.noun_chunks if chunk.root.dep_ in ["nsubj", "dobj", "pobj"])
    summary = summarizer(text, max_length=50, min_length=20, do_sample=False)[0]['summary_text']

    return {
        "sentiment": sentiment['label'],
        "confidence": sentiment['score'],
        "aspects": list(aspects),
        "summary": summary
    }


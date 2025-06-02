# Review Analyzer 


import logging
from typing import Dict, List, Any

import spacy
from transformers import pipeline
from transformers.pipelines.base import PipelineException

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# Load external NLP resources
try:
    nlp = spacy.load("en_core_web_sm")
except OSError as e:
    logger.error("spaCy model 'en_core_web_sm' not found. Run: python -m spacy download en_core_web_sm")
    raise e

try:
    sentiment_analyzer = pipeline("sentiment-analysis")
    summarizer = pipeline("summarization", model="sshleifer/distilbart-cnn-12-6")
except PipelineException as e:
    logger.error("Failed to load Hugging Face pipelines.")
    raise e


def extract_aspects(text: str) -> List[str]:
    """Extracts important noun-based aspects from the text."""
    doc = nlp(text)
    aspects = {
        chunk.text.lower()
        for chunk in doc.noun_chunks
        if chunk.root.dep_ in {"nsubj", "dobj", "pobj"}
    }
    return sorted(aspects)


def summarize_text(text: str, min_length: int = 20, max_length: int = 50) -> str:
    """Generates a summary for the provided text."""
    try:
        summary_output = summarizer(
            text,
            min_length=min_length,
            max_length=max_length,
            do_sample=False
        )
        return summary_output[0]["summary_text"]
    except Exception as e:
        logger.warning(f"Summarization failed: {e}")
        return ""


def analyze_review(text: str) -> Dict[str, Any]:
    """Performs sentiment analysis, aspect extraction, and summarization."""
    if not text or len(text.strip()) < 10:
        raise ValueError("Input text is too short for meaningful analysis.")

    sentiment_result = sentiment_analyzer(text)[0]
    aspects = extract_aspects(text)
    summary = summarize_text(text)

    return {
        "sentiment": sentiment_result["label"],
        "confidence": sentiment_result["score"],
        "aspects": aspects,
        "summary": summary
    }

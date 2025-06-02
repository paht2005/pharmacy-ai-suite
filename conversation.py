# Conversational Assistant for Basic Medical Help
# Phat Nguyen Cong - https://github.com/paht2005

from transformers import pipeline

# Set up text generation and summarization tools
try:
    reply_generator = pipeline("text-generation", model="tiiuae/falcon-7b-instruct")
except Exception as load_err:
    raise RuntimeError(f"Couldn't load response model: {load_err}")

try:
    chat_summarizer = pipeline("summarization", model="sshleifer/distilbart-cnn-12-6")
except Exception as load_err:
    raise RuntimeError(f"Summarization model couldn't be loaded: {load_err}")

def respond_to_conversation(history: str) -> str:
    """
    Create a follow-up message for the patient based on previous messages.
    """
    base_prompt = (
        "Act like a virtual health assistant.\n"
        "Given the conversation below, write a helpful follow-up or suggestion.\n\n"
        f"{history}\n\n"
        "Assistant:"
    )
    try:
        result = reply_generator(base_prompt, max_new_tokens=150)
        return result[0]['generated_text'].strip()
    except Exception as gen_err:
        return f"[Issue during response: {gen_err}]"

def get_chat_summary(history: str) -> str:
    """
    Shrink down the conversation to a short summary.
    """
    try:
        result = chat_summarizer(
            history,
            max_length=100,
            min_length=30,
            do_sample=False
        )
        return result[0]['summary_text'].strip()
    except Exception as sum_err:
        return f"[Issue during summarization: {sum_err}]"
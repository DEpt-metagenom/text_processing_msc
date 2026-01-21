import json
import csv
import time
from openai import OpenAI


INPUT_JSON = "patients_sessions_clean.json"
OUTPUT_CSV = "emotion_labels.csv"

MODEL = "gpt-4.1-mini"

LABELS = [
    "neutral",
    "positive",
    "negative",
    "anxiety",
    "anger",
    "sadness",
    "confidence",
    "confusion"
]

client = OpenAI()


def label_utterance(text):
    prompt = f"""
You are given a Hungarian patient utterance from a therapy session.

Choose exactly ONE emotional label from the following list:
{", ".join(LABELS)}

Utterance:
\"\"\"{text}\"\"\"

Respond with ONLY the label.
"""

    response = client.chat.completions.create(
        model=MODEL,
        messages=[
            {"role": "system", "content": "You label emotional tone of text."},
            {"role": "user", "content": prompt}
        ],
        temperature=0
    )

    return response.choices[0].message.content.strip().lower()


def main():
    with open(INPUT_JSON, "r", encoding="utf-8") as f:
        data = json.load(f)

    with open(OUTPUT_CSV, "w", encoding="utf-8", newline="") as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow([
            "patient_id",
            "session",
            "utterance_index",
            "text",
            "emotion_label"
        ])

        for patient_id, sessions in data.items():
            for session_id, utterances in sessions.items():
                for idx, text in enumerate(utterances):
                    label = label_utterance(text)

                    writer.writerow([
                        patient_id,
                        session_id,
                        idx,
                        text,
                        label
                    ])

                    
                    time.sleep(0.3)

    print("Kész:", OUTPUT_CSV)


if __name__ == "__main__":
    main()


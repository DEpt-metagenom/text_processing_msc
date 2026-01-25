import json
import time
from more_itertools import windowed
from openai import OpenAI

# ====== KONFIG ======

INPUT_JSON = "transcript_sessions_dialogue_with_meta.json"
OUTPUT_JSONL = "patient_turn_emotions.jsonl"

MODEL_NAME = "gpt-4.1-mini"
WINDOW_SIZE = 4          # aktuális turn + előző 3
CONTEXT_SIZE = 3         # csak az előzők mennek kontextusba
SLEEP_BETWEEN_CALLS = 0.5

EMOTION_LABELS = [
    "anxiety",
    "anger",
    "sadness",
    "disgust",
    "confusion",
    "confidence",
    "joy",
    "positive",
    "neutral"
]

client = OpenAI()


# ====== PROMPT ======

def build_prompt(context_turns, current_turn):
    context_text = ""

    for t in context_turns:
        if t is None:
            continue
        meta_part = f" (meta: {', '.join(t['meta'])})" if t["meta"] else ""
        context_text += f"- {t['speaker']}: {t['text']}{meta_part}\n"

    current_meta = (
        f" (meta: {', '.join(current_turn['meta'])})"
        if current_turn["meta"] else ""
    )

    prompt = f"""
You are analyzing a psychotherapy-related conversation involving:
- a psychiatrist/therapist,
- a VR avatar representing hallucinated voices,
- a psychiatric patient.

Below are the previous conversation turns (context), followed by the patient's current turn.

Your task:
- Assign EXACTLY ONE emotional label to the patient's current turn.
- Use the context AND the patient's text and meta information.
- If the emotional state is not clear, assign: neutral.
- Only use ONE of the allowed labels listed below.
- Do NOT explain your reasoning.
- Output ONLY the label.

Allowed labels:
{", ".join(EMOTION_LABELS)}

Context:
{context_text}

Patient current turn:
patient: {current_turn['text']}{current_meta}
"""
    return prompt.strip()


# ====== LLM HÍVÁS ======

def annotate_emotion(prompt):
    response = client.responses.create(
        model=MODEL_NAME,
        input=prompt
    )

    label = response.output_text.strip().lower()

    if label not in EMOTION_LABELS:
        return "neutral"

    return label


# ====== SEGÉDFÜGGVÉNY ======

def linearize_session(session):
    turns = []

    # pre_dialogue
    for t in session.get("pre_dialogue", []):
        turns.append(t)

    # dialogs
    for d in session.get("dialogs", []):
        for t in d.get("turns", []):
            turns.append(t)

    return turns


# ====== FŐ FUTÁS ======

def main():
    with open(INPUT_JSON, "r", encoding="utf-8") as f:
        data = json.load(f)

    with open(OUTPUT_JSONL, "w", encoding="utf-8") as out:

        for patient_id, patient in data["patients"].items():
            for session in patient["sessions"]:
                session_id = session["session"]

                all_turns = linearize_session(session)

                indexed_turns = [
                    {**t, "turn_index": i}
                    for i, t in enumerate(all_turns)
                ]

                windows = windowed(
                    indexed_turns,
                    WINDOW_SIZE,
                    fillvalue=None
                )

                for w in windows:
                    current_turn = w[-1]
                    if current_turn is None:
                        continue

                    if current_turn["speaker"] != "patient":
                        continue

                    context = [
                        t for t in w[:-1]
                        if t is not None
                    ]

                    prompt = build_prompt(context, current_turn)
                    emotion = annotate_emotion(prompt)

                    output_record = {
                        "patient_id": patient_id,
                        "session": session_id,
                        "turn_index": current_turn["turn_index"],
                        "text": current_turn["text"],
                        "meta": current_turn["meta"],
                        "emotion_label": emotion,
                        "model": MODEL_NAME,
                        "window_size": WINDOW_SIZE
                    }

                    out.write(
                        json.dumps(output_record, ensure_ascii=False)
                        + "\n"
                    )

                    time.sleep(SLEEP_BETWEEN_CALLS)

    print(f"Kész: {OUTPUT_JSONL}")


if __name__ == "__main__":
    main()


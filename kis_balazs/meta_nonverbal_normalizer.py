import json
import re

INPUT_JSON = "transcript_sessions_dialogue_normalized.json"
OUTPUT_JSON = "transcript_sessions_dialogue_with_meta.json"


PAREN_PATTERN = re.compile(r"\((.*?)\)")
LEADING_DASH_PATTERN = re.compile(r"^[\-\–\—]\s*")


def normalize_turn(turn):
    """
    Egy turn átalakítása:
    - zárójelezett részek -> meta
    - kötőjel eltávolítása a szöveg elejéről
    - text tisztítása
    """

    original_text = turn.get("text", "")
    meta = []

    # 1. Meta kigyűjtése zárójelekből
    meta_matches = PAREN_PATTERN.findall(original_text)
    if meta_matches:
        meta.extend(m.strip() for m in meta_matches if m.strip())

    # 2. Zárójelezett részek eltávolítása a szövegből
    text = PAREN_PATTERN.sub("", original_text)

    # 3. Vezető kötőjel eltávolítása
    text = LEADING_DASH_PATTERN.sub("", text)

    # 4. Whitespace normalizálás
    text = text.strip()

    # 5. Ha nincs verbális szöveg, legyen üres string
    if not text:
        text = ""

    return {
        "speaker": turn["speaker"],
        "text": text,
        "meta": meta
    }


def normalize_all(data):
    for patient in data["patients"].values():
        for session in patient["sessions"]:

            # --- pre_dialogue normalizálása ---
            if "pre_dialogue" in session and session["pre_dialogue"]:
                session["pre_dialogue"] = [
                    normalize_turn(turn)
                    for turn in session["pre_dialogue"]
                ]

            # --- dialogs normalizálása ---
            for dialog in session["dialogs"]:
                dialog["turns"] = [
                    normalize_turn(turn)
                    for turn in dialog["turns"]
                ]

    return data


def main():
    with open(INPUT_JSON, "r", encoding="utf-8") as f:
        data = json.load(f)

    normalized = normalize_all(data)

    with open(OUTPUT_JSON, "w", encoding="utf-8") as f:
        json.dump(normalized, f, ensure_ascii=False, indent=2)

    print(f"Kész: {OUTPUT_JSON}")


if __name__ == "__main__":
    main()


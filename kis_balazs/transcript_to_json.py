from docx import Document
import json
import re
from collections import defaultdict

# =========================
# REGEXEK (VÉDETT)
# =========================

PATIENT_ID_RE = re.compile(r"^\s*(0\d{2})\s*[:,]")
SESSION_RE    = re.compile(r"(\d+)\s*\.\s*ülés", re.IGNORECASE)
DIALOG_RE     = re.compile(r"(\d+)\s*\.\s*dialógus", re.IGNORECASE)
DATE_RE       = re.compile(r"\b\d{4}\.\s*\d{2}\.\s*\d{2}\b")

# =========================
# SEGÉDFÜGGVÉNYEK
# =========================

def extract_patient_id(text):
    m = PATIENT_ID_RE.search(text)
    return m.group(1) if m else None

def extract_session(text):
    m = SESSION_RE.search(text)
    return int(m.group(1)) if m else None

def extract_dialog(text):
    m = DIALOG_RE.search(text)
    return int(m.group(1)) if m else None

def extract_date(text):
    m = DATE_RE.search(text)
    return m.group(0) if m else None

def detect_speaker(paragraph):
    for run in paragraph.runs:
        if run.bold:
            return "avatar"
        if run.underline:
            return "patient"
        if run.italic:
            return "therapist"
    return "unknown"

# =========================
# MAIN
# =========================

def main():
    doc = Document("transcript_javitott.docx")

    patients = defaultdict(lambda: {"sessions": []})

    current_patient = None
    current_session = None
    current_dialog  = None

    for para in doc.paragraphs:
        text = para.text.strip()
        if not text:
            continue

        # -------- PÁCIENS / ÜLÉS METAADAT --------
        pid = extract_patient_id(text)
        sess = extract_session(text)
        date = extract_date(text)
        dial = extract_dialog(text)

        if pid or sess:
            # új páciens
            if pid:
                current_patient = patients[pid]
                current_session = None
                current_dialog  = None

            # új ülés
            if sess is not None:
                current_session = {
                    "session": sess,
                    "date": date,
                    "pre_dialogue": [],
                    "dialogs": []
                }
                current_patient["sessions"].append(current_session)
                current_dialog = None

            # új dialógus (ha jelölve van)
            if dial is not None and current_session is not None:
                current_dialog = {
                    "dialog": dial,
                    "turns": []
                }
                current_session["dialogs"].append(current_dialog)

            continue  # metaadat sort nem dolgozunk fel szövegként

        # -------- SZÖVEGES TARTALOM --------
        if current_session is None:
            continue  # nincs hova tenni → kihagyjuk

        speaker = detect_speaker(para)

        # ha még nincs dialógus → előszöveg
        if current_dialog is None:
            current_session["pre_dialogue"].append({
                "speaker": speaker,
                "text": text
            })
        else:
            current_dialog["turns"].append({
                "speaker": speaker,
                "text": text
            })

    # -------- KIMENTÉS --------
    with open("transcript_corrected.json", "w", encoding="utf-8") as f:
        json.dump({"patients": patients}, f, ensure_ascii=False, indent=2)

    print(f"{len(patients)} páciens feldolgozva.")

# =========================

if __name__ == "__main__":
    main()


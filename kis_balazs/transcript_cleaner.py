import json
from collections import defaultdict

INPUT_FILE = "transcript_corrected.json"
OUTPUT_FILE = "patients_sessions_clean.json"


def main():
    with open(INPUT_FILE, "r", encoding="utf-8") as f:
        data = json.load(f)

    result = {}

    for patient_id, patient_data in data["patients"].items():
        sessions_out = defaultdict(list)

        for session_block in patient_data["sessions"]:
            session_number = str(session_block["session"])

            # ---- pre_dialogue ----
            for turn in session_block.get("pre_dialogue", []):
                if turn.get("speaker") == "patient":
                    sessions_out[session_number].append(
                        turn.get("text", "").strip()
                    )

            # ---- dialogs ----
            for dialog in session_block.get("dialogs", []):
                for turn in dialog.get("turns", []):
                    if turn.get("speaker") == "patient":
                        sessions_out[session_number].append(
                            turn.get("text", "").strip()
                        )

        # csak akkor vesszük fel, ha van tényleges páciens szöveg
        if sessions_out:
            result[patient_id] = dict(sessions_out)

    with open(OUTPUT_FILE, "w", encoding="utf-8") as f:
        json.dump(result, f, ensure_ascii=False, indent=2)

    print(f"Kész. Feldolgozott páciensek: {len(result)}")
    print(f"Kimeneti fájl: {OUTPUT_FILE}")


if __name__ == "__main__":
    main()


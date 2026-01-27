import json
from collections import OrderedDict

INPUT_JSON = "transcript_corrected.json"
OUTPUT_JSON = "transcript_sessions_merged.json"


def merge_sessions_for_patient(sessions):
    """
    Azonos session számú session-objektumok összevonása,
    az eredeti sorrend megtartásával.
    """

    merged = OrderedDict()

    for session_obj in sessions:
        session_id = session_obj["session"]

        if session_id not in merged:
            merged[session_id] = {
                "session": session_id,
                "date": session_obj.get("date"),
                "pre_dialogue": [],
                "dialogs": []
            }

        # pre_dialogue hozzáfűzése (ha van)
        if session_obj.get("pre_dialogue"):
            merged[session_id]["pre_dialogue"].extend(
                session_obj["pre_dialogue"]
            )

        # dialogs hozzáfűzése (ha van)
        if session_obj.get("dialogs"):
            merged[session_id]["dialogs"].extend(
                session_obj["dialogs"]
            )

    # dialógusok újraszámozása session-ön belül
    for session in merged.values():
        for idx, dialog in enumerate(session["dialogs"], start=1):
            dialog["dialog"] = idx

    return list(merged.values())


def merge_all_patients(data):
    output = {"patients": {}}

    for patient_id, patient_data in data["patients"].items():
        sessions = patient_data["sessions"]
        merged_sessions = merge_sessions_for_patient(sessions)

        output["patients"][patient_id] = {
            "sessions": merged_sessions
        }

    return output


def main():
    with open(INPUT_JSON, "r", encoding="utf-8") as f:
        data = json.load(f)

    merged_data = merge_all_patients(data)

    with open(OUTPUT_JSON, "w", encoding="utf-8") as f:
        json.dump(merged_data, f, ensure_ascii=False, indent=2)

    print(f"Kész: {OUTPUT_JSON}")


if __name__ == "__main__":
    main()


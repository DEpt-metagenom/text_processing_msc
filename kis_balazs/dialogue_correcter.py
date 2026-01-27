import json

INPUT_JSON = "transcript_sessions_merged.json"
OUTPUT_JSON = "transcript_sessions_dialogue_normalized.json"


def normalize_session(session):
    """
    Ha egy session-nek nincs dialógusa, de van pre_dialogue-ja,
    akkor a pre_dialogue-ból csinálunk egyetlen (1.) dialógust.
    """

    if session.get("dialogs") and len(session["dialogs"]) > 0:
        # már dialógus-alapú, nem nyúlunk hozzá
        return session

    pre = session.get("pre_dialogue", [])

    if not pre:
        # se pre, se dialog – üres session (ritka, de hagyjuk békén)
        return session

    # pre_dialogue → 1. dialog
    session["dialogs"] = [
        {
            "dialog": 1,
            "turns": pre
        }
    ]

    session["pre_dialogue"] = []

    return session


def normalize_all(data):
    for patient in data["patients"].values():
        patient["sessions"] = [
            normalize_session(session)
            for session in patient["sessions"]
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


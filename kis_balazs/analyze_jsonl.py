import json
import pandas as pd
from collections import Counter

INPUT_JSONL = "patient_turn_emotions.jsonl"

NEGATIVE = {"anxiety", "anger", "sadness", "disgust", "confusion"}
POSITIVE = {"confidence", "joy", "positive"}
NEUTRAL = {"neutral"}


def load_jsonl(path):
    rows = []
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            rows.append(json.loads(line))
    return pd.DataFrame(rows)


def emotion_group(label):
    if label in NEGATIVE:
        return "negative"
    if label in POSITIVE:
        return "positive"
    return "neutral"


def main():
    df = load_jsonl(INPUT_JSONL)

    # extra oszlop: érzelemcsoport
    df["emotion_group"] = df["emotion_label"].apply(emotion_group)

    print("\n=== Globális érzelmi megoszlás ===")
    print(df["emotion_label"].value_counts(normalize=True).round(3))

    print("\n=== Páciensenként, session-önkénti összesítés ===")

    for patient_id in sorted(df["patient_id"].unique()):
        print(f"\n--- Páciens {patient_id} ---")

        patient_df = df[df["patient_id"] == patient_id]

        for session in sorted(patient_df["session"].unique()):
            session_df = patient_df[patient_df["session"] == session]

            total = len(session_df)
            label_counts = Counter(session_df["emotion_label"])
            group_counts = Counter(session_df["emotion_group"])

            print(f"\nSession {session} (n={total})")

            print(" Érzelmi címkék:")
            for label, cnt in label_counts.most_common():
                print(f"  {label:12s}: {cnt/total:.2%}")

            print(" Érzelemcsoportok:")
            for group, cnt in group_counts.items():
                print(f"  {group:8s}: {cnt/total:.2%}")

    print("\n=== Időbeli tendencia (negatív arány) ===")

    trend_rows = []

    for (patient_id, session), g in df.groupby(["patient_id", "session"]):
        neg_ratio = (g["emotion_group"] == "negative").mean()
        trend_rows.append({
            "patient_id": patient_id,
            "session": session,
            "negative_ratio": neg_ratio
        })

    trend_df = pd.DataFrame(trend_rows)
    print(trend_df.sort_values(["patient_id", "session"]))


if __name__ == "__main__":
    main()


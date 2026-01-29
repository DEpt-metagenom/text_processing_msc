import os
import xml.etree.ElementTree as ET
import ollama
from transformers import AutoTokenizer

# ==============================
# BEÁLLÍTÁSOK
# ==============================


MODEL = "qwen2.5:7b"
TOKENIZER_NAME = "Qwen/Qwen2.5-32B-Instruct"

XML_DIR = r"C:\Users\szk19\OneDrive\Asztali gép\Bioinformatika I félév\Új mappa\Szakirodalmi szövegbányászat\pmc_xml"
OUT_DIR = "myplace/summaries/qwen2.5_32b/"

MAX_CONTEXT_TOKENS = 8000        # biztonságos chunk méret
SUMMARY_TEMPERATURE = 0.2

# ==============================
# PROMPT
# ==============================

SYSTEM_PROMPT_PRE = (
    "Az alábbi szöveg egy kutatási interjú vagy dokumentum XML-ből kinyert tartalma.\n"
    "Feladatod, hogy tömören és strukturáltan összefoglald a lényeget.\n\n"
    "Kérlek, az összefoglalás tartalmazza:\n"
    "• fő témák\n"
    "• visszatérő minták\n"
    "• fontos megállapítások\n\n"
    "<szöveg>\n"
)

SYSTEM_PROMPT_POST = "\n</szöveg>"

# ==============================
# KÖNYVTÁRAK ELŐKÉSZÍTÉSE
# ==============================

os.makedirs(OUT_DIR, exist_ok=True)

# ==============================
# TOKENIZER
# ==============================

tokenizer = AutoTokenizer.from_pretrained(TOKENIZER_NAME)

def count_tokens(text: str) -> int:
    return len(tokenizer.encode(text))

# ==============================
# XML SZÖVEG KINYERÉS
# ==============================

def extract_text_from_xml(xml_path: str) -> str:
    tree = ET.parse(xml_path)
    root = tree.getroot()

    texts = []

    for elem in root.iter():
        if elem.text and elem.text.strip():
            texts.append(elem.text.strip())

    return "\n".join(texts)

# ==============================
# SZÖVEG DARABOLÁS
# ==============================

def chunk_text(text: str, max_tokens: int):
    words = text.split()
    chunks = []
    current = []

    for word in words:
        current.append(word)
        if count_tokens(" ".join(current)) >= max_tokens:
            chunks.append(" ".join(current))
            current = []

    if current:
        chunks.append(" ".join(current))

    return chunks

# ==============================
# OLLAMA ÖSSZEFOGLALÁS
# ==============================

def summarize_text(text: str) -> str:
    prompt = SYSTEM_PROMPT_PRE + text + SYSTEM_PROMPT_POST

    response = ollama.generate(
        model=MODEL,
        prompt=prompt,
        options={
            "temperature": SUMMARY_TEMPERATURE
        }
    )

    return response["response"]

# ==============================
# FŐ PROGRAM
# ==============================

def main():
    xml_files = sorted([f for f in os.listdir(XML_DIR) if f.endswith(".xml")])

    print(f"{len(xml_files)} XML fájl feldolgozása indul...\n")

    for i, fname in enumerate(xml_files, start=1):
        print(f"[{i}/{len(xml_files)}] {fname}")

        xml_path = os.path.join(XML_DIR, fname)
        text = extract_text_from_xml(xml_path)

        if not text.strip():
            print("  ⚠ Üres fájl, kihagyva\n")
            continue

        chunks = chunk_text(text, MAX_CONTEXT_TOKENS)
        print(f"  → {len(chunks)} részre bontva")

        partial_summaries = []

        for idx, chunk in enumerate(chunks, start=1):
            print(f"    • Rész {idx}/{len(chunks)} összefoglalása...")
            summary = summarize_text(chunk)
            partial_summaries.append(summary)

        # Ha több rész volt, készítünk egy végső összefoglalót
        if len(partial_summaries) > 1:
            print("  → Végső összefoglaló készítése...")
            combined = "\n".join(partial_summaries)
            final_summary = summarize_text(combined)
        else:
            final_summary = partial_summaries[0]

        out_file = fname.replace(".xml", "_summary.txt")
        out_path = os.path.join(OUT_DIR, out_file)

        with open(out_path, "w", encoding="utf-8") as f:
            f.write(final_summary)

        print(f"  ✔ Mentve: {out_file}\n")

    print("✅ Feldolgozás kész.")

# ==============================
# FUTTATÁS
# ==============================

if __name__ == "__main__":
    main()

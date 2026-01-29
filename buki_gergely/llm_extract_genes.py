import json
import time
from pathlib import Path
from collections import Counter

import requests
from bs4 import BeautifulSoup

# =========================
# CONFIG
# =========================
XML_PATH = Path("pubmed_articles.xml")
MODEL = "gemma3:12b"
OLLAMA_URL = "http://127.0.0.1:11434/api/generate"

OUT_JSON = Path("pmid_to_genes.json")
OUT_CSV = Path("gene_frequencies.csv")

CHECKPOINT_EVERY = 20
SLEEP_SEC = 0.2   # ha túl gyors, emeld 0.3–0.5-re

PROMPT_TEMPLATE = """You are given the title and abstract of a scientific article.

TASK:
Extract the gene symbol(s) that the abstract explicitly suggests as causal for the disorder/disease being studied.

RULES:
- Output ONLY official gene symbols that APPEAR VERBATIM in the provided text.
- Do NOT infer or guess genes that are not explicitly named.
- If multiple causal genes are explicitly mentioned, output a JSON array of gene symbols.
- If no causal gene is suggested, output: []
- Do NOT include any other text, explanations, keys, punctuation, or markdown.

TITLE:
{title}

ABSTRACT:
{abstract}
"""

# =========================
# FUNCTIONS
# =========================
def parse_pubmed_xml(xml_path: Path):
    with open(xml_path, "r", encoding="utf-8") as f:
        soup = BeautifulSoup(f, features="xml")

    articles = []
    for art in soup.find_all("PubmedArticle"):
        pmid_tag = art.find("PMID")
        if not pmid_tag:
            continue
        pmid = pmid_tag.get_text(strip=True)

        title_tag = art.find("ArticleTitle")
        title = title_tag.get_text(" ", strip=True) if title_tag else ""

        abstract_tag = art.find("Abstract")
        if abstract_tag:
            abstract = " ".join(
                t.get_text(" ", strip=True)
                for t in abstract_tag.find_all("AbstractText")
            )
        else:
            abstract = ""

        if title or abstract:
            articles.append({
                "pmid": pmid,
                "title": title,
                "abstract": abstract
            })

    return articles


def call_ollama(session: requests.Session, title: str, abstract: str):
    prompt = PROMPT_TEMPLATE.format(title=title, abstract=abstract)

    payload = {
        "model": MODEL,
        "prompt": prompt,
        "stream": False,
        "options": {
            "temperature": 0,
            "num_predict": 64,
            "num_ctx": 2048
        }
    }

    last_error = None
    for attempt in range(3):
        try:
            r = session.post(OLLAMA_URL, json=payload, timeout=180)

            if r.status_code >= 500:
                last_error = f"HTTP {r.status_code}"
                time.sleep(2 * (attempt + 1))
                continue

            r.raise_for_status()
            text = (r.json().get("response") or "").strip()

            try:
                genes = json.loads(text)
                if isinstance(genes, list):
                    return [g.strip() for g in genes if isinstance(g, str) and g.strip()]
                return []
            except json.JSONDecodeError:
                if text == "[]":
                    return []
                if text:
                    return [text.strip()]
                return []

        except Exception as e:
            last_error = str(e)
            time.sleep(2 * (attempt + 1))

    raise RuntimeError(f"Ollama failed after retries: {last_error}")


def load_existing(path: Path):
    if path.exists():
        try:
            return json.loads(path.read_text(encoding="utf-8"))
        except Exception:
            return {}
    return {}


def save_json(obj, path: Path):
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=False), encoding="utf-8")


def write_gene_frequencies(pmid_to_genes: dict, out_csv: Path):
    counter = Counter()
    for genes in pmid_to_genes.values():
        for g in genes:
            counter[g] += 1

    with open(out_csv, "w", encoding="utf-8") as f:
        f.write("gene,count\n")
        for gene, count in counter.most_common():
            f.write(f"{gene},{count}\n")


# =========================
# MAIN
# =========================
def main():
    if not XML_PATH.exists():
        raise FileNotFoundError(f"Missing XML file: {XML_PATH}")

    articles = parse_pubmed_xml(XML_PATH)
    print(f"Parsed {len(articles)} PubMed articles")

    pmid_to_genes = load_existing(OUT_JSON)
    done = set(pmid_to_genes.keys())

    session = requests.Session()
    session.trust_env = False  # ignore proxy completely

    new = 0
    total = len(articles)

    for idx, art in enumerate(articles, start=1):
        pmid = art["pmid"]
        if pmid in done:
            continue

        genes = call_ollama(session, art["title"], art["abstract"])
        pmid_to_genes[pmid] = genes
        new += 1

        print(f"{idx}/{total} PMID {pmid} -> {genes}")

        if new % CHECKPOINT_EVERY == 0:
            save_json(pmid_to_genes, OUT_JSON)
            print(f"Checkpoint saved ({len(pmid_to_genes)} records)")

        time.sleep(SLEEP_SEC)

    save_json(pmid_to_genes, OUT_JSON)
    write_gene_frequencies(pmid_to_genes, OUT_CSV)

    print("DONE")
    print(f"Results: {OUT_JSON}")
    print(f"Gene counts: {OUT_CSV}")


if __name__ == "__main__":
    main()

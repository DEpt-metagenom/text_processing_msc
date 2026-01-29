import os
import time
from urllib.parse import quote
from urllib.request import urlopen
from bs4 import BeautifulSoup

# ─── NCBI E-utilities beállítások ──────────────────────────────────────────────
SEARCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
FETCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"

DB_PUBMED = "pubmed"
DB_PMC = "pmc"

MAX_RESULTS = 50
DELAY = 0.34   # < 3 kérés / másodperc (NCBI ajánlás)

TOOL = "pmc_fulltext_xml_downloader"
EMAIL = "szucskrisztina6@gmail.com"   # ← SAJÁT EMAIL


# ─── 1. PubMed keresés (free full text) ───────────────────────────────────────
def pubmed_search(term, max_results=MAX_RESULTS):
    term = f"({term}) AND free full text[sb]"

    query = (
        f"{SEARCH_URL}?db={DB_PUBMED}"
        f"&retmax={max_results}"
        f"&term={quote(term)}"
        f"&tool={TOOL}"
        f"&email={EMAIL}"
    )

    with urlopen(query, timeout=20) as r:
        soup = BeautifulSoup(r, "lxml-xml")

    return [id_tag.text for id_tag in soup.find_all("Id")]


# ─── 2. PMID → PMCID kinyerése ────────────────────────────────────────────────
def pmid_to_pmcid(pmid):
    url = (
        f"{FETCH_URL}?db={DB_PUBMED}"
        f"&id={pmid}"
        f"&retmode=xml"
        f"&tool={TOOL}"
        f"&email={EMAIL}"
    )

    with urlopen(url, timeout=20) as r:
        soup = BeautifulSoup(r, "lxml-xml")

    pmc_id = soup.find("ArticleId", {"IdType": "pmc"})
    if pmc_id:
        return pmc_id.text.strip()

    return None


# ─── 3. PMC efetch → TELJES CIKK XML ──────────────────────────────────────────
def fetch_pmc_fulltext_xml(pmcid, out_dir="pmc_xml"):
    if not pmcid:
        return False

    os.makedirs(out_dir, exist_ok=True)

    url = (
        f"{FETCH_URL}?db={DB_PMC}"
        f"&id={pmcid}"
        f"&retmode=xml"
        f"&tool={TOOL}"
        f"&email={EMAIL}"
    )

    try:
        with urlopen(url, timeout=30) as r:
            xml = r.read()

        if b"<article" not in xml.lower():
            return False

        out_path = os.path.join(out_dir, f"{pmcid}.xml")
        with open(out_path, "wb") as f:
            f.write(xml)

        return True

    except Exception:
        return False


# ─── 4. Teljes pipeline ───────────────────────────────────────────────────────
def download_pmc_fulltexts(query):
    pmids = pubmed_search(query)
    print(f"{len(pmids)} PubMed találat")

    downloaded = 0

    for i, pmid in enumerate(pmids, 1):
        print(f"{i}/{len(pmids)} → PMID {pmid}")

        pmcid = pmid_to_pmcid(pmid)

        if pmcid:
            if fetch_pmc_fulltext_xml(pmcid):
                downloaded += 1
                print(f"   ✔ letöltve: {pmcid}")
            else:
                print(f"   ✖ PMC XML nem elérhető")
        else:
            print(f"   ✖ nincs PMCID")

        time.sleep(DELAY)

    print(f"\nKész. Letöltött teljes cikkek: {downloaded}")
    print("Mappa: ./pmc_xml/")


# ─── Futtatás ────────────────────────────────────────────────────────────────
if __name__ == "__main__":
    download_pmc_fulltexts(
        "BRCA tumor suppressor gene role in breast cancer"
    )

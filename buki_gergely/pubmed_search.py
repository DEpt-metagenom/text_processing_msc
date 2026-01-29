import os
from Bio import Entrez
import time

# Set your email for NCBI API access
Entrez.email = "buki.gergely@pte.hu"

def search_pubmed_articles(query, max_results=100):
    """Searches PubMed for articles and returns their PubMed IDs (PMIDs)."""
    try:
        handle = Entrez.esearch(db="pubmed", term=query, retmax=max_results)
        record = Entrez.read(handle)
        handle.close()
        return record.get("IdList", [])
    except Exception as e:
        print(f"Error during search: {e}")
        return []

def download_xml_records(article_ids, db="pubmed"):
    """Downloads XML records from NCBI EFetch given a list of PMIDs."""
    if not article_ids:
        print("No article IDs provided.")
        return None

    try:
        time.sleep(1)  # Avoid NCBI rate limits
        handle = Entrez.efetch(db=db, id=article_ids, rettype="xml", retmode="text")
        xml_data = handle.read()
        handle.close()

        xml_text = xml_data.decode("utf-8") if isinstance(xml_data, bytes) else xml_data
        print("Downloaded XML is None?", xml_text is None)
        print("Downloaded XML length:", len(xml_text) if xml_text else 0)

        return xml_text
    except Exception as e:
        print(f"Error downloading XML: {e}")
        return None

if __name__ == "__main__":
    print("Current working directory:", os.getcwd())

    # Step 1: Search for specific speech disorder related articles in PubMed
    query = (
        '("speech disorder"[Title/Abstract] OR "speech apraxia"[Title/Abstract] OR "childhood apraxia of speech"[Title/Abstract] '
        'OR "language impairment"[Title/Abstract]) '
        'AND (gene[Title/Abstract] OR genetic[Title/Abstract] OR mutation[Title/Abstract] OR variant[Title/Abstract])'
    )

    pubmed_ids = search_pubmed_articles(query, max_results=50)

    if pubmed_ids:
        print(f"Found {len(pubmed_ids)} articles.")
        print("First 5 PMIDs:", pubmed_ids[:5])

        # Step 2: Download XML records for PubMed IDs
        xml_data = download_xml_records(pubmed_ids)

        if xml_data:
            out_path = os.path.abspath("pubmed_articles.xml")
            with open(out_path, "w", encoding="utf-8") as f:
                f.write(xml_data)
            print("PubMed XML records saved to:", out_path)
        else:
            print("xml_data is empty/None, file not written.")
    else:
        print("No PubMed articles found.")

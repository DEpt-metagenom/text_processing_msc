import requests
from Bio import Entrez
import time
import xml.etree.ElementTree as ET

# Set your email for NCBI API access
Entrez.email = "your_email@example.com"  # Replace with your actual email

def search_pubmed_articles(keywords, max_results=100):
    """Searches PubMed for articles and returns their PubMed IDs (PMIDs)."""
    try:
        search_term = " AND ".join(keywords)
        handle = Entrez.esearch(db="pubmed", term=search_term, retmax=max_results)
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
        return xml_data.decode("utf-8") if isinstance(xml_data, bytes) else xml_data
    except Exception as e:
        print(f"Error downloading XML: {e}")
        return None

def extract_main_article_pmc_ids(xml_content):
    """Extracts PMC IDs from the main article (not references) in the XML."""
    pmc_ids = []
    try:
        root = ET.fromstring(xml_content)
        for pubmed_data in root.findall(".//PubmedData"):
            article_id_list = pubmed_data.find("ArticleIdList")
            if article_id_list:
                for article_id in article_id_list.findall("ArticleId"):
                    if article_id.attrib.get("IdType") == "pmc":
                        pmc_ids.append(article_id.text.replace("PMC", ""))  # Remove 'PMC' prefix
    except ET.ParseError as e:
        print(f"XML parsing error: {e}")
    return pmc_ids

def fetch_full_text_from_pmc(pmcid):
    """Fetches full-text XML from PMC via OAI-PMH API."""
    try:
        url = f"https://www.ncbi.nlm.nih.gov/pmc/oai/oai.cgi?verb=GetRecord&metadataPrefix=pmc&identifier=oai:pubmedcentral.nih.gov:{pmcid}"
        response = requests.get(url)

        if response.status_code == 200:
            return response.text  # Return full text XML
        else:
            print(f"Failed to fetch full text for PMC ID {pmcid} (Status: {response.status_code})")
            return None
    except Exception as e:
        print(f"Error fetching full text for PMC ID {pmcid}: {e}")
        return None

if __name__ == "__main__":
    # Step 1: Search for articles in PubMed
    keywords = ["wgs", "mycobacterium", "nanopore"]
    pubmed_ids = search_pubmed_articles(keywords, max_results=50)

    if pubmed_ids:
        print(f"Found {len(pubmed_ids)} articles: {', '.join(pubmed_ids)}")

        # Step 2: Download XML records for PubMed IDs
        xml_data = download_xml_records(pubmed_ids)

        if xml_data:
            # Save PubMed XML for reference
            with open("pubmed_articles.xml", "w", encoding="utf-8") as f:
                f.write(xml_data)
            print("PubMed XML records saved.")

            # Step 3: Extract PMC IDs (only from main articles)
            pmc_ids = extract_main_article_pmc_ids(xml_data)

            if pmc_ids:
                print(f"Extracted {len(pmc_ids)} PMC IDs: {', '.join(pmc_ids)}")

                # Step 4: Fetch full-text articles from PMC
                full_texts = []
                for pmcid in pmc_ids:
                    print(f"Fetching full-text for PMC ID: {pmcid}...")
                    full_text = fetch_full_text_from_pmc(pmcid)
                    if full_text:
                        full_texts.append(f"<h2>PMC{pmcid}</h2>\n{full_text}")
                        time.sleep(1)  # Respect NCBI rate limits

                # Step 5: Save all full-text articles to an HTML file
                if full_texts:
                    with open("full_texts.html", "w", encoding="utf-8") as f:
                        f.write("<html><body>\n" + "\n<hr>\n".join(full_texts) + "\n</body></html>")
                    print("Full-text articles saved to full_texts.html")
                else:
                    print("No full-text articles were retrieved.")
            else:
                print("No PMC IDs found.")
        else:
            print("Failed to download PubMed XML records.")
    else:
        print("No PubMed articles found.")

!pip install biopython
from Bio import Entrez
import time

def search_pubmed_central(keywords, max_results=100):
    """
    Searches PubMed for articles based on keywords using the ESearch API.

    Args:
        keywords (list): A list of keywords to search for.
        max_results (int): The maximum number of results to retrieve.

    Returns:
        list: A list of PubMed IDs of matching articles. Returns an empty list on error.
    """
    Entrez.email = "schmidtviki15@gmail.com"  # Replace with your actual email

    try:
        search_term = " AND ".join(keywords)  # Join keywords with AND
        handle = Entrez.esearch(db="pubmed", term=search_term, retmax=max_results)
        record = Entrez.read(handle)
        handle.close()

        pubmed_ids = record.get("IdList", [])  # Safely get IdList
        return pubmed_ids

    except Exception as e:
        print(f"An error occurred: {e}")
        return []  # Return an empty list on error.

if __name__ == "__main__":
    keywords = ["wgs", "mycobacterium", "nanopore"]

    # Search PubMed Central
    pubmed_ids = search_pubmed_central(keywords)

    if pubmed_ids:
        print(f"Found {len(pubmed_ids)} articles with the specified keywords:")
        for pubmed_id in pubmed_ids:
            print(pubmed_id)
    else:
        print("No articles found with the specified keywords.")
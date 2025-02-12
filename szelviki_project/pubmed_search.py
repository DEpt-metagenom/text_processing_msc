pip install requests beautifulsoup4
import requests
from bs4 import BeautifulSoup
from time import sleep
import re

SEARCH_URL = "https://pubmed.ncbi.nlm.nih.gov/?term={}&sort=relevance"
ARTICLE_URL = "https://pubmed.ncbi.nlm.nih.gov/{}"
RESULTS_MAX = 20  # Number of results to fetch
DELAY = 1  # Delay between requests

KEYWORDS = ["wgs", "mycobacterium", "nanopore"]

def get_pubmed_results():
    search_query = "+".join(KEYWORDS)
    url = SEARCH_URL.format(search_query)
    response = requests.get(url)
    if response.status_code != 200:
        print("Failed to fetch PubMed search results.")
        return []
    
    soup = BeautifulSoup(response.text, "html.parser")
    article_links = [a["href"] for a in soup.select(".docsum-title")][:RESULTS_MAX]
    article_ids = [re.findall(r'\d+', link)[0] for link in article_links]
    return article_ids


def get_article_details(article_id):
    url = ARTICLE_URL.format(article_id)
    response = requests.get(url)
    if response.status_code != 200:
        print(f"Failed to fetch article {article_id}.")
        return None
    
    soup = BeautifulSoup(response.text, "html.parser")
    title = soup.select_one("h1.heading-title").text.strip()
    doi_link = soup.select_one("a[href*='doi.org']")
    doi_url = doi_link["href"] if doi_link else None
    
    full_text_div = soup.find("div", class_="abstract-content")
    full_text = full_text_div.text.strip() if full_text_div else "No abstract available"
    
    return {
        "id": article_id,
        "title": title,
        "full_text": full_text,
        "doi_url": doi_url
    }


def get_methods_from_doi(doi_url):
    if not doi_url:
        return "DOI not available"
    
    response = requests.get(doi_url)
    if response.status_code != 200:
        return "Failed to fetch DOI page"
    
    soup = BeautifulSoup(response.text, "html.parser")
    
    possible_sections = ["methods", "methodology", "materials and methods"]
    for section in possible_sections:
        methods_section = soup.find(lambda tag: tag.name in ["h2", "h3"] and section in tag.text.lower())
        if methods_section:
            methods_content = []
            sibling = methods_section.find_next_sibling()
            while sibling and sibling.name not in ["h2", "h3"]:  # Stop at the next major section
                methods_content.append(sibling.text)
                sibling = sibling.find_next_sibling()
            return "\n".join(methods_content) if methods_content else "Methods section not found"
    
    return "Methods section not found"


def main():
    article_ids = get_pubmed_results()
    for article_id in article_ids:
        details = get_article_details(article_id)
        if not details:
            continue
        
        methods_text = get_methods_from_doi(details["doi_url"]) if details["doi_url"] else "No DOI available"
        
        print("\n=== ARTICLE ===")
        print(f"Title: {details['title']}")
        print(f"PubMed ID: {details['id']}")
        print(f"Full Text: {details['full_text'][:500]}...")  # Print first 500 chars
        print(f"Methods: {methods_text[:500]}...")  # Print first 500 chars
        print(f"DOI: {details['doi_url']}")
        sleep(DELAY)


if __name__ == "__main__":
    main()

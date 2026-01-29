from bs4 import BeautifulSoup

with open("pubmed_articles.xml") as f:
    soup = BeautifulSoup(f, "lxml")

articles = []

for art in soup.find_all("pubmedarticle"):
    pmid = art.find("pmid").text if art.find("pmid") else None

    title_tag = art.find("articletitle")
    title = title_tag.text.strip() if title_tag else ""
    
    abstract_tag = art.find("abstract")
    abstract = " ".join([t.text for t in abstract_tag.find_all("abstracttext")]) if abstract_tag else ""

    articles.append({
        "pmid": pmid,
        "title": title,
        "abstract": abstract
    })

print(f"Parsed {len(articles)} articles")
print("First article example:")
print(articles[0])
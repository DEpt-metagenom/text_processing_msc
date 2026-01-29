import csv
import re
import time
from Bio import Entrez

def osszes_pubmed_cikk_mentese(filenev="deltaF508_osszes_absztrakt.csv"):
    # 1. BEÁLLÍTÁSOK - EMAIL KÖTELEZŐ!
    Entrez.email = "csorba.gabriella@med.unideb.hu.hu"  # IDE ÍRD A SAJÁT EMAILED
    
    # Keresési feltételek
    query = '(cystic fibrosis[Title/Abstract]) AND (deltaF508[Title/Abstract] OR F508del[Title/Abstract] OR "ΔF508"[Title/Abstract])'
    mutacio_minta = r"delta\s*F508|F508\s*del|ΔF508|p\.?Phe508del|Phe508del"
    
    print(f"Keresés: {query}")
    print("Kezdődik a letöltés...")

    # 2. KERESÉS
    try:
        handle = Entrez.esearch(db="pubmed", term=query, usehistory="y", retmax=100000)
        search_results = Entrez.read(handle)
        handle.close()
    except Exception as e:
        print(f"Hiba a keresés során: {e}")
        return

    count = int(search_results["Count"])
    webenv = search_results["WebEnv"]
    query_key = search_results["QueryKey"]

    print(f"Összesen {count} találat.")
    if count == 0:
        return

    # 3. LETÖLTÉS
    batch_size = 200
    talalatok_szama = 0
    sikeres_batch = 0

    with open(filenev, mode='w', newline='', encoding='utf-8') as csv_file:
        writer = csv.writer(csv_file)
        writer.writerow(['PMID', 'Cím', 'Link', 'Absztrakt', 'Év', 'Szerzők', 'Folyóirat'])

        for start in range(0, count, batch_size):
            print(f"Batch {start//batch_size + 1}/{(count-1)//batch_size + 1}: {start}-{min(start+batch_size, count)}")
            
            try:
                fetch_handle = Entrez.efetch(
                    db="pubmed",
                    retstart=start,
                    retmax=batch_size,
                    webenv=webenv,
                    query_key=query_key,
                    rettype="abstract",
                    retmode="xml"
                )
                
                articles = Entrez.read(fetch_handle)
                fetch_handle.close()
                sikeres_batch += 1
                
            except Exception as e:
                print(f"  Hiba: {e}, kihagyom ezt a köteget")
                time.sleep(3)
                continue

            # FELDOLGOZÁS
            batch_talalat = 0
            for article in articles['PubmedArticle']:
                pmid = str(article['MedlineCitation']['PMID'])
                title = article['MedlineCitation']['Article']['ArticleTitle']
                
                # Absztrakt
                abstract_parts = article['MedlineCitation']['Article'].get('Abstract', {}).get('AbstractText', [])
                abstract_text = " ".join([str(resz) for resz in abstract_parts])
                
                # Metaadatok
                pub_date = article['MedlineCitation']['Article']['Journal']['JournalIssue'].get('PubDate', {})
                year = pub_date.get('Year', pub_date.get('MedlineDate', 'Ismeretlen')).split()[0]  # Csak év
                
                authors = article['MedlineCitation']['Article'].get('AuthorList', [])
                author_names = "; ".join([f"{author.get('LastName', '')} {author.get('ForeName', '')[:1]}." 
                                        for author in authors[:3]])  # Első 3 szerző
                
                journal = article['MedlineCitation']['Article']['Journal']['Title']
                link = f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/"
                
                # Szűrés a mutációra
                if re.search(mutacio_minta, abstract_text, re.IGNORECASE):
                    writer.writerow([pmid, title, link, abstract_text, year, author_names, journal])
                    talalatok_szama += 1
                    batch_talalat += 1
            
            print(f"  Ebben a batchben: {batch_talalat} új találat")
            time.sleep(0.5)  # Szünet a szerverek védelme érdekében

    print(f"\nKÉSZ! {talalatok_szama} releváns cikk mentve '{filenev}' fájlba.")
    print(f"{sikeres_batch} batch sikeresen lezárva.")

# Futtatás
if __name__ == "__main__":
    osszes_pubmed_cikk_mentese()
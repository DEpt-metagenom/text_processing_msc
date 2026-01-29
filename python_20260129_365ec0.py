import re
import pandas as pd
from typing import List, Dict
from process_text_files import process_text_files  # Az előző fájlból importáljuk

def analyze_txt_mutations(directory: str = "text_files") -> None:
    """
    Fő analízis függvény txt fájlokra.
    """
    
    print("="*60)
    print("DELTAF508 ÉS CFTRDELE2,3 MUTÁCIÓK KERESÉSE TXT FÁJLOKBAN")
    print("="*60)
    
    # Fájlok feldolgozása
    data = process_text_files(directory)
    
    if not data:
        print("Nincs feldolgozható adat!")
        return
    
    # Statisztikák számítása
    total_files = len(data)
    both_mutations = sum(1 for item in data if item['has_both'])
    only_deltaF508 = sum(1 for item in data if item['deltaF508_count'] > 0 and item['cftr_dele23_count'] == 0)
    only_cftr_dele23 = sum(1 for item in data if item['cftr_dele23_count'] > 0 and item['deltaF508_count'] == 0)
    neither = sum(1 for item in data if item['deltaF508_count'] == 0 and item['cftr_dele23_count'] == 0)
    
    # Eredmények kiírása
    print(f"\nÖsszes feldolgozott fájl: {total_files}")
    print(f"Két mutáció együtt: {both_mutations}")
    print(f"Csak deltaF508: {only_deltaF508}")
    print(f"Csak CFTRdele2,3: {only_cftr_dele23}")
    print(f"Egyik mutáció sem: {neither}")
    
    # Mindkét mutációt tartalmazó fájlok
    both_files = [item for item in data if item['has_both']]
    
    if both_files:
        print(f"\nKét mutációt tartalmazó fájlok ({len(both_files)} db):")
        for i, item in enumerate(both_files[:10], 1):
            print(f"{i}. {item['filename']}")
            print(f"   deltaF508: {item['deltaF508_count']}, CFTRdele2,3: {item['cftr_dele23_count']}")
        
        if len(both_files) > 10:
            print(f"... és még {len(both_files) - 10} további fájl")
        
        # Eredmények mentése
        df = pd.DataFrame(both_files)
        df.to_csv('mutation_results.csv', index=False, encoding='utf-8')
        print(f"\nEredmények elmentve: mutation_results.csv")

if __name__ == "__main__":
    analyze_txt_mutations()
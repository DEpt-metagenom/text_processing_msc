import os
import re
from typing import List, Dict

def process_text_files(directory: str = "text_files") -> List[Dict]:
    """
    Szöveges fájlok feldolgozása a megadott könyvtárból.
    
    Args:
        directory: A txt fájlokat tartalmazó könyvtár
        
    Returns:
        Lista a feldolgozott adatokkal
    """
    
    if not os.path.exists(directory):
        print(f"A '{directory}' könyvtár nem létezik!")
        return []
    
    results = []
    txt_files = [f for f in os.listdir(directory) if f.endswith('.txt')]
    
    print(f"{len(txt_files)} txt fájl található a '{directory}' könyvtárban")
    
    for filename in txt_files:
        filepath = os.path.join(directory, filename)
        
        try:
            with open(filepath, 'r', encoding='utf-8') as file:
                content = file.read()
                
                # Mutációk keresése
                deltaF508_count = len(re.findall(r'deltaF508|F508del|ΔF508|p\.Phe508del', content, re.IGNORECASE))
                cftr_dele23_count = len(re.findall(r'CFTRdele2,3|CFTR\s*dele2,3|dele2,3|del2,3', content, re.IGNORECASE))
                
                results.append({
                    'filename': filename,
                    'content': content,
                    'deltaF508_count': deltaF508_count,
                    'cftr_dele23_count': cftr_dele23_count,
                    'has_both': deltaF508_count > 0 and cftr_dele23_count > 0
                })
                
                print(f"Feldolgozva: {filename}")
                
        except Exception as e:
            print(f"Hiba a {filename} fájl olvasásakor: {e}")
    
    return results
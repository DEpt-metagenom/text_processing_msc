import os
import ollama # pip install ollama
from tqdm import tqdm

# --- ÚTVONALAK ---
# Ahol a Qwen-nel készült 47 db kis txt fájl van
IN_DIR = r"C:\Users\szk19\OneDrive\Asztali gép\Bioinformatika I félév\Új mappa\Szakirodalmi szövegbányászat\myplace\qwen2.5_32b"
# Ahol a végső, nagy hivatkozott jelentés lesz
FINAL_OUT = r"C:\Users\szk19\OneDrive\Asztali gép\Bioinformatika I félév\Új mappa\Szakirodalmi szövegbányászat\myplace\BRCA_Final_Synthesis_Qwen.md"

def main():
    # Fájlok listázása
    files = [f for f in os.listdir(IN_DIR) if f.endswith(".txt") and os.path.getsize(os.path.join(IN_DIR, f)) > 0]
    
    if not files:
        print(f"Hiba: Nem találok txt fájlokat itt: {IN_DIR}")
        return

    print(f"Indítás: {len(files)} fájl feldolgozása a Qwen2.5 modellel...")

    with open(FINAL_OUT, "w", encoding="utf-8") as f_out:
        f_out.write("# BRCA Tumor Szuppresszor Gének Szerepe az Emlőrákban\n")
        f_out.write("## Összefoglaló jelentés (Helyi Qwen szintézis)\n\n")

        for fname in tqdm(files):
            pmc_id = fname.replace("_summary.txt", "").replace(".txt", "")
            
            with open(os.path.join(IN_DIR, fname), "r", encoding="utf-8") as f_in:
                content = f_in.read()

            # Megkérjük a helyi Qwen-t, hogy formázza meg szakmaira
            prompt = (
                f"Te egy onkológiai bioinformatikus vagy. Az alábbi összefoglaló alapján írj egy "
                f"szakszerű magyar nyelvű bekezdést, amely bemutatja a kutatás lényegét. "
                f"A bekezdés végén kötelezően hivatkozz így: ({pmc_id}).\n\n"
                f"FORRÁS ADATOK:\n{content}"
            )

            try:
                # Az Ollama-t hívjuk meg (qwen2.5:7b-t javaslok a sebesség miatt)
                response = ollama.generate(model="qwen2.5:7b", prompt=prompt)
                
                f_out.write(response['response'] + "\n\n")
                f_out.flush() # Azonnal írja ki a lemezre

            except Exception as e:
                print(f"Hiba a {pmc_id} feldolgozásakor: {e}")

    print(f"\n✅ SIKER! A végleges szakmai jelentés elkészült: {FINAL_OUT}")

if __name__ == "__main__":
    main()

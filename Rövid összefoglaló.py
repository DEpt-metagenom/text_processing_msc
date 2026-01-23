import ollama
import os

INPUT_FILE = r"C:\Users\szk19\OneDrive\Asztali gép\Bioinformatika I félév\Új mappa\Szakirodalmi szövegbányászat\myplace\BRCA_Final_Synthesis_Qwen.md"
SHORT_REF_OUT = r"C:\Users\szk19\OneDrive\Asztali gép\Bioinformatika I félév\Új mappa\Szakirodalmi szövegbányászat\myplace\BRCA_Szakmai_Osszefoglalo_Sok_Hivatkozassal.md"

def generate_dense_summary():
    if not os.path.exists(INPUT_FILE):
        print("Hiba: A forrásfájl nem található!")
        return

    with open(INPUT_FILE, "r", encoding="utf-8") as f:
        data = f.read()

    # Szigorúbb, kényszerítő erejű prompt
    prompt = (
        "Feladat: Készíts egy professzionális tudományos összefoglalót a BRCA1/2 génekről.\n\n"
        "FORRÁSSZÖVEG (ebben keresd a PMC kódokat):\n"
        f"{data[:25000]}\n\n"
        "SZIGORÚ UTASÍTÁSOK:\n"
        "1. Írj magyarul egy sűrű, szakmai szöveget.\n"
        "2. MINDEN EGYES MONDAT végére tegyél legalább egy PMC hivatkozást a forrásszövegből (pl. PMC10055970).\n"
        "3. Ha egy állítás több forrásban is szerepel, sorold fel mindet (pl. PMC12345, PMC67890).\n"
        "4. Téma: BRCA szerepe a DNS-javításban (Homológ Rekombináció), mutációk hatása és PARP gátlás.\n"
        "5. Tilos hivatkozás nélküli állítást írni!"
    )

    print("Sűrű hivatkozásos szintézis generálása (Kényszerített hivatkozásokkal)...")
    try:
        response = ollama.generate(model="qwen2.5:7b", prompt=prompt)
        
        with open(SHORT_REF_OUT, "w", encoding="utf-8") as f_out:
            f_out.write("# BRCA1/2 Szerepe: Sűrített Tudományos Szintézis (Hivatkozásokkal)\n\n")
            f_out.write(response['response'])
        
        print(f"✅ KÉSZ! Ellenőrizd a fájlt: {SHORT_REF_OUT}")
    except Exception as e:
        print(f"Hiba történt: {e}")

if __name__ == "__main__":
    generate_dense_summary()

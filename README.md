# Szakirodalmi szövegbányászat gyakorlat Házifeladat Határidő 2026 Dec. 11

Önállóan, saját gépen gyakorolható feladatsor: parancssori szövegfeldolgozás (bash, awk, perl), Python szövegtisztítás, egy kis bioinformatikai mini projekt és egy házi feladat helyi nyelvi modellel.

## Hogyan dolgozzatok
1. **Fork** a repóról a saját GitHub fiókotokba, majd klónozás (részletek: [tasks/00_setup.md](tasks/00_setup.md)).
2. Dolgozzatok saját gépen (Linux/WSL), GitHub Codespacesben (`.devcontainer` benne van), vagy hasonló környezetben.
3. A feladatokat sorban oldjátok meg, a bemeneti fájlok a `samples_files/` mappában vannak (ne írjatok bele, a kimenetek az `out/` mappába menjenek).
4. Ellenőrzés: a `solutions/` mappában a várt kimeneteket (`*.expected`) találjátok. Először mindig saját megoldást készítsetek.

## Feladatok
| # | Téma | Fájl |
|---|---|---|
| 00 | Előkészületek (fork, környezet) | [tasks/00_setup.md](tasks/00_setup.md) |
| 01 | bash: grep, wc, sort, uniq, pipe | [tasks/01_bash.md](tasks/01_bash.md) |
| 02 | awk: számolás, naplók és táblázatok ellenőrzése | [tasks/02_awk.md](tasks/02_awk.md) |
| 03 | Reguláris kifejezések: grep -P, perl, sed | [tasks/03_regex_perl.md](tasks/03_regex_perl.md) |
| 04 | Python: beolvasás, tisztítás, tokenizálás | [tasks/04_python_text.md](tasks/04_python_text.md) |
| 05 | Mini projekt: FASTA → PHYLIP → fehérje FASTA | [tasks/05_fasta_mini_project.md](tasks/05_fasta_mini_project.md) |
| 06 | Házi feladat: saját projekt Ollama modellel | [tasks/06_homework_llm.md](tasks/06_homework_llm.md) |

## A repó tartalma
- `tasks/` – a feladatleírások
- `samples_files/` – bemeneti adatok: kis szövegfájlok, Project Gutenberg könyvek, NCBI baktérium-metaadatok, Nextflow naplók (egy szintetikus és egy valódi, anonimizált), phiX174 genom és génjei
- `solutions/` – várt kimenetek és a mini projekt mintamegoldása
- `llm_example/` – példa: weboldal feldolgozása és összefoglalás helyi LLM-mel
- `tools/` – oktatói segédszkript a szintetikus Nextflow napló újragenerálásához
- `szöveges dokumentumok 1.docx` – előadásanyag (fájltípusok, karakterkódolás)

## Számonkérés
Minden résztvevőnek nyitunk egy-egy issue-t. Lehetőleg ezt használjuk az óra teljesítéséhez kapcsolódó feladatok megbeszélésére: itt írhatjátok meg, hogy miből szeretnétek a házi feladatot (06) elkészíteni, és itt egyeztetünk a részletekről. A kész házi feladatot pull requestként a saját mappátokba (`vezeteknev_keresztnev/`) töltsétek fel.

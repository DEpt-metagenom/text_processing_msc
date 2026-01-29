import re

# File path
file_path = r"D:\Egyetem\bioinfo\szovegbanyaszat\pubmed_abstracts_10000.txt"
output_file = r"D:\Egyetem\bioinfo\szovegbanyaszat\filtered_abstracts_2.txt"

# Keyword lists
plasmid_keywords = {"plasmid", "plasmids", "hybrid plasmid", "mosaic plasmid" "kimeric plasmid"}

virulence_keywords = {
    "toxin", "virulence", "fimbriae", "adhesion", "biofilm", "capsule",
    "hemolysin", "siderophore", "invasion", "iroN", "sitABCD", "iucA", "hypervirulence"
}

resistance_keywords = {
    "antibiotic resistance", "resistance gene", "beta-lactamase", "bla",
    "carbapenemase", "aminoglycoside", "efflux pump", "multidrug resistance", "ESBL", "CTX-M"
}

# Function to check if text contains any of the keywords
def contains_keywords(text, keywords):
    text = text.lower()  # Case insensitive search
    return any(keyword in text for keyword in keywords)

# Read file
with open(file_path, "r", encoding="utf-8") as file:
    text = file.read()

# Split abstracts using regex (each abstract starts with a number followed by a dot)
abstracts = re.split(r"\n\d+\.", text)

# Filter abstracts
filtered_abstracts = []
for abstract in abstracts:
    if (
        contains_keywords(abstract, plasmid_keywords)
        and contains_keywords(abstract, virulence_keywords)
        and contains_keywords(abstract, resistance_keywords)
    ):
        filtered_abstracts.append(abstract.strip())

# Save filtered abstracts
with open(output_file, "w", encoding="utf-8") as file:
    file.write("\n\n".join(filtered_abstracts))

print(f"Filtered {len(filtered_abstracts)} abstracts containing all three categories.")

import re

# Read the Markdown file
with open("full_texts.md", "r", encoding="utf-8") as md_file:
    md_content = md_file.read()

# Define a regex pattern to extract all occurrences of "Materials and Methods" until "Results"
pattern = r"(?i)(materials and methods)(.*?)(?=\bresults\b|$)"  # Non-greedy match until "Results"

matches = re.findall(pattern, md_content, re.DOTALL)

# Extract all matched sections
if matches:
    extracted_texts = "\n\n".join(match[1].strip() for match in matches)  # Collect all sections

    # Save extracted content
    with open("materials_method.md", "w", encoding="utf-8") as output_file:
        output_file.write("Materials and Methods\n\n" + extracted_texts)

    print(f"Extracted {len(matches)} 'Materials and Methods' sections and saved as materials_method.md")
else:
    print("No 'Materials and Methods' section found.")

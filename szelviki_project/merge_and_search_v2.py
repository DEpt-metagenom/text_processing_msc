import os
from flashtext import KeywordProcessor

# Define folder paths
folder_path = "C:/Users/schmi/Desktop/extracted_materials_methods"  # Update with your actual path
merged_file = "merged_text.txt"
highlighted_output = "highlighted_results.txt"

# Keywords to search and highlight
keywords = ["ctab", "mechanical"]
keyword_processor = KeywordProcessor()
for kw in keywords:
    keyword_processor.add_keyword(kw)

# Step 1: Merge all .txt files into one
with open(merged_file, "w", encoding="utf-8", errors="ignore") as merged_f:
    for filename in os.listdir(folder_path):
        if filename.endswith(".txt"):
            file_path = os.path.join(folder_path, filename)
            with open(file_path, "r", encoding="utf-8", errors="ignore") as f:
                merged_f.write(f.read() + "\n\n")  # Add space between files

print(f"✅ Merged all .txt files into {merged_file}")

# Step 2: Extract keyword snippets with context
with open(merged_file, "r", encoding="utf-8", errors="ignore") as f:
    text = f.read()

matches = keyword_processor.extract_keywords(text, span_info=True)
highlighted_snippets = []

for match, start, end in matches:
    snippet_start = max(0, start - 100)  # 100 chars before
    snippet_end = min(len(text), end + 100)  # 100 chars after
    snippet = text[snippet_start:snippet_end]

    # Highlight the keyword
    highlighted_snippet = snippet.replace(match, f"**{match.upper()}**")
    highlighted_snippets.append(highlighted_snippet)

# Step 3: Save highlighted snippets to a file
if highlighted_snippets:
    with open(highlighted_output, "w", encoding="utf-8", errors="ignore") as f:
        f.write("\n\n---\n\n".join(highlighted_snippets))

    print(f"✅ Highlighted snippets saved to {highlighted_output}")
else:
    print("❌ No matches found.")

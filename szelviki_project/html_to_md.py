import html2text

# Read the HTML file
with open("full_texts.html", "r", encoding="utf-8") as html_file:
    html_content = html_file.read()

# Convert HTML to Markdown
markdown_text = html2text.html2text(html_content)

# Save the Markdown output
with open("full_texts.md", "w", encoding="utf-8") as md_file:
    md_file.write(markdown_text)

print("Conversion complete! Markdown saved as full_texts.md")

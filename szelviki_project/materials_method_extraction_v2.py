import os
from bs4 import BeautifulSoup

# Define folder paths
folder_path = r"C:/Users/schmi/Desktop/data_mining"  # Replace with your actual folder path
output_folder = "extracted_materials_methods"

# Create output folder if it doesn't exist
os.makedirs(output_folder, exist_ok=True)

def extract_materials_methods(html_content):
    """Extracts text from <p> tags in sections with sec-type='materials|methods' 
       or in sections where <title> contains 'method' or 'materials'.
    """
    soup = BeautifulSoup(html_content, "html.parser")
    extracted_texts = []

    # Extract from sections with sec-type="materials" or "methods"
    for section in soup.find_all("sec", {"sec-type": ["materials", "methods"]}):
        paragraphs = section.find_all("p")
        extracted_texts.append("\n".join(p.get_text() for p in paragraphs))
    
    # Extract from sections where <title> contains "method" or "materials"
    for title in soup.find_all("title"):
        if title.get_text().lower() in ["method", "methods", "material", "materials"]:
            section = title.find_parent("sec")  # Find the parent section of the title
            if section:
                paragraphs = section.find_all("p")
                extracted_texts.append("\n".join(p.get_text() for p in paragraphs))

    return "\n\n".join(extracted_texts) if extracted_texts else None

# Process all HTML files in the folder
for filename in os.listdir(folder_path):
    if filename.endswith(".html"):
        file_path = os.path.join(folder_path, filename)
        
        # Read file with UTF-8 and ignore errors
        with open(file_path, "r", encoding="utf-8", errors="ignore") as f:
            html_content = f.read()
        
        # Extract the text
        extracted_text = extract_materials_methods(html_content)
        
        if extracted_text:
            # Save extracted text to a new file
            output_file = os.path.join(output_folder, f"{filename}.txt")
            with open(output_file, "w", encoding="utf-8", errors="ignore") as out_f:
                out_f.write(extracted_text)
            
            print(f"✅ Extracted materials/methods section from {filename} → Saved to {output_file}")
        else:
            print(f"❌ No materials/methods section found in {filename}.")

import torch
import google.generativeai as genai
from google.colab import userdata
import os
import re
import requests
from bs4 import BeautifulSoup
from datetime import date
import time  # Import the time module for adding delays
from transformers import AutoTokenizer, AutoModelForSeq2SeqLM # Changed imports to AutoTokenizer, AutoModelForSeq2SeqLM

# --- GPU Check and Device Configuration ---
device = "cuda" if torch.cuda.is_available() else "cpu" # Determine device: GPU if available, else CPU
# print(f"Using device: {device}") # DEBUG: Print the device being used - Commented out

if device == "cuda":
    # gpu_info = !nvidia-smi # DEBUG: GPU Info - Commented out
    # gpu_info_str = '\n'.join(gpu_info) # DEBUG: GPU Info - Commented out
    # print(gpu_info_str) # DEBUG: GPU Info - Commented out
    pass # Keep GPU info retrieval for potential debugging, but don't print in default output
else:
    print("No GPU detected. Running on CPU.")

# Load tokenizer and model outside the function (for efficiency - load only once)
tokenizer_hunsum = AutoTokenizer.from_pretrained("SZTAKI-HLT/mT5-base-HunSum-1") # Using mT5-base-HunSum-1 now
model_hunsum = AutoModelForSeq2SeqLM.from_pretrained("SZTAKI-HLT/mT5-base-HunSum-1").to(device) # Using mT5-base-HunSum-1 now and move to device


def clean_text(text):
    """
    Remove known ad-like or repetitive phrases from the text.
    """
    ad_phrases = [
        "5 könyvTöbb mint 600 meghökkentő, érdekes és tanulságos történet!",
        "Kövesse az Indexet Facebookon is!"
    ]
    for phrase in ad_phrases:
        text = text.replace(phrase, "")
    return text.strip()

def get_search_results(query):
    """
    Constructs the search URL using dynamic dates and retrieves the HTML of the search results.
    """
    base_url = "https://index.hu/24ora/"
    today = date.today().isoformat()  # Today's date in YYYY-MM-DD format.
    params = {
        "word": "1",
        "pepe": "1",
        "tol": "1999-01-01",  # Fixed start date.
        "ig": today,          # End date set to today.
        "s": query
    }
    headers = {"User-Agent": "Mozilla/5.0"}
    response = requests.get(base_url, params=params, headers=headers)
    response.raise_for_status()
    return response.text

def parse_results(html):
    """
    Parses the HTML to extract up to the first 5 articles.
    For each article, extracts:
      - Date (from the <a> tag with class "datum")
      - Title and Link (from the <a> tag with a title attribute)
      - Short Description (from the <p> tag with class "ajanlo")
    """
    soup = BeautifulSoup(html, "html.parser")
    main_section = soup.find("main", id="talalatok")
    if not main_section:
        print("Main container with id 'talalatok' not found.")
        return []

    articles = main_section.find_all("article", class_="rovatajanlo", limit=5)
    results = []
    for article in articles:
        # Extract Date.
        date_tag = article.find("a", class_="datum")
        if date_tag:
            article_date = date_tag.get("data-date", date_tag.get_text(strip=True))
        else:
            article_date = "No Date Found"

        # Extract Title and Link.
        title_tag = article.find("a", title=True)
        if title_tag:
            title = title_tag.get("title")
            link = title_tag.get("href")
        else:
            title = "No Title Found"
            link = "No Link Found"

        # Extract and clean the Short Description.
        desc_tag = article.find("p", class_="ajanlo")
        description = desc_tag.get_text(strip=True) if desc_tag and desc_tag.get_text(strip=True) else "No Description Found"
        description = clean_text(description)

        results.append({
            "title": title,
            "link": link,
            "description": description,
            "date": article_date
        })
    return results

def get_article_text(link):
    """
    Follows the article link and attempts to extract the full article text.
    Tries several containers to account for variations in page layout.
    Returns the extracted text, but does NOT print it.
    """
    headers = {"User-Agent": "Mozilla/5.0"}
    try:
        response = requests.get(link, headers=headers)
        response.raise_for_status()
    except Exception as e:
        print(f"Failed to fetch article at {link}: {e}")
        return ""

    soup = BeautifulSoup(response.text, "html.parser")

    # Try different containers for the article text.
    article_body = soup.find("div", class_="cikk-torzs")
    if not article_body:
        article_body = soup.find("div", class_="cikk-tartalom")
    if not article_body:
        article_body = soup.select_one("div.text.clearafter")

    if article_body:
        # Extract text from all <p> tags.
        paragraphs = article_body.find_all("p")
        text = "\n".join(p.get_text(strip=True) for p in paragraphs)
        return clean_text(text)
    return ""

def generate_description_gemini(article_text):
    """
    Generates a concise Hungarian description (max 2 sentences) using the Gemini 1.5 Flash API.
    """
    try:
        model = genai.GenerativeModel('gemini-1.5-flash')
        prompt = f"""Magyarul foglalja össze maximum 2 mondatban a következő cikket:

        {article_text}
        """

        response = model.generate_content(prompt)

        if response.prompt_feedback and response.prompt_feedback.block_reason:
            return f"Gemini API blocked response: {response.prompt_feedback.block_reason_message}"
        if not response.text:
            return "Gemini API returned empty response."

        return response.text.strip()

    except Exception as e:
        return f"Error generating Gemini summary: {e} - {type(e).__name__}"


def generate_description_hunsum(article_text): # Function name updated to generate_description_hunsum
    """
    Generates a concise Hungarian description using locally loaded SZTAKI Bert2Bert-HunSum-1 model.
    Addresses RuntimeError by limiting input length.
    Experimenting with summarization parameters (Option 1) - Decreased length_penalty to reduce repetition.
    Explicitly using GPU if available.
    Refined Prompt (Hungarian - more directive for conciseness and key aspects).
    Removed 'prefix' argument from model.generate() as it caused ValueError.
    """
    try:  # <--- TRY BLOCK ADDED HERE
        # Tokenize input text with explicit max_length and truncation
        input_ids = tokenizer_hunsum.encode(
            article_text,
            max_length=512,  # Explicitly set max_length for tokenization
            truncation=True,
            return_tensors="pt",
        ).to(device) # Move input_ids to GPU if available

        # Refined Hungarian Prompt - More directive for conciseness and key aspects
        prompt = f"""Foglald össze magyarul a következő cikket maximum 3 mondatban, tömören és lényegre törően, kiemelve a legfontosabb tényeket és következtetéseket. Kerüld a felesleges ismétléseket: {article_text}"""

        # Generate summary with adjusted parameters (length_penalty decreased AND reverted max_length/min_length)
        summary_ids = model_hunsum.generate(
            input_ids,
            max_length=150,  # Reverted max_length back to 150 (from 350) - Shorter summaries
            min_length=30,  # Reverted min_length back to 30 (from 70) - Shorter summaries
            length_penalty=2.0, # length_penalty now at 2.0
            num_beams=5,      # num_beams kept at 5
            early_stopping=True,
        )

        # Decode summary
        summary_text = tokenizer_hunsum.decode(summary_ids[0], skip_special_tokens=True)
        return summary_text.strip()

    except Exception as e: # <--- EXCEPT BLOCK ADDED HERE - THIS WAS MISSING!
        return f"Error generating HunSum-1 summary locally: {e} - {type(e).__name__}"

def analyze_sentiment_gemini(article_text):
    """
    Analyzes sentiment, VERY ROBUST parsing, forces sentiment for ALL topics, handles formatting issues.
    """
    model = genai.GenerativeModel('gemini-pro')

    prompt = f"""Elemezze a következő cikkszöveget a tárgyalt fő témák azonosítása és az egyes témákhoz kapcsolódó érzelmek kifejezése céljából.
**FELTÉTLENÜL ÉRTÉKELJEN MINDEN TÉMÁT ÉRZELMILEG!** Ne hagyjon ki egyetlen témát sem az értékelésből.
Értékelje az egyes témákhoz tartozó érzelmeket egy 5 fokozatú skálán magyarul: Nagyon Negatív, Negatív, Vegyes, Pozitív, Nagyon Pozitív. **MINDEN TÉMÁHOZ adjon meg egy rövid, zárójeles magyarázatot is az érzésről, MÉG AKKOR IS, HA AZ ÉRZÉS NEHEZEN MEGFOGHATÓ, SEMLEGES VAGY ÖSSZETETT. HA NEM BIZTOS AZ ÉRZELMI ÉRTÉKELÉSBEN, AKKOR IS PRÓBÁLJON MEGADNI EGYET!**

Cikkszöveg:{article_text}

Az eredményt témák listájaként adja meg az érzelmi értékelésükkel, **AZ ÉRZÉST ZÁRÓJELBEN FELTÜNTETVE. HA EGY TÉMÁHOZ NEM LEHET EGYÉRTELMŰ ÉRZELMET RENDELNI, AKKOR HASZNÁLJA A "SEMLEGES" ÉRTÉKELÉST, ÉS MAGYARÁZZA MEG RÖVIDEN, MIÉRT SEMLEGES.  FONTOS, HOGY MINDEN TÉMA ÉRTÉKELVE LEGYEN!** Például:

Témák és Érzelmek:
- 1. Téma: [Érzelmi Értékelés] **(Fő érzés vagy magyarázat)**
- 2. Téma: [Érzelmi Értékelés] **(Fő érzés vagy magyarázat)**
- 3. Téma: [Érzelmi Értékelés] **(Fő érzés vagy magyarázat)**
...
"""

    try:
        response = model.generate_content(prompt)
        response.resolve()

        if response.text:
            sentiment_results = {}
            lines = response.text.strip().split('\n')
            topic_number = 1 # Initialize topic number counter
            for line in lines:
                if line.startswith('- ') or re.match(r'^\d+\.\s*Téma:\s*\*\*', line): # Handle both '-' and numbered formats, and "Téma: **"
                    line_text = line.lstrip('- ').lstrip() # Remove leading '-', and spaces
                    parts = line_text.split(':') # Split at the first ":"
                    if len(parts) >= 2: # Ensure there are at least topic and sentiment parts
                        topic_part = parts[0].replace('Téma', '').replace('*', '').replace('.', '').replace('*', '').strip() # EVEN MORE aggressive cleaning, remove extra '*'
                        sentiment_part = ":".join(parts[1:]).strip() # Reconstruct sentiment part if split occurred within it
                        sentiment_rating = sentiment_part.split('(')[0].strip().replace('[', '').replace(']', '').replace('*', '').strip().replace('SEMLEGES', 'Semleges') # Extract rating, remove extra '*' too, normalize "SEMLEGES"
                        feeling = ""
                        if '(' in sentiment_part and ')' in sentiment_part:
                            feeling = sentiment_part.split('(')[1].split(')')[0].strip() # Extract feeling

                        topic_key = f"{topic_number}. Téma: {topic_part.strip()}" if re.match(r'^\d+\.\s*Téma:\s*\*\*', line) else topic_part.strip() # Keep numbered format if detected, REMOVE BOLDING FROM TOPIC KEY

                        sentiment_results[topic_key] = f"{sentiment_rating} ({feeling})" if feeling else sentiment_rating
                        topic_number += 1 # Increment topic number

                    elif len(parts) == 1: # Handle lines with only topic, no rating or colon
                        topic_part = parts[0].strip()
                        sentiment_results[topic_part] = "Semleges (Nincs értékelés)" # Assign neutral sentiment

            return sentiment_results
        else:
            print("Gemini API returned empty response text for sentiment analysis.")
            return None

    except Exception as e:
        print(f"Error during Gemini sentiment analysis: {e}")
        return None


def main():
    # --- API Key Configuration ---
    # Gemini API Key
    try:
        GOOGLE_API_KEY = userdata.get('GOOGLE_API_KEY')
    except KeyError:
        GOOGLE_API_KEY = userdata.get('GOOGLE_API_KEY')

    if not GOOGLE_API_KEY:
        print("Error: Gemini API key not found. Please set it as a Colab secret or environment variable 'GOOGLE_API_KEY'.")
        return
    genai.configure(api_key=GOOGLE_API_KEY)
    #print("Gemini API Key configured successfully.")

    # Hugging Face API Token - Not needed for local model
    HUGGINGFACE_API_TOKEN = None # Set to None as it's not used

    use_hunsum = True # Enable HunSum-1 (local)
    use_sentiment = True # Enable Sentiment Analysis

    query = input("Enter your search keyphrase: ").strip()

    try:
        html = get_search_results(query)
    except Exception as e:
        print("An error occurred while fetching the search results:", e)
        return

    articles = parse_results(html)
    if not articles:
        print("No articles found.")
        return


    for i, article in enumerate(articles, start=1):
        print("\n" + "=" * 40)
        # Display Title and scraped short description.
        print(f"Title: {article['title']}")
        print(f"Scraped Short Description: {article['description']}")

        # Initialize sentiment_analysis_result here
        sentiment_analysis_result = None

        # Fetch the full article text to generate a new description.
        article_text = get_article_text(article["link"])
        if not article_text:
            gemini_summary = "Failed to fetch full article text."
            hunsum_summary = "Failed to fetch full article text." # Variable name updated to hunsum_summary
            # sentiment_analysis_result = None <- No need to initialize again here
        else:
            time.sleep(1) # Delay before Gemini call
            gemini_summary = generate_description_gemini(article_text)

            hunsum_summary = "HunSum-1 summarization (local)." # Message updated to local HunSum-1
            if use_hunsum: # Variable name updated to use_hunsum
                time.sleep(1) # Delay before HunSum-1 call (optional, but good practice)
                hunsum_summary = generate_description_hunsum(article_text) # Function name updated to generate_description_hunsum

            if use_sentiment: # Only call sentiment analysis if enabled and article text is available
                time.sleep(1) # Delay before sentiment analysis call
                sentiment_analysis_result = analyze_sentiment_gemini(article_text)


        # Display the generated summaries, URL, and date.
        print("Generated Gemini Description:")
        if isinstance(gemini_summary, str): # Check if it's a string before splitting
            for sentence in gemini_summary.split('. '): # Split into sentences and print each
                sentence = sentence.strip()
                if not sentence.endswith('.'): # Check if sentence ends with a period
                    print(f"  {sentence}.") # Add period if missing
                else:
                    print(f"  {sentence}") # Print as is if period exists
        else:
            print(f"  {gemini_summary}") # Print error message if not a string


        if use_hunsum: # Variable name updated to use_hunsum
            print("Generated HunSum-1 Description (Local):")
            if isinstance(hunsum_summary, str): # Check if it's a string before splitting
                for sentence in hunsum_summary.split('. '): # Split into sentences and print each
                    sentence = sentence.strip()
                    if not sentence.endswith('.'): # Check if sentence ends with a period
                        print(f"  {sentence}.") # Add period if missing
                    else:
                        print(f"  {sentence}") # Print as is if period exists
        else:
             print(f"Generated HunSum-1 Description: {hunsum_summary}") # Output label for HunSum-1 disabled state

        if use_sentiment and sentiment_analysis_result: # Now sentiment_analysis_result is guaranteed to be initialized
            print("\nSentiment Analysis (Gemini):")
            if sentiment_analysis_result:
                for topic, sentiment in sentiment_analysis_result.items():
                    print(f"  {topic}: {sentiment}") # Removed '-' prefix and ** marks
            else:
                print("Sentiment analysis failed.")


        print(f"URL: {article['link']}")
        print(f"Date: {article['date']}")
        print("=" * 40)

if __name__ == "__main__":
    main()
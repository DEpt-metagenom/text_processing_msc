Az első scripttel (transcript_to_json.py) a docx formátumú nyers átiratot json formátummá alakítottam. Először ezt adtam át az LLM-nek, de így nem működött,
a dialogue/predialogue elkülönítések összekuszálták, így egy újabb scripttel (transcript_cleaner.py) a json file-ból eldobattam a terapeuta és az avatar hozzászólásokat,
a dialogue/predialogue elkülönítésekkel együtt. Ezt a módosított json file-t töltöttem fel a gpt-4.1-mininek (emotion_labeller.py), érzelmi töltet szerinti címkézésre.

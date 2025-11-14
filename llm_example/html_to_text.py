from bs4 import BeautifulSoup
import ollama

with open("C:\\Users\\lg_gr\\Downloads\\article.htm", encoding="utf-8") as file:
    content = file.read()

soup = BeautifulSoup(content, "html.parser")
paras = soup.find_all("p")
for i, para in enumerate(paras[1:89]):
    print(f"Paragraph {i}: {para.get_text()}")

response = ollama.generate(model="gpt-oss:20b",
                           prompt="Summarise the research approach proposed in this article, its advantages and limitations:\n" + "\n".join(para.get_text() for para in paras[1:89])
)
print(response.text)

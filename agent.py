import requests
import time
import os
import html
from Bio import Entrez
import google.generativeai as genai
from google.generativeai.types import HarmCategory, HarmBlockThreshold

# --- НАСТРОЙКИ ---
Entrez.email = "tvoj_email@example.com" 

TELEGRAM_TOKEN = os.environ.get("TG_TOKEN")
TELEGRAM_CHAT_ID = os.environ.get("TG_CHAT_ID")
GEMINI_API_KEY = os.environ.get("GEMINI_API_KEY")

HISTORY_FILE = "history.txt"

# --- СМЯГЧЕННЫЙ ФИЛЬТР ---
# Добавили просто "Review", чтобы находить больше интересных статей
QUALITY_FILTER = " AND (Meta-Analysis[ptyp] OR Randomized Controlled Trial[ptyp] OR Systematic Review[ptyp] OR Review[ptyp])"

if GEMINI_API_KEY:
    genai.configure(api_key=GEMINI_API_KEY)

# --- ТВОЯ БАЗА ЗНАНИЙ ---
RAW_QUERIES = {
    "🧬 1. Фундамент": [
        "(mitochondrial biogenesis) AND (exercise)",
        "(NAD+) AND (aging)",
        "(AMPK) AND (fasting)",
        "(metabolic flexibility) AND (fat oxidation)",
        "(hormesis) AND (sauna)"
    ],
    "💊 2. Нутрицевтика": [
        "(creatine) AND (brain)",
        "(magnesium) AND (sleep)",
        "(omega-3) AND (recovery)",
        "(vitamin D) AND (performance)",
        "(caffeine) AND (performance)",
        "(beta-alanine) AND (exercise)",
        "(nitrate) AND (blood flow)",
        "(Ashwagandha) AND (stress)",
        "(Rhodiola) AND (fatigue)",
        "(Lion's mane) AND (cognitive)",
        "(Cordyceps) AND (endurance)",
        "(Tongkat Ali) AND (testosterone)"
    ],
    "💪 3. Сила и Механика": [
        "(hypertrophy) AND (volume)",
        "(eccentric training) AND (tendon)",
        "(blood flow restriction) AND (hypertrophy)",
        "(plyometric training) AND (sprint)",
        "(velocity-based training)"
    ],
    "🫁 4. Кардио": [
        "(VO2max) AND (longevity)",
        "(heart rate variability) AND (recovery)",
        "(respiratory muscle training)"
    ],
    "🧠 5. Мозг": [
        "(mental toughness) AND (sport)",
        "(flow state) AND (performance)",
        "(neuroplasticity) AND (exercise)",
        "(tDCS) AND (sport)"
    ],
    "💤 6. Сон": [
        "(circadian rhythm) AND (performance)",
        "(sleep deprivation) AND (recovery)",
        "(deep sleep) AND (recovery)"
    ]
}

def load_history():
    if not os.path.exists(HISTORY_FILE): return set()
    with open(HISTORY_FILE, "r") as f: return set(line.strip() for line in f)

def save_history(new_ids):
    with open(HISTORY_FILE, "a") as f:
        for pmid in new_ids: f.write(f"{pmid}\n")

def analyze_abstract_with_gemini(title, abstract):
    if not GEMINI_API_KEY: return "⚠️ Нет ключа Gemini"

    prompt = f"""
    Analyze this abstract for a biohacker.
    Title: {title}
    Abstract: {abstract}
    
    Task:
    1. Summarize finding in ONE sentence in RUSSIAN.
    2. Format: "✅ [Supplement/Method] -> [Result] (value/%)."
    """
    
    safety_settings = {
        HarmCategory.HARM_CATEGORY_HARASSMENT: HarmBlockThreshold.BLOCK_NONE,
        HarmCategory.HARM_CATEGORY_HATE_SPEECH: HarmBlockThreshold.BLOCK_NONE,
        HarmCategory.HARM_CATEGORY_SEXUALLY_EXPLICIT: HarmBlockThreshold.BLOCK_NONE,
        HarmCategory.HARM_CATEGORY_DANGEROUS_CONTENT: HarmBlockThreshold.BLOCK_NONE,
    }

    try:
        model = genai.GenerativeModel('gemini-1.5-flash')
        response = model.generate_content(prompt, safety_settings=safety_settings)
        return response.text.strip() if response.text else "⚠️ Пустой ответ AI"
    except Exception as e:
        return f"Ошибка AI: {str(e)}"

def search_pubmed(query, days=None, retmax=5, sort="date"):
    full_query = query + QUALITY_FILTER
    try:
        params = {"db": "pubmed", "term": full_query, "retmax": retmax, "sort": sort}
        if days:
            params["reldate"] = days
            params["datetype"] = "pdat"
        handle = Entrez.esearch(**params)
        record = Entrez.read(handle)
        handle.close()
        return record["IdList"]
    except Exception as e:
        print(f"❌ Ошибка поиска '{query}': {e}")
        return []

def fetch_details_and_analyze(id_list):
    if not id_list: return []
    ids = ",".join(id_list)
    try:
        handle = Entrez.efetch(db="pubmed", id=ids, retmode="xml")
        records = Entrez.read(handle)
        handle.close()
        papers = []
        for article in records['PubmedArticle']:
            try:
                pmid = article['MedlineCitation']['PMID']
                title = article['MedlineCitation']['Article']['ArticleTitle']
                abstract_list = article['MedlineCitation']['Article'].get('Abstract', {}).get('AbstractText', [])
                abstract = " ".join(abstract_list) if abstract_list else ""
                
                if not abstract:
                    summary = "Нет абстракта"
                else:
                    summary = analyze_abstract_with_gemini(title, abstract)
                
                link = f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/"
                pub_date = article['MedlineCitation']['Article']['Journal']['JournalIssue']['PubDate']
                year = pub_date.get('Year', 'N/A')
                
                papers.append({'summary': summary, 'link': link, 'id': str(pmid), 'year': year})
            except: continue
        return papers
    except: return []

def send_telegram_message(message):
    if not TELEGRAM_TOKEN or not TELEGRAM_CHAT_ID: return
    url = f"https://api.telegram.org/bot{TELEGRAM_TOKEN}/sendMessage"
    data = {"chat_id": TELEGRAM_CHAT_ID, "text": message, "parse_mode": "HTML", "disable_web_page_preview": True}
    requests.post(url, data=data)

def main():
    print("📢 Запуск агента v4.1 (LOGGING MODE)...")
    seen_ids = load_history()
    all_papers = []
    new_seen_ids = []

    # --- ЭТАП 1: СВЕЖЕЕ ---
    print("\n🔎 Ищу свежие статьи (24ч)...")
    found_fresh = False
    for category, query_list in RAW_QUERIES.items():
        for q in query_list:
            ids = search_pubmed(q, days=1, retmax=1)
            unique_ids = [i for i in ids if i not in seen_ids]
            
            # ЛОГИРОВАНИЕ: Пишем, сколько нашли
            if len(unique_ids) > 0:
                print(f"   ✅ {q[:20]}... -> Найдено: {len(unique_ids)}")
                details = fetch_details_and_analyze(unique_ids)
                for paper in details:
                    paper['category'] = category
                    paper['type'] = 'fresh'
                    all_papers.append(paper)
                    seen_ids.add(paper['id'])
                    new_seen_ids.append(paper['id'])
                found_fresh = True
            else:
                # Если 0, молчим, чтобы не засорять лог
                pass
            time.sleep(0.5)

    # --- ЭТАП 2: АРХИВ (ЕСЛИ МАЛО) ---
    if len(all_papers) < 10:
        print(f"\n📚 Мало свежего ({len(all_papers)}). Ищу в архиве (10 лет!)...")
        needed = 15 - len(all_papers)
        
        for category, query_list in RAW_QUERIES.items():
            if needed <= 0: break
            for q in query_list:
                # Ищем за 3650 дней (10 лет) и берем по 2 самые релевантные
                ids = search_pubmed(q, days=3650, retmax=2, sort="relevance")
                unique_ids = [i for i in ids if i not in seen_ids]
                
                print(f"   🔎 Архив {q[:20]}... -> Новых: {len(unique_ids)}")
                
                if unique_ids:
                    # Берем только 1, чтобы разнообразить
                    to_take = unique_ids[:1]
                    details = fetch_details_and_analyze(to_take)
                    for paper in details:
                        paper['category'] = category
                        paper['type'] = 'archive'
                        all_papers.append(paper)
                        seen_ids.add(paper['id'])
                        new_seen_ids.append(paper['id'])
                        needed -= 1
                    time.sleep(1)

    if not all_papers:
        print("\n❌ ИТОГ: Ничего не найдено. Проверь запросы или фильтры.")
        # Отправим отладочное сообщение в ТГ, чтобы ты знал
        send_telegram_message("⚠️ Агент завершил работу, но не нашел ни одной статьи по фильтрам.")
        return

    # --- ОТПРАВКА ---
    print(f"\n📨 Подготовка к отправке {len(all_papers)} статей...")
    buffer_message = "<b>🧠 Biohack Digest (v4.1)</b>\n\n"
    current_category = ""
    
    # Сортировка по категориям
    all_papers.sort(key=lambda x: x.get('category', ''))

    for paper in all_papers:
        article_text = ""
        if paper.get('category') != current_category:
            article_text += f"<b>🔹 {paper.get('category')}</b>\n"
            current_category = paper.get('category')
        
        icon = "🔥" if paper['type'] == 'fresh' else "📜"
        clean_summary = html.escape(paper['summary']).replace("**", "")
        
        article_text += f"{icon} <a href='{paper['link']}'>Источник</a> ({paper['year']})\n{clean_summary}\n\n"
        
        if len(buffer_message) + len(article_text) > 3000:
            send_telegram_message(buffer_message)
            buffer_message = article_text
        else:
            buffer_message += article_text

    if buffer_message:
        send_telegram_message(buffer_message)

    if new_seen_ids:
        save_history(new_seen_ids)
        print("✅ Успех. История обновлена.")

if __name__ == "__main__":
    main()

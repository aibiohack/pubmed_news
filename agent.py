import requests
import time
import os
import html
from Bio import Entrez
import google.generativeai as genai

# --- НАСТРОЙКИ ---
Entrez.email = "tvoj_email@example.com" 

TELEGRAM_TOKEN = os.environ.get("TG_TOKEN")
TELEGRAM_CHAT_ID = os.environ.get("TG_CHAT_ID")
GEMINI_API_KEY = os.environ.get("GEMINI_API_KEY")

HISTORY_FILE = "history.txt"
QUALITY_FILTER = " AND (Meta-Analysis[ptyp] OR Randomized Controlled Trial[ptyp] OR Systematic Review[ptyp])"

# Настройка Gemini
if GEMINI_API_KEY:
    genai.configure(api_key=GEMINI_API_KEY)

RAW_QUERIES = {
    "Метаболизм и восстановление": [
        "(mitochondrial biogenesis) AND (exercise)",
        "(NAD+) AND (exercise performance) OR (skeletal muscle)",
        "(lactate metabolism) OR (lactate as signaling molecule) AND (training)",
        "(sleep extension) OR (sleep quality) AND (athletic performance)",
        "(cold exposure) OR (cryotherapy) AND (recovery) OR (inflammation)",
        "(heat acclimation) OR (sauna) AND (performance) OR (heat shock proteins)"
    ],
    "Нутрицевтика": [
        "(ketogenic diet) AND (endurance performance) NOT (epilepsy)",
        "(exogenous ketones) AND (endurance) OR (recovery)",
        "(nitrate supplementation) AND (beetroot juice) AND (exercise efficiency)",
        "(caffeine timing) OR (low-dose caffeine) AND (performance)",
        "(creatine supplementation) AND (cognitive function) OR (recovery)",
        "(beta-alanine) AND (high-intensity exercise)",
        "(collagen peptides) AND (tendon) OR (ligament)",
        "(\"time-restricted eating\") OR (\"intermittent fasting\") AND (body composition) OR (performance)"
    ],
    "Нейромышечный контроль": [
        "(post-activation potentiation) AND (PAP) protocols",
        "(blood flow restriction) OR (KAATSU training) AND (hypertrophy) OR (rehabilitation)",
        "(velocity-based training) AND (strength)",
        "(tendon stiffness) AND (performance) OR (injury prevention)",
        "(electromyostimulation) OR (EMS) AND (recovery) OR (performance)"
    ],
    "Мониторинг и персонализация": [
        "(wearable devices) AND (heart rate variability) HRV AND (recovery monitoring)",
        "(muscle oximetry) OR (NIRS) AND (training load)",
        "(genetic polymorphisms) AND (response to exercise) OR (injury risk)",
        "(microbiome) AND (athlete) OR (immune function)",
        "(\"omics\" in sports) (metabolomics, proteomics, transcriptomics)"
    ],
    "Ментальный биохакинг": [
        "(transcranial direct current stimulation) tDCS AND (motor learning) OR (endurance)",
        "(neurofeedback) AND (sports performance)",
        "(mindfulness) OR (meditation) AND (sport) OR (recovery)",
        "(vagus nerve stimulation) AND (recovery)"
    ]
}

def load_history():
    if not os.path.exists(HISTORY_FILE):
        return set()
    with open(HISTORY_FILE, "r") as f:
        return set(line.strip() for line in f)

def save_history(new_ids):
    with open(HISTORY_FILE, "a") as f:
        for pmid in new_ids:
            f.write(f"{pmid}\n")

# --- МОДУЛЬ АНАЛИЗА (С ЗАЩИТОЙ ОТ СБОЕВ) ---
def analyze_abstract_with_gemini(title, abstract):
    if not GEMINI_API_KEY:
        return "⚠️ Нет ключа Gemini"

    prompt = f"""
    Задача: Проанализируй научный абстракт как спортивный физиолог.
    
    Заголовок: {title}
    Текст: {abstract}

    Требования:
    1. Напиши ОДНО предложение на русском языке.
    2. Формат: "✅ [Суть вмешательства] на [Кол-во людей/животных] -> [Результат/Вывод] (цифры/проценты если есть)."
    3. Будь предельно краток.
    """

    # Список моделей: сначала пробуем новую, если нет - старую надежную
    models_to_try = ['gemini-1.5-flash', 'gemini-pro']

    for model_name in models_to_try:
        try:
            model = genai.GenerativeModel(model_name)
            response = model.generate_content(prompt)
            
            if response.text:
                return response.text.strip()
        except Exception:
            # Если не вышло с этой моделью, пробуем следующую
            continue
    
    # Если все модели не сработали
    return f"Заголовок: {title} (Не удалось проанализировать через AI)"

# --- ПОИСК ---
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
        print(f"Ошибка поиска: {e}")
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
                title_en = article['MedlineCitation']['Article']['ArticleTitle']
                
                abstract_parts = article['MedlineCitation']['Article'].get('Abstract', {}).get('AbstractText', [])
                full_abstract = " ".join(abstract_parts) if abstract_parts else ""

                if not full_abstract:
                    summary = f"Заголовок: {title_en} (Нет текста статьи)"
                else:
                    summary = analyze_abstract_with_gemini(title_en, full_abstract)
                
                link = f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/"
                pub_date = article['MedlineCitation']['Article']['Journal']['JournalIssue']['PubDate']
                year = pub_date.get('Year', 'N/A')
                
                papers.append({'summary': summary, 'link': link, 'id': str(pmid), 'year': year})
            except Exception as e:
                print(f"Ошибка статьи {pmid}: {e}")
                continue
        return papers
    except Exception as e:
        print(f"Ошибка скачивания: {e}")
        return []

# --- TELEGRAM ---
def send_telegram_message(message):
    if not TELEGRAM_TOKEN or not TELEGRAM_CHAT_ID:
        print("❌ Токены не настроены")
        return
    url = f"https://api.telegram.org/bot{TELEGRAM_TOKEN}/sendMessage"
    data = {"chat_id": TELEGRAM_CHAT_ID, "text": message, "parse_mode": "HTML", "disable_web_page_preview": True}
    requests.post(url, data=data)

# --- MAIN ---
def main():
    print("Запуск агента v3.2 (Gemini + AutoFix)...")
    
    seen_ids = load_history()
    all_papers = []
    new_seen_ids = []

    # 1. Свежее
    print("Этап 1: Поиск свежих...")
    for category, query_list in RAW_QUERIES.items():
        for q in query_list:
            ids = search_pubmed(q, days=1, retmax=2)
            unique_ids = [i for i in ids if i not in seen_ids]
            if unique_ids:
                details = fetch_details_and_analyze(unique_ids)
                for paper in details:
                    paper['category'] = category
                    paper['type'] = 'fresh'
                    all_papers.append(paper)
                    seen_ids.add(paper['id'])
                    new_seen_ids.append(paper['id'])
            time.sleep(1)

    # 2. Архив
    if len(all_papers) < 10:
        print("Этап 2: Поиск в архиве...")
        needed = 10 - len(all_papers)
        for category, query_list in RAW_QUERIES.items():
            if needed <= 0: break
            for q in query_list:
                ids = search_pubmed(q, days=1825, retmax=5, sort="relevance")
                candidates = [i for i in ids if i not in seen_ids]
                if candidates:
                    to_take = candidates[:1]
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
        print("Ничего нового.")
        return

    all_papers.sort(key=lambda x: x['type'], reverse=True)

    # 3. ОТПРАВКА
    buffer_message = "<b>🧠 Biohack Digest (AI)</b>\n\n"
    current_category = ""
    
    for paper in all_papers:
        article_text = ""
        if paper['category'] != current_category:
            article_text += f"<b>🔹 {paper['category']}</b>\n"
            current_category = paper['category']
        
        icon = "⚡️" if paper['type'] == 'fresh' else "🔬"
        
        clean_summary = html.escape(paper['summary'])
        # Чистим Markdown, который иногда любит Gemini
        clean_summary = clean_summary.replace("**", "").replace("##", "")
        
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

if __name__ == "__main__":
    main()

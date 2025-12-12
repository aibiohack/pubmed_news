import requests
import time
import os
import html
from Bio import Entrez
import google.generativeai as genai # <--- Библиотека Gemini

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
    # Используем быструю и экономную модель Flash
    model = genai.GenerativeModel('gemini-1.5-flash')

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

# --- МОДУЛЬ АНАЛИЗА (GEMINI) ---
def analyze_abstract_with_gemini(title, abstract):
    """Отправляет абстракт в Gemini для извлечения сути на русском."""
    if not GEMINI_API_KEY:
        return f"Заголовок: {title} (Нет API ключа Gemini)"

    # Промпт (инструкция) для модели
    prompt = f"""
    Ты профессиональный спортивный физиолог. Проанализируй этот научный абстракт.
    
    Заголовок: {title}
    Текст: {abstract}

    Твоя задача:
    1. Перевести суть на русский язык.
    2. Сформулировать ОДНО емкое предложение по схеме:
       "✅ [Вмешательство/Добавка] ([Дозировка/Схема]) на [Кол-во людей] -> [Конкретный результат] (изменение на X% или p-value, если есть)."
    
    Если данных мало, просто напиши вывод. Не используй вводные слова типа "Исследование показало". Сразу к делу.
    """

    try:
        response = model.generate_content(prompt)
        return response.text.strip()
    except Exception as e:
        print(f"Ошибка Gemini: {e}")
        # Если ИИ ошибся, просто возвращаем английский заголовок, чтобы не потерять статью
        return f"Paper: {title} (Ошибка анализа)"

# --- ПОИСК PUBMED ---
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
                
                # Достаем абстракт
                abstract_parts = article['MedlineCitation']['Article'].get('Abstract', {}).get('AbstractText', [])
                full_abstract = " ".join(abstract_parts) if abstract_parts else ""

                if not full_abstract:
                    summary = f"Заголовок: {title_en} (Нет текста статьи)"
                else:
                    # GEMINI АНАЛИЗ
                    summary = analyze_abstract_with_gemini(title_en, full_abstract)
                
                link = f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/"
                pub_date = article['MedlineCitation']['Article']['Journal']['JournalIssue']['PubDate']
                year = pub_date.get('Year', 'N/A')
                
                papers.append({'summary': summary, 'link': link, 'id': str(pmid), 'year': year})
            except Exception as e:
                print(f"Ошибка обработки статьи {pmid}: {e}")
                continue
        return papers
    except Exception as e:
        print(f"Ошибка скачивания деталей: {e}")
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
    print("Запуск агента v3.1 (Gemini Analysis)...")
    
    seen_ids = load_history()
    all_papers = []
    new_seen_ids = []

    # 1. Свежее (берем по 2 статьи на тему, чтобы не перегрузить)
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
            time.sleep(1) # Чуть большая пауза для Gemini API (Rate limits)

    # 2. Архив (если мало)
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
    buffer_message = "<b>🧠 Biohack Digest (by Gemini)</b>\n\n"
    current_category = ""
    
    for paper in all_papers:
        article_text = ""
        if paper['category'] != current_category:
            article_text += f"<b>🔹 {paper['category']}</b>\n"
            current_category = paper['category']
        
        icon = "⚡️" if paper['type'] == 'fresh' else "🔬"
        
        # summary уже очищаем для HTML
        clean_summary = html.escape(paper['summary'])
        # Жирный шрифт для ключевых цифр (Gemini иногда ставит звездочки для Markdown, уберем их)
        clean_summary = clean_summary.replace("**", "")
        
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
        print(f"Сохранено {len(new_seen_ids)} статей.")

if __name__ == "__main__":
    main()

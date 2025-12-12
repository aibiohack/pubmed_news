import requests
import time
import os
from datetime import datetime
from Bio import Entrez

# --- НАСТРОЙКИ ---
Entrez.email = "energy17429@gmail.com"  # <--- ЗАМЕНИ НА СВОЙ EMAIL

# Берем из секретов или переменных
TELEGRAM_TOKEN = os.environ.get("TG_TOKEN")
TELEGRAM_CHAT_ID = os.environ.get("TG_CHAT_ID")

# Файл памяти (создастся сам)
HISTORY_FILE = "history.txt"

# Жесткий фильтр типов статей (только наука высокого качества)
QUALITY_FILTER = " AND (Meta-Analysis[ptyp] OR Randomized Controlled Trial[ptyp] OR Systematic Review[ptyp])"

# Твои категории и запросы
# Мы будем автоматически добавлять к ним фильтр качества при поиске
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

# --- МОДУЛЬ ПАМЯТИ ---
def load_history():
    """Загружает список ID статей, которые мы уже видели."""
    if not os.path.exists(HISTORY_FILE):
        return set()
    with open(HISTORY_FILE, "r") as f:
        return set(line.strip() for line in f)

def save_history(new_ids):
    """Дописывает новые ID в файл."""
    with open(HISTORY_FILE, "a") as f:
        for pmid in new_ids:
            f.write(f"{pmid}\n")

# --- ПОИСК ---
def search_pubmed(query, days=None, retmax=5, sort="date"):
    """Ищет статьи. Если days=None, ищет без ограничения по дате (но топ релевантных)."""
    # Добавляем фильтр качества к запросу
    full_query = query + QUALITY_FILTER
    
    try:
        params = {
            "db": "pubmed",
            "term": full_query,
            "retmax": retmax,
            "sort": sort
        }
        if days:
            params["reldate"] = days
            params["datetype"] = "pdat" # Дата публикации
        
        handle = Entrez.esearch(**params)
        record = Entrez.read(handle)
        handle.close()
        return record["IdList"]
    except Exception as e:
        print(f"Ошибка поиска '{query}': {e}")
        return []

def fetch_details(id_list):
    if not id_list: return []
    ids = ",".join(id_list)
    try:
        handle = Entrez.efetch(db="pubmed", id=ids, retmode="xml")
        records = Entrez.read(handle)
        handle.close()
        papers = []
        for article in records['PubmedArticle']:
            try:
                title = article['MedlineCitation']['Article']['ArticleTitle']
                pmid = article['MedlineCitation']['PMID']
                link = f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/"
                # Пробуем найти год публикации
                pub_date = article['MedlineCitation']['Article']['Journal']['JournalIssue']['PubDate']
                year = pub_date.get('Year', 'N/A')
                
                papers.append({'title': title, 'link': link, 'id': str(pmid), 'year': year})
            except:
                continue
        return papers
    except:
        return []

# --- TELEGRAM ---
def send_telegram_message(message):
    if not TELEGRAM_TOKEN or not TELEGRAM_CHAT_ID:
        print("Telegram токены не настроены. Пропуск отправки.")
        return
    url = f"https://api.telegram.org/bot{TELEGRAM_TOKEN}/sendMessage"
    data = {
        "chat_id": TELEGRAM_CHAT_ID,
        "text": message,
        "parse_mode": "HTML",
        "disable_web_page_preview": True
    }
    requests.post(url, data=data)

# --- ГЛАВНАЯ ЛОГИКА ---
def main():
    print("Запуск агента v2.0...")
    seen_ids = load_history()
    all_papers = []
    new_seen_ids = []

    # 1. Попытка найти СВЕЖЕЕ (за 24 часа) - самое важное
    print("Этап 1: Поиск свежих статей за 24 часа...")
    for category, query_list in RAW_QUERIES.items():
        for q in query_list:
            ids = search_pubmed(q, days=1, retmax=3)
            # Фильтруем то, что уже видели (маловероятно за 24ч, но вдруг)
            unique_ids = [i for i in ids if i not in seen_ids]
            
            if unique_ids:
                details = fetch_details(unique_ids)
                for paper in details:
                    paper['category'] = category
                    paper['type'] = 'fresh'
                    all_papers.append(paper)
                    seen_ids.add(paper['id'])
                    new_seen_ids.append(paper['id'])
            time.sleep(0.3)

    # 2. Если набрали меньше 15 статей, добиваем "Золотым фондом" (лучшее за 5 лет)
    # Но показываем только то, чего не было в истории
    if len(all_papers) < 15:
        print("Мало свежего. Этап 2: Поиск в архиве (5 лет)...")
        needed = 20 - len(all_papers)
        
        for category, query_list in RAW_QUERIES.items():
            if needed <= 0: break
            for q in query_list:
                # Ищем топ-10 самых релевантных за 5 лет (reldate=1825 дней)
                ids = search_pubmed(q, days=1825, retmax=10, sort="relevance")
                
                # Самое важное: берем только те ID, которых НЕТ в seen_ids
                candidates = [i for i in ids if i not in seen_ids]
                
                if candidates:
                    # Берем по 1-2 штуки с запроса, чтобы не забить все одной темой
                    to_take = candidates[:1] 
                    details = fetch_details(to_take)
                    for paper in details:
                        paper['category'] = category
                        paper['type'] = 'archive'
                        all_papers.append(paper)
                        seen_ids.add(paper['id'])
                        new_seen_ids.append(paper['id'])
                        needed -= 1
                time.sleep(0.3)

    # 3. Отправка и сохранение
    if not all_papers:
        print("Ничего нового не найдено.")
        # Можно не отправлять сообщение в ТГ, чтобы не спамить "пустотой"
        return

    # Сортировка: сначала свежие, потом архив
    all_papers.sort(key=lambda x: x['type'], reverse=True) 

    message = "<b>🧬 Biohack Daily Digest</b>\n"
    message += "<i>Только РКИ, Мета-анализы и Обзоры</i>\n\n"
    
    current_category = ""
    for paper in all_papers:
        if paper['category'] != current_category:
            message += f"<b>🔹 {paper['category']}</b>\n"
            current_category = paper['category']
        
        icon = "🔥" if paper['type'] == 'fresh' else "📚" # Огонь для новых, Книги для архива
        title = paper['title'].replace("<", "").replace(">", "")
        year = paper.get('year', '')
        
        message += f"{icon} <a href='{paper['link']}'>{title}</a> ({year})\n\n"

    # Режем сообщение если длинное
    chunks = [message[i:i+4096] for i in range(0, len(message), 4096)]
    for chunk in chunks:
        send_telegram_message(chunk)

    # СОХРАНЯЕМ ИСТОРИЮ
    if new_seen_ids:
        save_history(new_seen_ids)
        print(f"Сохранено {len(new_seen_ids)} новых статей в историю.")

if __name__ == "__main__":
    main()
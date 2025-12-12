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
QUALITY_FILTER = " AND (Meta-Analysis[ptyp] OR Randomized Controlled Trial[ptyp] OR Systematic Review[ptyp])"

# Настройка Gemini
if GEMINI_API_KEY:
    genai.configure(api_key=GEMINI_API_KEY)

RAW_QUERIES = {
    "🧬 1. Фундамент: Метаболизм и Долголетие": [
        # Митохондрии и Энергия
        "(mitochondrial biogenesis) AND (exercise) OR (zone 2)",
        "(NAD+) OR (NMN) OR (NR) AND (aging) OR (muscle performance)",
        "(AMPK activation) OR (sirtuins) AND (fasting) OR (exercise)",
        "(metabolic flexibility) AND (fat oxidation) OR (insulin sensitivity)",
        "(autophagy) AND (exercise) OR (time-restricted eating)",
        "(hormesis) AND (sauna) OR (cold exposure) OR (hypoxia)"
    ],

    "💊 2. Биохимия: Нутрицевтика и Адаптогены": [
        # База
        "(creatine monohydrate) AND (brain) OR (muscle) OR (recovery)",
        "(magnesium) AND (sleep) OR (muscle relaxation) OR (stress)",
        "(omega-3) OR (fish oil) AND (inflammation) OR (recovery) OR (concussion)",
        "(vitamin D) AND (athletic performance) OR (strength) OR (testosterone)",
        
        # Эргогены (Усилители)
        "(caffeine) AND (endurance) OR (power) OR (cognitive performance)",
        "(beta-alanine) AND (high intensity) OR (tactical athlete)",
        "(nitrate) OR (beetroot juice) AND (blood flow) OR (efficiency)",
        "(ketogenic diet) OR (exogenous ketones) AND (metabolism) OR (endurance)",
        
        # Адаптогены (Травы и Грибы)
        "(Ashwagandha) AND (cortisol) OR (strength) OR (testosterone)",
        "(Rhodiola rosea) AND (fatigue) OR (mental performance)",
        "(Lion's mane) OR (Hericium) AND (nerve growth factor) OR (cognition)",
        "(Cordyceps) AND (VO2max) OR (ATP production)",
        "(Tongkat Ali) OR (Eurycoma) AND (hormonal profile) OR (stress)",
        "(Shilajit) AND (mitochondria) OR (muscle strength)"
    ],

    "💪 3. Тело: Сила, Гипертрофия и Механика": [
        # Мышцы
        "(hypertrophy) AND (volume) OR (frequency) OR (mechanical tension)",
        "(eccentric training) AND (tendon) OR (strength gains)",
        "(blood flow restriction) OR (BFR) AND (rehabilitation) OR (growth)",
        
        # Взрывная сила и Плиометрика
        "(\"rate of force development\") OR (RFD) AND (explosive strength)",
        "(plyometric training) AND (sprint speed) OR (jumping)",
        "(\"stretch-shortening cycle\") AND (performance) OR (efficiency)",
        "(velocity-based training) AND (power) OR (autoregulation)",
        
        # Связки и профилактика
        "(tendon stiffness) AND (injury prevention) OR (energy return)",
        "(mobility) OR (flexibility) AND (performance) NOT (elderly)"
    ],

    "🫁 4. Двигатель: Кардио и Дыхание": [
        "(VO2max) AND (longevity) OR (performance)",
        "(heart rate variability) OR (HRV) AND (recovery) OR (readiness)",
        "(stroke volume) OR (cardiac output) AND (athlete's heart)",
        "(respiratory muscle training) OR (IMT) AND (endurance) OR (breathing)",
        "(nasal breathing) OR (mouth taping) AND (sleep) OR (exercise)"
    ],

    "🧠 5. Разум: Нейроатлетика и Психология": [
        # Психология победителя
        "(\"flow state\") AND (sport) OR (peak performance)",
        "(mental toughness) OR (resilience) AND (anxiety) OR (competition)",
        "(visualization) OR (motor imagery) AND (strength) OR (skill)",
        "(self-talk) AND (endurance) OR (effort perception)",
        
        # Нейрофизиология
        "(neuroplasticity) AND (exercise) OR (motor learning)",
        "(stroboscopic training) OR (visual training) AND (reaction time)",
        "(dopamine) AND (exercise motivation) OR (reward system)",
        "(transcranial direct current stimulation) OR (tDCS) AND (sport)"
    ],

    "💤 6. Восстановление: Сон и Ритмы": [
        "(circadian rhythm) OR (chronotype) AND (performance)",
        "(slow wave sleep) OR (deep sleep) AND (physical recovery)",
        "(REM sleep) AND (motor memory) OR (mental health)",
        "(sleep deprivation) AND (testosterone) OR (injury risk)",
        "(glymphatic system) AND (sleep) OR (brain clearance)"
    ],

    "🍃 7. Среда и Технологии": [
        "(wearable technology) AND (accuracy) OR (load monitoring)",
        "(blue light) AND (sleep quality) OR (alertness)",
        "(grounding) OR (earthing) AND (inflammation) OR (recovery)",
        "(music) OR (binaural beats) AND (focus) OR (relaxation)",
        "(continuous glucose monitoring) AND (athlete) OR (fueling)"
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
    if not GEMINI_API_KEY:
        print("❌ ОШИБКА: Нет GEMINI_API_KEY")
        return "⚠️ Нет ключа Gemini"

    prompt = f"""
    You are a sports physiologist. Analyze this abstract.
    Title: {title}
    Abstract: {abstract}
    
    Task:
    1. Summarize the key finding in ONE sentence in RUSSIAN.
    2. Format: "✅ [Action/Supplement] on [Subjects] -> [Result] (change % or value)."
    """

    # ОБНОВЛЕННЫЙ СПИСОК МОДЕЛЕЙ
    # Убрали старую 'gemini-pro', оставили только актуальные версии 1.5
    models_to_try = ['gemini-1.5-flash', 'gemini-1.5-pro']

    # Настройки безопасности (разрешаем всё, чтобы не блочило медицинские термины)
    safety_settings = {
        HarmCategory.HARM_CATEGORY_HARASSMENT: HarmBlockThreshold.BLOCK_NONE,
        HarmCategory.HARM_CATEGORY_HATE_SPEECH: HarmBlockThreshold.BLOCK_NONE,
        HarmCategory.HARM_CATEGORY_SEXUALLY_EXPLICIT: HarmBlockThreshold.BLOCK_NONE,
        HarmCategory.HARM_CATEGORY_DANGEROUS_CONTENT: HarmBlockThreshold.BLOCK_NONE,
    }

    last_error = ""

    for model_name in models_to_try:
        try:
            # print(f"🔍 Пробуем модель: {model_name}...") # Можно раскомментировать для отладки
            model = genai.GenerativeModel(model_name)
            
            response = model.generate_content(
                prompt, 
                safety_settings=safety_settings
            )
            
            if response.text:
                return response.text.strip()
            
        except Exception as e:
            # Если модель не сработала, молча пробуем следующую
            last_error = str(e)
            continue
    
    # Если перепробовали все и ничего не вышло
    return f"Заголовок: {title} (Ошибка AI: {last_error})"
    
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
    print("Запуск агента v3.3 (DEBUG MODE)...")
    
    seen_ids = load_history()
    all_papers = []
    new_seen_ids = []

    # 1. Свежее
    print("Этап 1: Поиск свежих...")
    for category, query_list in RAW_QUERIES.items():
        for q in query_list:
            ids = search_pubmed(q, days=1, retmax=1) # Берем по 1 для быстрого теста
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
    if len(all_papers) < 5:
        print("Этап 2: Поиск в архиве...")
        needed = 5 - len(all_papers)
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
    buffer_message = "<b>🧠 Biohack Digest (DEBUG)</b>\n\n"
    current_category = ""
    
    for paper in all_papers:
        article_text = ""
        if paper['category'] != current_category:
            article_text += f"<b>🔹 {paper['category']}</b>\n"
            current_category = paper['category']
        
        icon = "⚡️" if paper['type'] == 'fresh' else "🔬"
        
        clean_summary = html.escape(paper['summary'])
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

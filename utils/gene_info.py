import requests
from bs4 import BeautifulSoup
from googletrans import Translator
from concurrent.futures import ThreadPoolExecutor, as_completed

headers = {
    "User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64)"
}
BASE_URL = "http://wheatomics.sdau.edu.cn/cgi-bin/geneDetail.py?search="

def fetch_single_gene_raw(gene_id):
    """底层函数：抓取单个基因信息"""
    gid = gene_id.strip()
    if not gid:
        return gid, "空白", "空白"
    try:
        r = requests.get(BASE_URL + gid, headers=headers, timeout=10)
        r.encoding = "utf-8"
        soup = BeautifulSoup(r.text, "html.parser")
        
        # 三代 ID
        third_id = "未找到"
        table = soup.find("table", {"id": "genetable"})
        if table:
            td = table.find("td")
            if td:
                lines = list(td.stripped_strings)
                if len(lines) >= 3:
                    third_id = lines[2].strip()

        # Description
        desc_en = "未找到"
        for row in soup.find_all("tr"):
            th, td = row.find("th"), row.find("td")
            if th and td and th.text.strip() == "Description:":
                desc_en = td.text.strip()
                break
        return gid, third_id, desc_en
    except:
        return gid, "异常", "异常"

def translate_text(text):
    """底层函数：翻译"""
    if text in ("未找到", "异常", "空白"):
        return text
    try:
        translator = Translator()
        return translator.translate(text, src="en", dest="zh-cn").text
    except:
        return "翻译失败"

# ---------------- 导出给 app.py 使用的接口 ----------------

def fetch_gene_info(gene_id):
    """保持原签名，兼容旧逻辑"""
    _, tid, desc = fetch_single_gene_raw(gene_id)
    return tid, desc

def translate_description(text):
    """保持原签名，兼容旧逻辑"""
    return translate_text(text)

def batch_process_gene_info(gene_ids, max_workers=10, progress_callback=None):
    """新增加的批量处理接口"""
    total = len(gene_ids)
    results_map = {}
    
    # 阶段 1: 并发爬取
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        f_to_g = {executor.submit(fetch_single_gene_raw, g): g for g in gene_ids}
        for i, f in enumerate(as_completed(f_to_g), 1):
            gid, tid, desc = f.result()
            results_map[gid] = {"tid": tid, "en": desc}
            if progress_callback:
                progress_callback(i / (total * 2), f"🧬 爬取中: {gid} ({i}/{total})")

    # 阶段 2: 并发翻译
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        f_to_t = {executor.submit(translate_text, info["en"]): g 
                  for g, info in results_map.items()}
        for i, f in enumerate(as_completed(f_to_t), 1):
            gid = f_to_t[f]
            results_map[gid]["zh"] = f.result()
            if progress_callback:
                progress_callback((total + i) / (total * 2), f"🌏 翻译中: {gid} ({i}/{total})")

    return [[g, results_map[g]["tid"], results_map[g]["en"], results_map[g]["zh"]] for g in gene_ids]
import requests
from bs4 import BeautifulSoup
from googletrans import Translator

headers = {
    "User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64)"
}
base_url = "http://wheatomics.sdau.edu.cn/cgi-bin/geneDetail.py?search="

# 全局 translator（避免每次都重新初始化，速度会更快）
translator = Translator()


def fetch_gene_info(gene_id):
    """获取三代基因号 + 英文功能描述"""
    try:
        r = requests.get(base_url + gene_id, headers=headers, timeout=10)
        r.encoding = "utf-8"
        if r.status_code != 200:
            return "未找到", "未找到"
    except Exception:
        return "未找到", "未找到"

    soup = BeautifulSoup(r.text, "html.parser")

    # 三代基因号
    third_id = "未找到"
    table = soup.find("table", {"id": "genetable"})
    if table:
        td = table.find("td")
        if td:
            gene_lines = list(td.stripped_strings)
            if len(gene_lines) >= 3:
                third_id = gene_lines[2].strip()

    # Description
    description = "未找到"
    for row in soup.find_all("tr"):
        th, td = row.find("th"), row.find("td")
        if th and td and th.text.strip() == "Description:":
            description = td.text.strip()
            break

    return third_id, description


def translate_description(text):
    """翻译英文描述 -> 中文"""
    if text == "未找到":
        return "未找到"
    try:
        translation = translator.translate(text, src="en", dest="zh-cn")
        return translation.text
    except Exception as e:
        return f"翻译失败: {e}"

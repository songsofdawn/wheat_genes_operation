import requests
from bs4 import BeautifulSoup
import re
import textwrap
import time
import random

BASE_URL = "http://wheatomics.sdau.edu.cn/cgi-bin/get_fasta_bedtools.py?database=all_gene&ID="

HEADERS = {"User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64)"}

def wrap_sequence(seq, width=60):
    """把序列按指定宽度换行"""
    return "\n".join(textwrap.wrap(seq, width))

def fetch_cdna_and_cds(gene_id):
    """返回 cDNA 和 CDS 序列，未获取到返回 'NA'"""
    if not gene_id.endswith(".1"):
        gene_id += ".1"
    url = BASE_URL + gene_id

    for attempt in range(3):
        try:
            r = requests.get(url, headers=HEADERS, timeout=20)
            r.raise_for_status()
            soup = BeautifulSoup(r.text, "html.parser")
            seq_div = soup.find("div", {"id": "seq"})
            if not seq_div:
                return "NA", "NA"

            content = seq_div.get_text("\n").strip()
            lines = [line.strip() for line in content.splitlines() if line.strip()]
            if len(lines) < 2:
                return "NA", "NA"

            header = lines[0]
            seq = "".join(lines[1:])
            cdna_seq = wrap_sequence(seq)

            # 提取 CDS
            cds_match = re.search(r"CDS=(\d+)-(\d+)", header)
            if cds_match:
                start, end = map(int, cds_match.groups())
                cds_seq = wrap_sequence(seq[start - 1:end])
            else:
                cds_seq = "NA"

            return cdna_seq, cds_seq

        except Exception as e:
            wait = random.uniform(2, 5) * (attempt + 1)
            time.sleep(wait)

    return "NA", "NA"

import requests
from bs4 import BeautifulSoup
import re
import textwrap
import time
import random

# 标准遗传密码子表
CODON_TABLE = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M', 'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K', 'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L', 'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q', 'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V', 'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E', 'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S', 'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'_', 'TAG':'_', 'TGC':'C', 'TGT':'C', 'TGA':'_', 'TGG':'W',
}

def translate_cds(cds_raw_seq):
    """内部函数：将原始碱基序列翻译为蛋白质序列"""
    if not cds_raw_seq or cds_raw_seq == "NA":
        return "NA"
    cds = cds_raw_seq.upper().replace("\n", "").replace("\r", "").strip()
    protein = []
    for i in range(0, len(cds) - 2, 3):
        codon = cds[i:i+3]
        aa = CODON_TABLE.get(codon, 'X') # 无法识别则显示 X
        if aa == '_': # 遇到终止密码子停止
            break
        protein.append(aa)
    return "".join(protein)

def wrap_sequence(seq, width=60):
    """格式化序列换行"""
    if seq == "NA": return "NA"
    return "\n".join(textwrap.wrap(seq, width))

def fetch_cdna_cds_protein(gene_id):
    """
    适配 Streamlit APP 的主函数
    返回: cdna_seq, cds_seq, protein_seq (均为格式化后的字符串)
    """
    BASE_URL = "http://wheatomics.sdau.edu.cn/cgi-bin/get_fasta_bedtools.py?database=all_gene&ID="
    HEADERS = {"User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64)"}
    
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
                return "NA", "NA", "NA"

            content = seq_div.get_text("\n").strip()
            lines = [line.strip() for line in content.splitlines() if line.strip()]
            if len(lines) < 2:
                return "NA", "NA", "NA"

            header = lines[0]
            raw_full_seq = "".join(lines[1:])
            
            # 1. 提取 cDNA
            cdna_formatted = wrap_sequence(raw_full_seq)

            # 2. 提取并翻译 CDS
            cds_match = re.search(r"CDS=(\d+)-(\d+)", header)
            if cds_match:
                start, end = map(int, cds_match.groups())
                raw_cds = raw_full_seq[start - 1 : end]
                cds_formatted = wrap_sequence(raw_cds)
                
                # 翻译蛋白质
                raw_protein = translate_cds(raw_cds)
                protein_formatted = wrap_sequence(raw_protein)
            else:
                cds_formatted = "NA"
                protein_formatted = "NA"

            return cdna_formatted, cds_formatted, protein_formatted

        except Exception:
            time.sleep(random.uniform(1, 3))
            
    return "NA", "NA", "NA"
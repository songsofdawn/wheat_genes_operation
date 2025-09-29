import os
import requests
from bs4 import BeautifulSoup
import re

UPSTREAM_LEN = 2000

# 序列反向互补
def reverse_complement(seq):
    complement = str.maketrans("ACGTacgt", "TGCAtgca")
    return seq.translate(complement)[::-1]

# 从 GFF 文件提取基因坐标
def extract_coordinates(gff_file, gene_ids):
    coords = []
    gene_set = set(gene_ids)
    if not os.path.exists(gff_file):
        raise FileNotFoundError(f"GFF 文件不存在: {gff_file}")

    with open(gff_file, 'r', encoding='utf-8') as gff:
        for line in gff:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.strip().split("\t")
            if len(parts) < 9:
                continue
            chrom, feature_type, start, end, strand, attributes = parts[0], parts[2], parts[3], parts[4], parts[6], parts[8]
            if feature_type != "mRNA":
                continue
            attr_dict = dict(item.split("=") for item in attributes.split(";") if "=" in item)
            mrna_id = attr_dict.get("ID", "")
            if mrna_id in gene_set:
                coords.append((mrna_id, chrom, int(start), int(end), strand))
    return coords

# 构建上游启动子查询
def build_queries(coords, upstream_len=UPSTREAM_LEN):
    queries, records = [], []
    for gene_id, chrom, start, end, strand in coords:
        if strand == '+':
            q_start = max(0, start - upstream_len)
            q_end = start - 1
        else:
            q_start = end + 1
            q_end = end + upstream_len
        coord_str = f"{chrom}_Fielder:{q_start}-{q_end}"
        queries.append(coord_str)
        records.append((gene_id, strand))
    return queries, records

# 请求 Fielder 获取启动子序列
def fetch_promoters(queries):
    if not queries:
        return None
    url = "http://wheatomics.sdau.edu.cn/cgi-bin/get_fasta_bedtools.py"
    headers = {
        "User-Agent": "Mozilla/5.0",
        "Content-Type": "application/x-www-form-urlencoded"
    }
    data = {
        "database": "Fielder.genome",
        "ID": "\n".join(queries)
    }
    try:
        response = requests.post(url, headers=headers, data=data, timeout=30)
        response.encoding = 'utf-8'
        soup = BeautifulSoup(response.text, 'html.parser')
        seq_div = soup.find('div', id='seq')
        if not seq_div:
            return None
        return seq_div.text.strip()
    except Exception as e:
        print(f"❌ 请求失败: {e}")
        return None

# 主功能：返回 TXT 格式启动子序列（负链反向互补处理）
def fetch_promoter_sequences(gene_ids, gff_file):
    coords = extract_coordinates(gff_file, gene_ids)
    txt_lines = []

    # 处理每个基因，如果找不到坐标，也要写入提示
    gene_coord_map = {gid: None for gid in gene_ids}
    for gene_id, chrom, start, end, strand in coords:
        gene_coord_map[gene_id] = (chrom, start, end, strand)

    queries, records = [], []
    for gene_id in gene_ids:
        val = gene_coord_map[gene_id]
        if val:
            chrom, start, end, strand = val
            if strand == '+':
                q_start = max(0, start - UPSTREAM_LEN)
                q_end = start - 1
            else:
                q_start = end + 1
                q_end = end + UPSTREAM_LEN
            queries.append(f"{chrom}_Fielder:{q_start}-{q_end}")
            records.append((gene_id, strand))
        else:
            # 如果没有找到坐标
            txt_lines.append(f">{gene_id} 未找到坐标\nNA\n")

    # 请求序列
    if queries:
        fasta_text = fetch_promoters(queries)
        if fasta_text:
            fasta_blocks = fasta_text.strip().split('>')
            fasta_blocks = [b for b in fasta_blocks if b.strip()]
            for block, (gene_id, strand) in zip(fasta_blocks, records):
                lines = block.splitlines()
                seq = "".join(lines[1:]).replace(" ", "").replace("\n", "")
                if strand == '-':
                    seq = reverse_complement(seq)
                txt_lines.append(f">{gene_id} ({strand} strand)\n{seq}\n")
        else:
            # 如果没有返回序列
            for gene_id, strand in records:
                txt_lines.append(f">{gene_id} ({strand} strand)\n未找到启动子序列\n")
    return "\n".join(txt_lines)

import requests
from bs4 import BeautifulSoup

BASE_URL_GENE = "http://wheatomics.sdau.edu.cn/cgi-bin/geneDetail.py?search="
BASE_URL_FASTA = "http://wheatomics.sdau.edu.cn/cgi-bin/get_fasta_bedtools.py"

# -------------------- 基因坐标获取 --------------------
def fetch_gene_info_wheatomics(gene_id):
    """
    返回格式：[gene_id, chrom, start, end, strand]
    """
    headers = {"User-Agent": "Mozilla/5.0"}
    try:
        response = requests.get(BASE_URL_GENE + gene_id, headers=headers, timeout=10)
        response.encoding = 'utf-8'
        soup = BeautifulSoup(response.text, 'html.parser')
        table = soup.find('table', id='genetable')
        if not table:
            return [gene_id, "NOTFOUND", "NOTFOUND", "NOTFOUND", "NOTFOUND"]

        location = strand = None
        for row in table.find_all('tr'):
            th = row.find('th')
            td = row.find('td')
            if not th or not td:
                continue
            key, val = th.text.strip(), td.text.strip()
            if key == 'Location:':
                location = val
            elif key == 'Strand:':
                strand = val

        if location and strand:
            chrom_part, pos_part = location.split(":")
            chrom = chrom_part.replace("Chinese_Spring1.0_", "")
            start, end = map(int, pos_part.replace(",", "").split(" - "))
            return [gene_id, chrom, str(start), str(end), strand]
        else:
            return [gene_id, "NA", "NA", "NA", "NA"]

    except Exception as e:
        print(f"[ERROR] 获取 {gene_id} 出错: {e}")
        return [gene_id, "ERROR", "ERROR", "ERROR", "ERROR"]

# -------------------- 构建查询坐标 --------------------
def build_query_list(gene_infos, upstream_len=2000):
    queries = []
    gene_info_list = []

    for info in gene_infos:
        gene_id, chrom, start, end, strand = info
        if "ERROR" in info or "NA" in info or "NOTFOUND" in info:
            continue
        start, end = int(start), int(end)
        if strand == "+":
            query_start = max(0, start - upstream_len - 1)
            query_end = start - 1
        elif strand == "-":
            query_start = end + 1
            query_end = end + upstream_len + 1
        else:
            continue
        coord = f"{chrom}_Chinese_Spring1.0:{query_start}-{query_end}"
        queries.append(coord)
        gene_info_list.append((gene_id, strand))
    return queries, gene_info_list

# -------------------- 获取启动子序列 --------------------
def fetch_promoters_from_wheatomics(queries):
    headers = {
        "User-Agent": "Mozilla/5.0",
        "Content-Type": "application/x-www-form-urlencoded"
    }
    data = {
        "database": "all_genomes",
        "ID": "\n".join(queries)
    }
    response = requests.post(BASE_URL_FASTA, headers=headers, data=data, timeout=20)
    response.encoding = 'utf-8'
    soup = BeautifulSoup(response.text, 'html.parser')
    seq_div = soup.find('div', id='seq')
    if not seq_div:
        raise RuntimeError("未找到返回的序列，请检查请求格式或网站状态")
    return seq_div.text.strip()

# -------------------- 反向互补 --------------------
def reverse_complement(seq):
    comp_table = str.maketrans("ACGTacgt", "TGCAtgca")
    return seq.translate(comp_table)[::-1]

# -------------------- 格式化FASTA --------------------
def format_fasta(fasta_text, gene_info_list):
    fasta_blocks = [b for b in fasta_text.strip().split('>') if b.strip()]
    records = []
    for block, (gene_id, strand) in zip(fasta_blocks, gene_info_list):
        lines = block.splitlines()
        seq = "".join(lines[1:]).replace(" ", "").replace("\n", "")
        if strand == "-":
            seq = reverse_complement(seq)
        # 每60个碱基换行
        formatted_seq = [seq[i:i+60] for i in range(0, len(seq), 60)]
        record = f">{gene_id} ({strand} strand)\n" + "\n".join(formatted_seq)
        records.append(record)
    return records

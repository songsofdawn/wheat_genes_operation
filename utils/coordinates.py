import requests
from bs4 import BeautifulSoup
import re

headers = {"User-Agent": "Mozilla/5.0 (Windows NT 10.0; Win64; x64)"}
base_url = "http://wheatomics.sdau.edu.cn/cgi-bin/geneDetail.py?search="
cds_url_template = (
    "http://wheatomics.sdau.edu.cn/cgi-bin/geneDetail_get_sequence.py?"
    "genome_db=all_genomes&chrom={chrom}_Chinese_Spring1.0&start={start}&end={end}"
    "&gene_db=all_gene&gene_id={gene_id}.1&protein_db=all_protein&protein_id={gene_id}.1"
)

def fetch_gene_coordinate(gene_id):
    try:
        r = requests.get(base_url + gene_id, headers=headers, timeout=10)
        r.encoding = "utf-8"
        soup = BeautifulSoup(r.text, "html.parser")
        table = soup.find("table", id="genetable")
        if not table:
            return [gene_id, "NOTFOUND", "NOTFOUND", "NOTFOUND", "NOTFOUND"]
        location = strand = None
        for row in table.find_all("tr"):
            th, td = row.find("th"), row.find("td")
            if not th or not td:
                continue
            if th.text.strip() == "Location:":
                location = td.text.strip()
            elif th.text.strip() == "Strand:":
                strand = td.text.strip()
        if location and strand:
            chrom_part, pos_part = location.split(":")
            chrom = chrom_part.replace("Chinese_Spring1.0_", "")
            start, end = map(int, pos_part.replace(",", "").split(" - "))
            return [gene_id, chrom, str(start), str(end), strand]
        else:
            return [gene_id, "NA", "NA", "NA", "NA"]
    except Exception:
        return [gene_id, "ERROR", "ERROR", "ERROR", "ERROR"]

def fetch_cds_fasta_block(gene_id, chrom, start, end):
    try:
        url = cds_url_template.format(chrom=chrom, start=start, end=end, gene_id=gene_id)
        r = requests.get(url, headers=headers, timeout=10)
        r.encoding = "utf-8"
        text = r.text
        pattern = re.compile(rf"(>{gene_id}\.1 gene={gene_id}.*?)(?=>\S|\Z)", re.DOTALL)
        match = pattern.search(text)
        if match:
            return match.group(1).strip()
    except Exception:
        return None
    return None

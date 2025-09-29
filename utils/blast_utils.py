import requests
from requests_toolbelt.multipart.encoder import MultipartEncoder
from bs4 import BeautifulSoup
import re
import time

blast_url = "http://wheatomics.sdau.edu.cn/blast/blastresult2.php"
blast_result_url_template = "http://wheatomics.sdau.edu.cn/blast/data/{}.blast2.html"
headers = {"User-Agent": "Mozilla/5.0"}

def submit_blast(seq):
    try:
        m = MultipartEncoder(
            fields={
                "querySeq": f">{seq[:10]}\n{seq}",
                "program": "blastn",
                "blastpath": "./blast+/bin",
                "group": "all_genes",
                "NWheat_IWGSC_RefSeq_v1.0_chromosomesdb[]": "Fielder_gene",
                "dbType": "N",
                "patientIDarrayChoix[]": "Fielder_gene",
                "patientIDarray[]": "Fielder_gene",
                "blast_flag": "1",
                "expect": "0.0001",
                "wordSize": "15",
                "targetSeqs": "50",
                "mmScore": "2,-3",
                "gapCost": "Existence: 5, Extension: 2",
                "filter": "T",
                "softMask": "m",
                "outFmt": "5",
                "OTHER_ADVANCED": "",
                "searchType": "basic",
            }
        )
        headers_local = {
            "User-Agent": "Mozilla/5.0",
            "Referer": "http://wheatomics.sdau.edu.cn/blast/blast.html",
            "Origin": "http://wheatomics.sdau.edu.cn",
            "Content-Type": m.content_type,
        }
        res = requests.post(blast_url, data=m, headers=headers_local, timeout=30)
        match = re.search(r"job id is (\d+)", res.text)
        if match:
            return match.group(1)
    except Exception:
        return None
    return None

def get_blast_result(jobid):
    try:
        for _ in range(10):
            res = requests.get(blast_result_url_template.format(jobid), headers=headers, timeout=10)
            if res.status_code == 200:
                soup = BeautifulSoup(res.text, "html.parser")
                tbody = soup.find("tbody")
                if tbody:
                    rows = tbody.find_all("tr")
                    if len(rows) >= 2:
                        tds = rows[1].find_all("td")
                        if tds and tds[-1].find("a"):
                            return tds[-1].find("a").text.strip()
            time.sleep(3)
    except Exception:
        return None
    return None

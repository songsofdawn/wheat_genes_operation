import time
import os
import pandas as pd
import streamlit as st

# 功能 1~4 导入
from utils.gene_info import fetch_gene_info, translate_description, batch_process_gene_info
from utils.sequence_fetcher import fetch_cdna_cds_protein
from utils.fasta_utils import parse_fasta
from utils.blast_utils import submit_blast, get_blast_result
from utils.fielder_promoter import fetch_promoter_sequences

# 新功能: 中国春启动子抓取
from utils.promoter_fetcher import (
    fetch_gene_info_wheatomics,
    build_query_list,
    fetch_promoters_from_wheatomics,
    format_fasta
)

# -------------------- 工具函数 --------------------
def read_gene_ids(uploaded_file, manual_input):
    if uploaded_file:
        return uploaded_file.read().decode("utf-8").splitlines()
    elif manual_input.strip():
        return manual_input.strip().splitlines()
    else:
        return []

# -------------------- 主程序 --------------------
def main():
    st.set_page_config(page_title="🌾 小麦基因批量处理工具", layout="wide")
    st.title("🌾 小麦基因批量处理工具")

    tool = st.sidebar.radio(
        "选择功能，初次使用请阅读ReadMe",
        [
            "ReadMe",
            "基因功能注释及三代基因号转换",
            "基因cDNA & CDS & protein sequences下载",
            "中国春 → Fielder 同源基因",
            "Fielder 基因 → 启动子序列",
            "中国春启动子抓取"
        ]
    )   

    # -------------------- ReadMe --------------------
    if tool == "ReadMe":
        st.header("📖 小麦基因批量处理工具 - 使用说明")
        st.markdown("""
        本网站基于WheatOmics开发，网址：http://wheatomics.sdau.edu.cn/
        **功能说明：**
        1. **基因功能注释及三代基因号转换**  
           因为是基于wheatomics，所以基因号不要加“ .1”。

        2. **基因 cDNA & CDS & protein下载**  
           上传基因号 TXT 文件（这个加不加.1无所谓），原理是基于wheatomics的GetSequence功能，可下载 cDNA 序列、CDS 序列和 protein 序列。

        3. **中国春 → Fielder 同源基因**  
           提交 CDS 序列 BLAST，可获取 Fielder 同源基因号。

        4. **Fielder 基因 → 启动子序列**  
           原理为：由基因号获取基因所在坐标，然后上溯2000bp，已考虑+-链问题，所得序列方向为：设TSS为0，从左到右为-2000到-1。

        **免责声明：**  
        - 请仔细核对网站结果，虽然是基于算法应该具有普适性，但是不排除代码出错，请仔细核对。

        **开发人员：**
        NWAFU科研楼2303wyz；GitHub：songofdawn（项目已开源）
        """)
        st.info("📌 请从左侧选择功能开始使用")
        return

    # -------------------- 功能 1 --------------------
    if tool == "基因功能注释及三代基因号转换":
            st.header("🔍 基因功能注释 (多线程加速版)")
            uploaded_file = st.file_uploader(
                "上传 TXT 文件（基因号一行一个，不要加“.1”）",
                type=["txt"], key="file_gene_info"
            )
            manual_input = st.text_area("或者手动输入基因号（每行一个）", key="input_gene_info")
            
            gene_ids = read_gene_ids(uploaded_file, manual_input)
            gene_ids = [g.strip() for g in gene_ids if g.strip()] 

            if not gene_ids:
                st.info("请上传文件或输入基因号")
                st.stop()

            if st.button("开始高速查询", key="btn_gene_info"):
                progress_bar = st.progress(0)
                status_text = st.empty()

                def update_ui(current_progress, text):
                    progress_bar.progress(current_progress)
                    status_text.text(text)

                # 调用新的批量处理函数
                results = batch_process_gene_info(
                    gene_ids, 
                    max_workers=10, 
                    progress_callback=update_ui
                )

                df = pd.DataFrame(results, columns=["输入基因号", "三代基因号", "功能描述（英文）", "功能描述（中文）"])
                st.success("✅ 查询完成！")
                st.dataframe(df, use_container_width=True)
                st.download_button("📥 下载结果 CSV", df.to_csv(index=False).encode("utf-8-sig"), "gene_info.csv", "text/csv")
                
                status_text.empty()
                progress_bar.empty()
    # -------------------- 功能 2 --------------------
    elif tool == "基因cDNA & CDS & protein sequences下载":
        st.header("📍 cDNA / CDS / Protein 下载")

        uploaded_file = st.file_uploader("上传 TXT 文件（基因号）", type=["txt"], key="file_sequences")
        manual_input = st.text_area("或者手动输入基因号（每行一个）", key="input_sequences")

        gene_ids = read_gene_ids(uploaded_file, manual_input)

        if not gene_ids:
            st.info("请上传文件或输入基因号")
            st.stop()

        if st.button("获取序列（cDNA / CDS / Protein）", key="btn_sequences"):

            cdna_records = []
            cds_records = []
            protein_records = []
            failed_genes = []

            progress = st.progress(0)
            status_text = st.empty()

            for idx, gene_id in enumerate(gene_ids, 1):
                gene_id = gene_id.strip()
                status_text.text(f"正在处理: {gene_id} ({idx}/{len(gene_ids)})")

                cdna_seq, cds_seq, protein_seq = fetch_cdna_cds_protein(gene_id)

                # ---------- cDNA ----------
                if cdna_seq == "NA":
                    cdna_records.append(f">{gene_id}\nNA\n")
                else:
                    cdna_records.append(f">{gene_id}\n{cdna_seq}\n")

                # ---------- CDS ----------
                if cds_seq == "NA":
                    cds_records.append(f">{gene_id}\nNA\n")
                else:
                    cds_records.append(f">{gene_id}\n{cds_seq}\n")

                # ---------- Protein ----------
                if protein_seq == "NA":
                    protein_records.append(f">{gene_id}\nNA\n")
                else:
                    protein_records.append(f">{gene_id}\n{protein_seq}\n")

                # ---------- 失败判断 ----------
                if cdna_seq == "NA" and cds_seq == "NA" and protein_seq == "NA":
                    failed_genes.append(gene_id)

                progress.progress(idx / len(gene_ids))

            # ---------- 下载按钮 ----------
            st.download_button(
                "📥 下载 cDNA FASTA",
                "\n".join(cdna_records),
                "cdna_sequences.fasta",
                "text/plain"
            )

            st.download_button(
                "📥 下载 CDS FASTA",
                "\n".join(cds_records),
                "cds_sequences.fasta",
                "text/plain"
            )

            st.download_button(
                "📥 下载 Protein FASTA ⭐",
                "\n".join(protein_records),
                "protein_sequences.fasta",
                "text/plain"
            )

            # ---------- 错误提示 ----------
            if failed_genes:
                st.warning(f"⚠️ 以下基因未获取到序列: {', '.join(failed_genes)}")

            status_text.empty()

# -------------------- 功能 3 --------------------
    elif tool == "中国春 → Fielder 同源基因":
        st.header("🧬 中国春基因号 → Fielder 同源基因 (BLAST)")
        
        # 1. 文件上传
        uploaded_file = st.file_uploader(
            "上传 CDS FASTA 或 TXT 文件", 
            type=["fasta", "fa", "txt"], 
            key="file_blast"
        )
        
        # 2. 文本输入（现在直接在下方）
        manual_input = st.text_area(
            "或者手动输入序列 (FASTA格式)", 
            height=200, 
            key="input_blast",
            placeholder=">GeneID\nATGC..."
        )

        # 整合输入源逻辑
        fasta_str = ""
        if uploaded_file:
            fasta_str = uploaded_file.read().decode("utf-8")
        elif manual_input.strip():
            fasta_str = manual_input.strip()

        if not fasta_str:
            st.info("💡 请上传文件或在上方输入框粘贴序列")
            st.stop()

        # 解析序列
        seq_records = parse_fasta(fasta_str)
        st.info(f"已解析到 {len(seq_records)} 条序列")
        
        if st.button("开始 BLAST", key="btn_blast"):
            results, progress, status_text = [], st.progress(0), st.empty()
            for idx, (name, seq) in enumerate(seq_records, 1):
                status_text.text(f"正在提交 BLAST: {name} ({idx}/{len(seq_records)})")
                jobid = submit_blast(seq)
                gene_id = get_blast_result(jobid) if jobid else "提交失败"
                results.append([name, gene_id if gene_id else "未找到"])
                
                # 更新进度
                progress.progress(idx / len(seq_records))
                # 适当延时防止被服务器屏蔽
                time.sleep(1) 
            
            df = pd.DataFrame(results, columns=["中国春基因", "Fielder 同源基因"])
            st.success("✅ BLAST 完成！")
            st.dataframe(df, use_container_width=True)
            st.download_button(
                "📥 下载结果 CSV", 
                df.to_csv(index=False).encode("utf-8-sig"), 
                "fielder_homologs.csv", 
                "text/csv"
            )
            status_text.empty()

 # -------------------- 功能 4 --------------------
    elif tool == "Fielder 基因 → 启动子序列":
        st.header("🌱 Fielder 基因 → 启动子序列抓取")
        
        uploaded_file = st.file_uploader("上传基因列表 TXT", type=["txt"], key="file_promoter")
        manual_input = st.text_area("或者手动输入基因号（每行一个）", key="input_promoter")
        
        gene_ids = read_gene_ids(uploaded_file, manual_input)
        gene_ids = [g.strip() for g in gene_ids if g.strip()]

        if not gene_ids:
            st.info("请上传文件或输入基因号")
            st.stop()

        if st.button("抓取启动子序列", key="btn_promoter"):
            BASE_DIR = os.path.dirname(os.path.abspath(__file__))
            
            # 建议直接指向压缩包名，或者保持 Fielder.gff 也可以（因为 utils 里有兼容逻辑）
            gff_path_default = os.path.join(BASE_DIR, "data", "Fielder.gff.gz") 

            if not os.path.exists(gff_path_default):
                st.error(f"GFF 压缩文件不存在: {gff_path_default}")
                st.stop()

            progress, status_text = st.progress(0), st.empty()
            txt_lines = []
            for idx, gene_id in enumerate(gene_ids, 1):
                status_text.text(f"正在抓取: {gene_id} ({idx}/{len(gene_ids)})")
                seq_txt = fetch_promoter_sequences([gene_id], gff_path_default)
                txt_lines.append(seq_txt)
                progress.progress(idx / len(gene_ids))
            
            st.download_button("📥 下载启动子 TXT", "\n".join(txt_lines), "promoter_sequences.txt", "text/plain")
            st.success("✅ 启动子抓取完成")
            status_text.empty()

# -------------------- 功能 5 --------------------
    elif tool == "中国春启动子抓取":
        st.header("🌱 中国春基因启动子抓取")
        
        uploaded_file = st.file_uploader("上传基因列表 TXT", type=["txt"], key="file_cs_promoter")
        manual_input = st.text_area("或者手动输入基因号（每行一个）", key="input_cs_promoter")
        upstream_len = st.number_input("上游长度(bp)", min_value=100, max_value=5000, value=2000, step=100)
        
        gene_ids = read_gene_ids(uploaded_file, manual_input)
        gene_ids = [g.strip() for g in gene_ids if g.strip()]

        if not gene_ids:
            st.info("请上传文件或输入基因号")
            st.stop()

        st.info(f"待处理基因数: {len(gene_ids)}")
        if st.button("开始抓取", key="btn_cs_promoter"):
            progress, status_text = st.progress(0), st.empty()
            gene_infos = []
            
            # 1. 获取坐标
            for idx, gene in enumerate(gene_ids, 1):
                status_text.text(f"获取坐标: {gene} ({idx}/{len(gene_ids)})")
                gene_infos.append(fetch_gene_info_wheatomics(gene))
                progress.progress((idx / len(gene_ids)) * 0.5) # 前 50% 进度
            
            # 2. 构建查询并下载
            queries, gene_info_list = build_query_list(gene_infos, upstream_len)
            fasta_records = []
            for idx, query in enumerate(queries, 1):
                status_text.text(f"抓取序列: {idx}/{len(queries)}")
                fasta_text = fetch_promoters_from_wheatomics([query])
                fasta_records.extend(format_fasta(fasta_text, [gene_info_list[idx-1]]))
                progress.progress(0.5 + (idx / len(queries)) * 0.5) # 后 50% 进度
            
            st.download_button("📥 下载启动子序列", "\n".join(fasta_records), "cs_promoter_sequences.txt", "text/plain")
            st.success("🎉 启动子抓取完成！")
            status_text.empty()


if __name__ == "__main__":
    main()

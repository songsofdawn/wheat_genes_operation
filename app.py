import os
import pandas as pd
import streamlit as st

from utils.db_query import (
    get_primary_gene_id,
    get_gene_annotation,
    get_gene_core,
    get_sequences,
    get_promoter,
    get_fielder_best_hit,
    get_fielder_all_hits,
    get_fielder_promoter,
    get_cs_self_best_hit,
    get_cs_self_all_hits
)

from utils.go_enrichment import run_go_enrichment, create_go_barplot_bytes
from utils.kegg_enrichment import (
    run_kegg_enrichment,
    create_kegg_bubbleplot_bytes,
    create_kegg_barplot_bytes
)

# -------------------- 工具函数 --------------------
def read_gene_ids(uploaded_file, manual_input):
    if uploaded_file:
        return uploaded_file.read().decode("utf-8").splitlines()
    elif manual_input.strip():
        return manual_input.strip().splitlines()
    else:
        return []

# -------------------- 打赏提示 --------------------
def show_tip_box():
    """
    在侧边栏显示打赏二维码。
    图片路径：
    data/TipCode.jpg
    """
    BASE_DIR = os.path.dirname(os.path.abspath(__file__))
    tip_code_path = os.path.join(BASE_DIR, "data", "TipCode.jpg")

    st.sidebar.markdown("---")
    with st.sidebar.expander("☕ 打赏支持 / Buy me a coffee", expanded=False):
        st.markdown(
            """
            如果这个工具帮到了你，可以扫码支持一下开发维护。
            """
        )

        if os.path.exists(tip_code_path):
            st.image(
                tip_code_path,
                caption="感谢支持 ❤️",
                use_container_width=True
            )
        else:
            st.warning("未找到打赏二维码图片，请检查路径：")
            st.code(tip_code_path)

# -------------------- 主程序 --------------------
def main():
    st.set_page_config(page_title="🌾 WheatGeneToolkit小麦基因批量处理工具", layout="wide")
    st.title("🌾 WheatGeneToolkit小麦基因批量处理工具")

    tool = st.sidebar.radio(
        "选择功能，初次使用请阅读ReadMe",
        [
            "ReadMe",
            "基因功能注释及三代基因号转换",
            "基因cDNA & CDS & protein sequences下载",
            "中国春同源基因检索（自身同源 + Fielder）",
            "Fielder 基因 → 启动子序列",
            "中国春启动子抓取",
            "GO富集分析",
            "KEGG富集分析"
        ]
    )
    show_tip_box()

    # -------------------- ReadMe --------------------
    if tool == "ReadMe":
        st.header("📖 小麦基因批量处理工具 - 使用说明")
        st.markdown("""
本工具面向小麦基因组相关分析，整合了**基因注释查询、基因号转换、序列下载、同源基因检索、启动子序列提取和 GO 富集分析、KEGG富集分析**等常用功能。  

---

## 功能说明

### 1. 基因功能注释及三代基因号转换
输入中国春基因号，可获取：
- 主键
- 三代基因号
- 基因功能描述（英文）

说明：
- 当前功能基于本地数据库查询
- 建议输入标准中国春基因号
- 一般不需要添加 `.1`

### 2. 基因 cDNA / CDS / protein sequences 下载
输入中国春基因号后，可批量获取：
- cDNA 序列
- CDS 序列
- protein 序列

说明：
- 当前序列来源于本地构建的中国春序列数据库
- 返回结果可能包含同一基因对应的多个 transcript 序列
- 下载文件为标准 FASTA 格式

### 3. 中国春同源基因检索（自身同源 + Fielder）
输入中国春基因号后，可同时查询：
- 中国春自身同源基因
- Fielder 同源基因

说明：
- 当前功能基于本地同源映射表，不再调用在线 BLAST
- 默认显示**最可信的对应基因**
- 同时保留全部候选同源基因，便于进一步人工判断
- 同源类型包括 `RBH` 和 `SBH`，其中 `RBH` 通常可信度更高

### 4. Fielder 基因 → 启动子序列
输入 Fielder 基因号后，可直接获取其启动子序列。

说明：
- 当前版本固定返回 **ATG 上游 2000 bp** 启动子序列
- 已考虑正负链方向
- 对于负链基因，输出序列已按转录方向进行反向互补处理
- 输出序列方向统一为启动子 5'→3' 方向,也就是说,你下的序列后面跟着的就是ATG所在的序列了

### 5. 中国春启动子抓取
输入中国春基因号后，可直接获取其启动子序列,同上

说明：
- 当前版本固定返回 **ATG 上游 2000 bp** 启动子序列
- 已考虑正负链方向
- 对于负链基因，输出序列已按转录方向进行反向互补处理
- 输出结果可直接下载为 FASTA 文件

### 6. GO富集分析
输入 DEG 列表（一行一个基因号），可进行 GO 富集分析，并输出：
- 富集结果表
- 显著结果表
- GO 富集条形图
- 分析摘要表

说明：
- 当前功能基于本地 GO 注释映射进行超几何检验和 BH 校正
- 推荐输入中国春 protein-coding gene
- 例如：`TraesCS1A02G000100`

### 7. KEGG 富集分析

输入 DEG 列表，一行一个基因号，可进行 KEGG 富集分析，并输出：

- KEGG 富集结果总表
- 显著富集通路表
- KEGG 富集气泡图
- KEGG 富集条形图
- KEGG 分析摘要表

说明：

- 当前功能基于本地 KEGG 注释文件，不依赖在线 KEGG 查询
- 使用本地 `gene-KO` 映射表将小麦基因转换为 KO 编号
- 再根据 `KO-pathway` 映射关系统计每个通路中的 KO 富集情况
- 富集检验使用超几何检验
- 支持自定义 P-value 阈值
- 支持设置最小和最大 pathway KO 数，过滤过小或过大的通路
- 支持绘制 KEGG 气泡图和条形图
- 推荐输入中国春基因号，例如：`TraesCS1A02G000100`

## 输入建议

- 基因号建议一行一个
- 中国春基因通常使用 `TraesCS...` 格式
- Fielder 基因通常使用 `TraesFLD...` 格式
- 若批量上传，建议使用 TXT 文件

---

## 注意事项

- 请仔细核对分析结果，本工具旨在提高批量处理效率，但不能替代人工审查
- 同源关系结果为计算推断结果，建议结合基因结构、共线性及研究背景进一步判断
- 启动子序列当前统一为 2000 bp 版本
- 若输入基因号格式不规范，可能导致查询失败

---

## 开发说明

本工具最初用于解决小麦基因组分析中重复查询、批量下载和多平台切换的问题。  
当前版本正在逐步从在线抓取逻辑重构为本地数据库驱动，以提高稳定性、查询速度和可维护性。

---

## 开发人员

NWAFU 科研楼 2303 wyz  
GitHub：songofdawn

---

## 免责声明

本工具结果仅供科研参考，请结合原始数据库、文献证据及具体研究背景进行综合判断。
        """)
        st.info("📌 请从左侧选择功能开始使用")
        return

    # -------------------- 功能 1 --------------------
    if tool == "基因功能注释及三代基因号转换":
        st.header("🔍 基因功能注释及三代基因号转换（本地数据库版）")

        uploaded_file = st.file_uploader(
            "上传 TXT 文件（基因号一行一个，不要加“.1”）",
            type=["txt"],
            key="file_gene_info"
        )
        manual_input = st.text_area(
            "或者手动输入基因号（每行一个）",
            key="input_gene_info"
        )

        gene_ids = read_gene_ids(uploaded_file, manual_input)
        gene_ids = [g.strip() for g in gene_ids if g.strip()]

        if not gene_ids:
            st.info("请上传文件或输入基因号")
            st.stop()

        st.info(f"待查询基因数：{len(gene_ids)}")

        if st.button("开始查询", key="btn_gene_info"):
            progress_bar = st.progress(0)
            status_text = st.empty()

            results = []

            for idx, input_gene_id in enumerate(gene_ids, 1):
                status_text.text(f"正在查询：{input_gene_id} ({idx}/{len(gene_ids)})")

                primary_gene_id = get_primary_gene_id(input_gene_id)
                anno_df = get_gene_annotation(input_gene_id)
                core_df = get_gene_core(input_gene_id)

                third_id = None
                description_en = None
                description_zh = None

                if not anno_df.empty:
                    description_en = anno_df.iloc[0].get("description_en", None)
                    description_zh = anno_df.iloc[0].get("description_zh", None)

                if not core_df.empty:
                    third_id = core_df.iloc[0].get("gene_id_v3", None)

                results.append([
                    input_gene_id,
                    primary_gene_id if primary_gene_id is not None else "未找到",
                    third_id if pd.notna(third_id) else "",
                    description_en if pd.notna(description_en) else "",
                    description_zh if pd.notna(description_zh) else ""
                ])

                progress_bar.progress(idx / len(gene_ids))

            df = pd.DataFrame(
                results,
                columns=["输入基因号", "统一主键", "三代基因号", "功能描述（英文）", "功能描述（中文）"]
            )

            st.success("✅ 查询完成！")
            st.dataframe(df, use_container_width=True)

            st.download_button(
                "📥 下载结果 CSV",
                df.to_csv(index=False).encode("utf-8-sig"),
                "gene_info_from_sqlite.csv",
                "text/csv"
            )

            status_text.empty()
            progress_bar.empty()

    # -------------------- 功能 2 --------------------
    elif tool == "基因cDNA & CDS & protein sequences下载":
        st.header("📍 cDNA / CDS / Protein 下载（本地数据库版）")

        uploaded_file = st.file_uploader("上传 TXT 文件（基因号）", type=["txt"], key="file_sequences")
        manual_input = st.text_area("或者手动输入基因号（每行一个）", key="input_sequences")

        gene_ids = read_gene_ids(uploaded_file, manual_input)
        gene_ids = [g.strip() for g in gene_ids if g.strip()]

        if not gene_ids:
            st.info("请上传文件或输入基因号")
            st.stop()

        st.info(f"待查询基因数：{len(gene_ids)}")

        if st.button("获取序列（cDNA / CDS / Protein）", key="btn_sequences"):
            cdna_records = []
            cds_records = []
            protein_records = []
            failed_genes = []
            summary_rows = []

            progress = st.progress(0)
            status_text = st.empty()

            for idx, input_gene_id in enumerate(gene_ids, 1):
                status_text.text(f"正在处理: {input_gene_id} ({idx}/{len(gene_ids)})")

                primary_gene_id = get_primary_gene_id(input_gene_id)
                seq_df = get_sequences(input_gene_id)

                if seq_df.empty:
                    failed_genes.append(input_gene_id)
                    summary_rows.append([input_gene_id, "", 0, 0, 0])
                    progress.progress(idx / len(gene_ids))
                    continue

                cdna_df = seq_df[seq_df["sequence_type"] == "cdna"].copy()
                cds_df = seq_df[seq_df["sequence_type"] == "cds"].copy()
                protein_df = seq_df[seq_df["sequence_type"] == "protein"].copy()

                summary_rows.append([
                    input_gene_id,
                    primary_gene_id if primary_gene_id else "",
                    len(cdna_df),
                    len(cds_df),
                    len(protein_df)
                ])

                for _, row in cdna_df.iterrows():
                    tx_id = row.get("transcript_id", "")
                    seq = row.get("sequence", "")
                    header = f">{input_gene_id}|{primary_gene_id}|{tx_id}|cdna"
                    cdna_records.append(f"{header}\n{seq}\n")

                for _, row in cds_df.iterrows():
                    tx_id = row.get("transcript_id", "")
                    seq = row.get("sequence", "")
                    header = f">{input_gene_id}|{primary_gene_id}|{tx_id}|cds"
                    cds_records.append(f"{header}\n{seq}\n")

                for _, row in protein_df.iterrows():
                    tx_id = row.get("transcript_id", "")
                    seq = row.get("sequence", "")
                    header = f">{input_gene_id}|{primary_gene_id}|{tx_id}|protein"
                    protein_records.append(f"{header}\n{seq}\n")

                if len(cdna_df) == 0 and len(cds_df) == 0 and len(protein_df) == 0:
                    failed_genes.append(input_gene_id)

                progress.progress(idx / len(gene_ids))

            summary_df = pd.DataFrame(
                summary_rows,
                columns=["输入基因号", "统一主键", "cDNA条数", "CDS条数", "Protein条数"]
            )

            st.success("✅ 序列查询完成！")
            st.dataframe(summary_df, use_container_width=True)

            st.download_button("📥 下载 cDNA FASTA", "".join(cdna_records), "cdna_sequences.fasta", "text/plain")
            st.download_button("📥 下载 CDS FASTA", "".join(cds_records), "cds_sequences.fasta", "text/plain")
            st.download_button("📥 下载 Protein FASTA ", "".join(protein_records), "protein_sequences.fasta", "text/plain")
            st.download_button(
                "📥 下载序列统计 CSV",
                summary_df.to_csv(index=False).encode("utf-8-sig"),
                "sequence_summary.csv",
                "text/csv"
            )

            if failed_genes:
                st.warning(f"⚠️ 以下基因未获取到序列: {', '.join(failed_genes)}")

            status_text.empty()
            progress.empty()

    # -------------------- 功能 3 --------------------
    elif tool == "中国春同源基因检索（自身同源 + Fielder）":
        st.header("🧬 中国春同源基因检索（自身同源 + Fielder）")

        uploaded_file = st.file_uploader("上传 TXT 文件（基因号一行一个）", type=["txt"], key="file_blast")
        manual_input = st.text_area("或者手动输入中国春基因号（每行一个）", height=200, key="input_blast")

        gene_ids = read_gene_ids(uploaded_file, manual_input)
        gene_ids = [g.strip() for g in gene_ids if g.strip()]

        if not gene_ids:
            st.info("💡 请上传文件或输入中国春基因号")
            st.stop()

        st.info(f"待查询基因数：{len(gene_ids)}")

        if st.button("开始查询同源基因", key="btn_blast"):
            progress = st.progress(0)
            status_text = st.empty()

            self_best_rows = []
            self_all_rows = []
            fld_best_rows = []
            fld_all_rows = []
            failed_genes = []

            for idx, input_gene_id in enumerate(gene_ids, 1):
                status_text.text(f"正在查询: {input_gene_id} ({idx}/{len(gene_ids)})")

                self_best_df = get_cs_self_best_hit(input_gene_id)
                self_all_df_one = get_cs_self_all_hits(input_gene_id)

                fld_best_df = get_fielder_best_hit(input_gene_id)
                fld_all_df_one = get_fielder_all_hits(input_gene_id)

                if self_best_df.empty and fld_best_df.empty:
                    failed_genes.append(input_gene_id)

                if not self_best_df.empty:
                    r = self_best_df.iloc[0]
                    self_best_rows.append([
                        input_gene_id,
                        r.get("cs_gene_id", ""),
                        r.get("self_homolog_gene_id", ""),
                        r.get("homolog_type", ""),
                        r.get("confidence", ""),
                        r.get("same_subgenome", ""),
                        r.get("same_chr_num", ""),
                        r.get("priority_score", "")
                    ])
                    tmp = self_all_df_one.copy()
                    tmp.insert(0, "input_gene_id", input_gene_id)
                    self_all_rows.append(tmp)

                if not fld_best_df.empty:
                    r = fld_best_df.iloc[0]
                    fld_best_rows.append([
                        input_gene_id,
                        r.get("cs_gene_id", ""),
                        r.get("fielder_gene_id", ""),
                        r.get("homolog_type", ""),
                        r.get("confidence", ""),
                        r.get("same_subgenome", ""),
                        r.get("same_chr_num", ""),
                        r.get("priority_score", "")
                    ])
                    tmp = fld_all_df_one.copy()
                    tmp.insert(0, "input_gene_id", input_gene_id)
                    fld_all_rows.append(tmp)

                progress.progress(idx / len(gene_ids))

            self_best_result_df = pd.DataFrame(
                self_best_rows,
                columns=["输入基因号", "中国春主键", "默认中国春自身同源基因", "同源类型", "可信度", "同亚基因组", "同染色体号", "优先级分数"]
            )

            fld_best_result_df = pd.DataFrame(
                fld_best_rows,
                columns=["输入基因号", "中国春主键", "默认Fielder同源基因", "同源类型", "可信度", "同亚基因组", "同染色体号", "优先级分数"]
            )

            self_all_df = pd.concat(self_all_rows, ignore_index=True) if self_all_rows else pd.DataFrame()
            fld_all_df = pd.concat(fld_all_rows, ignore_index=True) if fld_all_rows else pd.DataFrame()

            st.success("✅ 查询完成！")

            st.subheader("中国春自身同源：默认最可信结果")
            st.dataframe(self_best_result_df, use_container_width=True)

            st.subheader("中国春自身同源：全部候选")
            if not self_all_df.empty:
                st.dataframe(self_all_df, use_container_width=True)
                st.download_button(
                    "📥 下载中国春自身同源全部候选 CSV",
                    self_all_df.to_csv(index=False).encode("utf-8-sig"),
                    "cs_self_homolog_all_hits.csv",
                    "text/csv"
                )

            st.subheader("Fielder 同源：默认最可信结果")
            st.dataframe(fld_best_result_df, use_container_width=True)

            st.subheader("Fielder 同源：全部候选")
            if not fld_all_df.empty:
                st.dataframe(fld_all_df, use_container_width=True)
                st.download_button(
                    "📥 下载 Fielder 同源全部候选 CSV",
                    fld_all_df.to_csv(index=False).encode("utf-8-sig"),
                    "cs_to_fielder_all_hits.csv",
                    "text/csv"
                )

            if failed_genes:
                st.warning(f"⚠️ 以下基因未找到同源结果: {', '.join(failed_genes)}")

            status_text.empty()
            progress.empty()

    # -------------------- 功能 4 --------------------
    elif tool == "Fielder 基因 → 启动子序列":
        st.header("🌱 Fielder 基因 → 启动子序列（本地数据库版）")

        uploaded_file = st.file_uploader("上传 Fielder 基因列表 TXT", type=["txt"], key="file_fielder_promoter")
        manual_input = st.text_area("或者手动输入 Fielder 基因号（每行一个）", key="input_fielder_promoter")

        gene_ids = read_gene_ids(uploaded_file, manual_input)
        gene_ids = [g.strip() for g in gene_ids if g.strip()]

        if not gene_ids:
            st.info("请上传文件或输入 Fielder 基因号")
            st.stop()

        st.info("当前数据库内置的是 promoter_2000，因此这里固定返回 ATG 上游 2000 bp 启动子序列。")
        st.info(f"待处理基因数：{len(gene_ids)}")

        if st.button("开始抓取 Fielder 启动子", key="btn_fielder_promoter"):
            progress = st.progress(0)
            status_text = st.empty()

            fasta_records = []
            summary_rows = []
            failed_genes = []

            for idx, input_gene_id in enumerate(gene_ids, 1):
                status_text.text(f"正在查询: {input_gene_id} ({idx}/{len(gene_ids)})")

                promoter_df = get_fielder_promoter(input_gene_id)

                if promoter_df.empty:
                    failed_genes.append(input_gene_id)
                    summary_rows.append([input_gene_id, "", "", "", "", "", 0])
                    progress.progress(idx / len(gene_ids))
                    continue

                row = promoter_df.iloc[0]

                primary_id = row.get("primary_gene_id", "")
                chrom = row.get("chromosome", "")
                strand = row.get("strand", "")
                promoter_start = row.get("promoter_start", "")
                promoter_end = row.get("promoter_end", "")
                promoter_length = row.get("promoter_length", "")
                promoter_seq = row.get("promoter_sequence", "")

                header = f">{input_gene_id}|{primary_id}|promoter_2000|{chrom}:{promoter_start}-{promoter_end}|{strand}"
                fasta_records.append(f"{header}\n{promoter_seq}\n")

                summary_rows.append([
                    input_gene_id,
                    primary_id if pd.notna(primary_id) else "",
                    chrom if pd.notna(chrom) else "",
                    promoter_start if pd.notna(promoter_start) else "",
                    promoter_end if pd.notna(promoter_end) else "",
                    strand if pd.notna(strand) else "",
                    promoter_length if pd.notna(promoter_length) else 0
                ])

                progress.progress(idx / len(gene_ids))

            summary_df = pd.DataFrame(
                summary_rows,
                columns=["输入基因号", "Fielder主键", "染色体", "启动子起点", "启动子终点", "链方向", "启动子长度"]
            )

            st.success("🎉 Fielder 启动子抓取完成！")
            st.dataframe(summary_df, use_container_width=True)

            st.download_button(
                "📥 下载 Fielder 启动子 FASTA",
                "".join(fasta_records),
                "fielder_promoter_sequences.fasta",
                "text/plain"
            )

            st.download_button(
                "📥 下载 Fielder 启动子统计 CSV",
                summary_df.to_csv(index=False).encode("utf-8-sig"),
                "fielder_promoter_summary.csv",
                "text/csv"
            )

            if failed_genes:
                st.warning(f"⚠️ 以下基因未获取到启动子序列: {', '.join(failed_genes)}")

            status_text.empty()
            progress.empty()

    # -------------------- 功能 5 --------------------
    elif tool == "中国春启动子抓取":
        st.header("🌱 中国春基因启动子抓取（本地数据库版）")

        uploaded_file = st.file_uploader("上传基因列表 TXT", type=["txt"], key="file_cs_promoter")
        manual_input = st.text_area("或者手动输入基因号（每行一个）", key="input_cs_promoter")

        gene_ids = read_gene_ids(uploaded_file, manual_input)
        gene_ids = [g.strip() for g in gene_ids if g.strip()]

        if not gene_ids:
            st.info("请上传文件或输入基因号")
            st.stop()

        st.info("当前数据库内置的是 promoter_2000，因此这里固定返回 ATG 上游 2000 bp 启动子序列。")
        st.info(f"待处理基因数: {len(gene_ids)}")

        if st.button("开始抓取", key="btn_cs_promoter"):
            progress = st.progress(0)
            status_text = st.empty()

            fasta_records = []
            summary_rows = []
            failed_genes = []

            for idx, input_gene_id in enumerate(gene_ids, 1):
                status_text.text(f"正在查询: {input_gene_id} ({idx}/{len(gene_ids)})")

                promoter_df = get_promoter(input_gene_id)

                if promoter_df.empty:
                    failed_genes.append(input_gene_id)
                    summary_rows.append([input_gene_id, "", "", "", "", "", 0])
                    progress.progress(idx / len(gene_ids))
                    continue

                row = promoter_df.iloc[0]

                primary_id = row.get("primary_gene_id", "")
                chrom = row.get("chromosome", "")
                strand = row.get("strand", "")
                promoter_start = row.get("promoter_start", "")
                promoter_end = row.get("promoter_end", "")
                promoter_length = row.get("promoter_length", "")
                promoter_seq = row.get("promoter_sequence", "")

                header = f">{input_gene_id}|{primary_id}|promoter_2000|{chrom}:{promoter_start}-{promoter_end}|{strand}"
                fasta_records.append(f"{header}\n{promoter_seq}\n")

                summary_rows.append([
                    input_gene_id,
                    primary_id if pd.notna(primary_id) else "",
                    chrom if pd.notna(chrom) else "",
                    promoter_start if pd.notna(promoter_start) else "",
                    promoter_end if pd.notna(promoter_end) else "",
                    strand if pd.notna(strand) else "",
                    promoter_length if pd.notna(promoter_length) else 0
                ])

                progress.progress(idx / len(gene_ids))

            summary_df = pd.DataFrame(
                summary_rows,
                columns=["输入基因号", "统一主键", "染色体", "启动子起点", "启动子终点", "链方向", "启动子长度"]
            )

            st.success("🎉 启动子抓取完成！")
            st.dataframe(summary_df, use_container_width=True)

            st.download_button(
                "📥 下载启动子 FASTA",
                "".join(fasta_records),
                "cs_promoter_sequences.fasta",
                "text/plain"
            )

            st.download_button(
                "📥 下载启动子统计 CSV",
                summary_df.to_csv(index=False).encode("utf-8-sig"),
                "cs_promoter_summary.csv",
                "text/csv"
            )

            if failed_genes:
                st.warning(f"⚠️ 以下基因未获取到启动子序列: {', '.join(failed_genes)}")

            status_text.empty()
            progress.empty()

    # -------------------- 功能 6 --------------------
    elif tool == "GO富集分析":
        st.header("📊 GO 富集分析")
        st.caption("输入 DEG 列表（一行一个基因号），输出 GO 富集结果表和条形图。")

        uploaded_file = st.file_uploader("上传 DEG TXT 文件（一行一个基因号）", type=["txt"], key="file_go_deg")
        manual_input = st.text_area("或者手动输入 DEG（每行一个）", height=200, key="input_go_deg")

        gene_ids = read_gene_ids(uploaded_file, manual_input)
        gene_ids = [g.strip() for g in gene_ids if g.strip()]

        if not gene_ids:
            st.info("请上传 DEG 文件或手动输入基因号")
            st.stop()

        col1, col2, col3 = st.columns(3)
        with col1:
            top_n = st.number_input("每个大类展示前 N 个 term", min_value=3, max_value=30, value=10, step=1)
        with col2:
            padj_cutoff = st.number_input("FDR 阈值", min_value=0.0001, max_value=1.0, value=0.05, step=0.01, format="%.4f")
        with col3:
            plot_metric = st.selectbox(
                "作图横轴",
                options=["fdr", "pvalue"],
                format_func=lambda x: "-log10(FDR)" if x == "fdr" else "-log10(P-value)"
            )

        min_size = st.number_input("最小 GO 基因集大小", min_value=1, max_value=50, value=3, step=1)
        max_size = st.number_input("最大 GO 基因集大小", min_value=10, max_value=10000, value=2000, step=10)

        BASE_DIR = os.path.dirname(os.path.abspath(__file__))
        mapping_dir = os.path.join(BASE_DIR, "data", "go_mapping")

        term2gene_path = os.path.join(mapping_dir, "TERM2GENE_protein_coding.tsv")
        term2name_path = os.path.join(mapping_dir, "TERM2NAME_protein_coding.tsv")
        metadata_path = os.path.join(mapping_dir, "wheat_go_metadata.tsv")
        background_path = os.path.join(mapping_dir, "wheat_protein_coding_genes.tsv")

        required_files = [term2gene_path, term2name_path, metadata_path, background_path]
        missing_files = [fp for fp in required_files if not os.path.exists(fp)]

        if missing_files:
            st.error("以下 GO 注释文件不存在，请检查 data/go_mapping/ 目录：")
            for fp in missing_files:
                st.code(fp)
            st.stop()

        st.info(f"待分析基因数: {len(gene_ids)}")

        if st.button("开始 GO 富集分析", key="btn_go_enrichment"):
            with st.spinner("正在进行 GO 富集分析，请稍候..."):
                try:
                    results_df, sig_df, summary_df = run_go_enrichment(
                        gene_list=gene_ids,
                        term2gene_path=term2gene_path,
                        term2name_path=term2name_path,
                        metadata_path=metadata_path,
                        background_path=background_path,
                        min_size=min_size,
                        max_size=max_size,
                        padj_cutoff=padj_cutoff
                    )

                    if results_df.empty:
                        st.warning("没有得到任何 GO 富集结果，请检查输入基因 ID 是否与背景一致。")
                        st.stop()

                    st.success("✅ GO 富集分析完成！")

                    st.subheader("分析摘要")
                    st.dataframe(summary_df, use_container_width=True)

                    st.subheader("显著富集结果")
                    if sig_df.empty:
                        st.warning("当前 FDR 阈值下没有显著富集条目，下面展示全部结果。")
                        st.dataframe(results_df, use_container_width=True)
                    else:
                        st.dataframe(sig_df, use_container_width=True)

                    plot_df = sig_df if not sig_df.empty else results_df
                    plot_bytes = create_go_barplot_bytes(
                        df=plot_df,
                        top_n=top_n,
                        plot_metric=plot_metric
                    )

                    if plot_bytes is not None:
                        st.subheader("GO 富集条形图")
                        st.image(plot_bytes, caption="GO enrichment barplot", use_container_width=True)
                        st.download_button(
                            "📥 下载富集图 PNG",
                            data=plot_bytes,
                            file_name="GO_enrichment_barplot.png",
                            mime="image/png"
                        )

                    st.download_button(
                        "📥 下载全部结果 TSV",
                        data=results_df.to_csv(sep="\t", index=False).encode("utf-8-sig"),
                        file_name="GO_enrichment_results.tsv",
                        mime="text/tab-separated-values"
                    )

                    st.download_button(
                        "📥 下载显著结果 TSV",
                        data=sig_df.to_csv(sep="\t", index=False).encode("utf-8-sig"),
                        file_name="GO_enrichment_results_sig.tsv",
                        mime="text/tab-separated-values"
                    )

                    st.download_button(
                        "📥 下载分析摘要 TSV",
                        data=summary_df.to_csv(sep="\t", index=False).encode("utf-8-sig"),
                        file_name="GO_enrichment_summary.tsv",
                        mime="text/tab-separated-values"
                    )

                except Exception as e:
                    st.error(f"GO 富集分析失败：{e}")

    # -------------------- 功能 7 --------------------
    elif tool == "KEGG富集分析":
        st.header("📊 KEGG 富集分析")
        st.caption("输入 DEG 列表（一行一个基因号），基于本地 gene-KO 和 KEGG KO-pathway 映射进行 KEGG 富集分析。")

        uploaded_file = st.file_uploader(
            "上传 DEG TXT 文件（一行一个基因号）",
            type=["txt"],
            key="file_kegg_deg"
        )

        manual_input = st.text_area(
            "或者手动输入 DEG（每行一个）",
            height=200,
            key="input_kegg_deg"
        )

        gene_ids = read_gene_ids(uploaded_file, manual_input)
        gene_ids = [g.strip() for g in gene_ids if g.strip()]

        if not gene_ids:
            st.info("请上传 DEG 文件或手动输入基因号")
            st.stop()

        st.info(f"待分析基因数: {len(gene_ids)}")

        col1, col2, col3 = st.columns(3)

        with col1:
            top_n = st.number_input(
                "图中展示前 N 个通路",
                min_value=3,
                max_value=50,
                value=20,
                step=1
            )

        with col2:
            pvalue_cutoff = st.number_input(
                "P-value 阈值",
                min_value=0.0001,
                max_value=1.0,
                value=0.05,
                step=0.01,
                format="%.4f"
            )

        with col3:
            use_sig_only = st.checkbox(
                "优先展示显著通路",
                value=True,
                help="如果存在 pvalue 小于阈值的通路，则优先使用显著通路绘图；否则展示全部结果中 pvalue 最小的通路。"
            )

        col4, col5 = st.columns(2)

        with col4:
            min_size = st.number_input(
                "最小 pathway KO 数",
                min_value=1,
                max_value=100,
                value=3,
                step=1
            )

        with col5:
            max_size = st.number_input(
                "最大 pathway KO 数",
                min_value=10,
                max_value=5000,
                value=500,
                step=10
            )

        with st.expander("高级绘图参数"):
            col6, col7, col8 = st.columns(3)

            with col6:
                clip_minus_log10_p = st.checkbox(
                    "截断极端 -log10(P-value)",
                    value=True,
                    help="建议开启，避免极端小 pvalue 把横轴拉得过长。"
                )

            with col7:
                clip_quantile = st.number_input(
                    "截断分位数",
                    min_value=0.50,
                    max_value=1.00,
                    value=0.95,
                    step=0.01,
                    format="%.2f"
                )

            with col8:
                label_wrap_width = st.number_input(
                    "通路名称换行宽度",
                    min_value=20,
                    max_value=80,
                    value=42,
                    step=2
                )

            col9, col10 = st.columns(2)

            with col9:
                bubble_min_size = st.number_input(
                    "最小气泡大小",
                    min_value=5,
                    max_value=200,
                    value=25,
                    step=5
                )

            with col10:
                bubble_max_size = st.number_input(
                    "最大气泡大小",
                    min_value=50,
                    max_value=800,
                    value=220,
                    step=10
                )

        BASE_DIR = os.path.dirname(os.path.abspath(__file__))
        mapping_dir = os.path.join(BASE_DIR, "data", "kegg_mapping")

        gene2ko_path = os.path.join(mapping_dir, "gene2ko_clean.tsv")
        ko2pathway_path = os.path.join(mapping_dir, "kegg_ko2pathway.tsv")
        pathway2name_path = os.path.join(mapping_dir, "kegg_pathway2name.tsv")

        required_files = [
            gene2ko_path,
            ko2pathway_path,
            pathway2name_path
        ]

        missing_files = [fp for fp in required_files if not os.path.exists(fp)]

        if missing_files:
            st.error("以下 KEGG 注释文件不存在，请检查 data/kegg_mapping/ 目录：")
            for fp in missing_files:
                st.code(fp)
            st.stop()

        with st.expander("当前使用的 KEGG 本地注释文件"):
            st.code(gene2ko_path)
            st.code(ko2pathway_path)
            st.code(pathway2name_path)

        if st.button("开始 KEGG 富集分析", key="btn_kegg_enrichment"):
            with st.spinner("正在进行 KEGG 富集分析，请稍候..."):
                try:
                    results_df, sig_df, summary_df = run_kegg_enrichment(
                        gene_list=gene_ids,
                        gene2ko_path=gene2ko_path,
                        ko2pathway_path=ko2pathway_path,
                        pathway2name_path=pathway2name_path,
                        min_size=min_size,
                        max_size=max_size,
                        pvalue_cutoff=pvalue_cutoff
                    )

                    if results_df.empty:
                        st.warning("没有得到任何 KEGG 富集结果，请检查输入基因 ID 是否与 gene2ko_clean.tsv 中的基因 ID 一致。")
                        st.subheader("分析摘要")
                        st.dataframe(summary_df, use_container_width=True)
                        st.stop()

                    st.success("✅ KEGG 富集分析完成！")

                    st.subheader("分析摘要")
                    st.dataframe(summary_df, use_container_width=True)

                    st.subheader("显著富集结果")
                    if sig_df.empty:
                        st.warning("当前 P-value 阈值下没有显著富集通路，下面展示全部结果。")
                        st.dataframe(results_df, use_container_width=True)
                    else:
                        st.dataframe(sig_df, use_container_width=True)

                    if use_sig_only and not sig_df.empty:
                        plot_df = sig_df
                    else:
                        plot_df = results_df

                    bubble_bytes = create_kegg_bubbleplot_bytes(
                        df=plot_df,
                        top_n=top_n,
                        bubble_min_size=bubble_min_size,
                        bubble_max_size=bubble_max_size,
                        clip_minus_log10_p=clip_minus_log10_p,
                        clip_mode="quantile",
                        clip_quantile=clip_quantile,
                        clip_fixed_value=30,
                        min_x_cap=5,
                        label_wrap_width=label_wrap_width
                    )

                    barplot_bytes = create_kegg_barplot_bytes(
                        df=plot_df,
                        top_n=top_n,
                        clip_minus_log10_p=clip_minus_log10_p,
                        clip_mode="quantile",
                        clip_quantile=clip_quantile,
                        clip_fixed_value=30,
                        min_x_cap=5,
                        label_wrap_width=label_wrap_width
                    )

                    if bubble_bytes is not None:
                        st.subheader("KEGG 富集气泡图")
                        st.image(
                            bubble_bytes,
                            caption="KEGG enrichment bubble plot",
                            use_container_width=True
                        )

                        st.download_button(
                            "📥 下载 KEGG 气泡图 PNG",
                            data=bubble_bytes,
                            file_name="KEGG_enrichment_bubbleplot.png",
                            mime="image/png"
                        )

                    if barplot_bytes is not None:
                        st.subheader("KEGG 富集条形图")
                        st.image(
                            barplot_bytes,
                            caption="KEGG enrichment barplot",
                            use_container_width=True
                        )

                        st.download_button(
                            "📥 下载 KEGG 条形图 PNG",
                            data=barplot_bytes,
                            file_name="KEGG_enrichment_barplot.png",
                            mime="image/png"
                        )

                    st.download_button(
                        "📥 下载全部 KEGG 富集结果 TSV",
                        data=results_df.to_csv(sep="\t", index=False).encode("utf-8-sig"),
                        file_name="KEGG_enrichment_results_all.tsv",
                        mime="text/tab-separated-values"
                    )

                    st.download_button(
                        "📥 下载显著 KEGG 富集结果 TSV",
                        data=sig_df.to_csv(sep="\t", index=False).encode("utf-8-sig"),
                        file_name="KEGG_enrichment_results_sig.tsv",
                        mime="text/tab-separated-values"
                    )

                    st.download_button(
                        "📥 下载 KEGG 分析摘要 TSV",
                        data=summary_df.to_csv(sep="\t", index=False).encode("utf-8-sig"),
                        file_name="KEGG_enrichment_summary.tsv",
                        mime="text/tab-separated-values"
                    )

                except Exception as e:
                    st.error(f"KEGG 富集分析失败：{e}")

if __name__ == "__main__":
    main()
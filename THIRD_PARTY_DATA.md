# Third-party data sources, attribution, and terms of use

## Scope

WheatGeneToolkit contains original software as well as third-party biological data and
locally generated derivative datasets. The Apache License 2.0 in `LICENSE` applies to
the original software code and project documentation only. It does **not** replace or
override the licenses, database rights, citation requirements, or terms of use of the
data providers listed below.

Users who download, redistribute, deploy, or commercially use a copy of this repository
are responsible for checking the current terms of each provider. Attribution is not a
substitute for permission where a provider requires a separate licence.

## 数据来源与使用说明

WheatGeneToolkit 同时包含原创程序、第三方生物学数据以及由第三方数据生成的本地衍生数据。
根目录 `LICENSE` 中的 Apache License 2.0 **仅适用于本项目原创代码与项目文档**，不改变、
替代或覆盖下列数据提供方的许可证、数据库权利、引用要求或使用条款。

下载、再分发、公开部署或商业使用本仓库前，使用者应自行核对各数据提供方的最新条款。
在数据提供方要求另行授权时，仅注明来源并不等同于获得授权。

## 1. Chinese Spring reference genome / 中国春参考基因组

- **Source:** International Wheat Genome Sequencing Consortium (**IWGSC**),
  Chinese Spring `IWGSC RefSeq v1.1` gene annotation and its associated reference
  assembly.
- **Used for:** Chinese Spring gene models, gene/transcript/CDS/protein sequences,
  genomic structure and locally extracted promoter sequences.
- **Relevant repository content:** Chinese Spring-related records under `data/db/`,
  including core annotation, sequence, structure and promoter database shards.
- **Primary citation:** [Reference 1](REFERENCES.md#ref-1), International Wheat Genome Sequencing Consortium (IWGSC).
  *Shifting the limits in wheat research and breeding using a fully annotated
  reference genome.* Science 361, eaar7191 (2018).
  https://doi.org/10.1126/science.aar7191
- **Resource:** https://urgi.versailles.inrae.fr/download/iwgsc/

The upstream IWGSC/URGI terms and citation requirements continue to apply. The local
SQLite representation and promoter extraction performed by this project do not grant
new rights over the underlying genome or annotation data.

> 中文说明：本项目所称“中国春 IWGSC v1.1”主要指 IWGSC RefSeq v1.1 基因注释及其配套
> 参考序列。由这些序列提取的启动子和转换生成的 SQLite 数据库仍属于衍生数据，不能仅因
> 格式转换而被重新授权为 Apache-2.0。

## 2. Fielder reference genome / Fielder 参考基因组

- **Source:** chromosome-scale assembly of *Triticum aestivum* cultivar Fielder,
  `wheat_cv_fielder_v1_assembly`.
- **Assembly accession:** NCBI Assembly `GCA_907166925.1`.
- **Used for:** Fielder gene records and locally extracted 2,000 bp promoter sequences.
- **Relevant repository content:** Fielder-related records under `data/db/`, especially
  `fielder_gene_core.db` and `fielder_promoter_sequence/`.
- **Primary citation:** [Reference 2](REFERENCES.md#ref-2), Sato K. et al. *Chromosome-scale genome assembly of the
  transformation-amenable common wheat cultivar 'Fielder'.* DNA Research 28(3),
  dsab008 (2021). https://doi.org/10.1093/dnares/dsab008
- **NCBI resource:** https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_907166925.1/

The terms attached to the source assembly, annotation downloads and their hosting
repositories continue to apply to the corresponding local and derived data.

## 3. Homology relationships / 同源基因关系

- **Source:** Triticeae-GeneTribe (**TGT**), China Agricultural University.
- **Used for:** Chinese Spring self-homology and Chinese Spring-to-Fielder homology
  lookup tables.
- **Relevant repository content:** `data/db/homolog/`.
- **Primary citation:** [Reference 3](REFERENCES.md#ref-3), Chen Y. et al. *A Collinearity-incorporating Homology Inference
  Strategy for Connecting Emerging Assemblies in Triticeae Tribe as a Pilot Practice
  in the Plant Pangenomic Era.* Molecular Plant 13, 1694-1708 (2020).
  https://doi.org/10.1016/j.molp.2020.09.019
- **Resource:** https://wheat.cau.edu.cn/TGT/

TGT requests citation of the above publication. No blanket open-data licence for a
redistributed copy of the TGT database has been identified in the resource page.
Consequently, TGT-derived tables are excluded from this project's Apache-2.0 grant and
remain subject to TGT's current terms and any permission required by its maintainers.

## 4. Gene annotation and gene-ID conversion / 基因功能注释与基因号转换

- **Source:** **WheatOmics 1.0**, Shandong Agricultural University.
- **Used for:** gene functional descriptions, alias mapping and conversion among wheat
  gene-ID generations.
- **Relevant repository content:** annotation and alias records under `data/db/core/`.
- **Primary citation:** [Reference 4](REFERENCES.md#ref-4), Ma S. et al. *WheatOmics: A platform combining multiple omics
  data to accelerate functional genomics studies in wheat.* Molecular Plant 14,
  1965-1968 (2021). https://doi.org/10.1016/j.molp.2021.10.006
- **Resource:** https://wheatomics.sdau.edu.cn/

WheatOmics requests citation of the above publication. No explicit blanket licence
allowing this project to relicense the underlying database under Apache-2.0 has been
identified. WheatOmics-derived records are therefore attributed but excluded from the
project software licence.

## 5. JASPAR motif data / JASPAR 启动子 motif 数据

- **Source:** **JASPAR 2026 CORE Plants**, non-redundant position frequency matrices
  and metadata.
- **Used for:** local PFM-to-PWM conversion and promoter motif scanning.
- **Relevant repository content:** `data/motif_db/jaspar_plants/`.
- **Resource and licence:** https://jaspar.elixir.no/about/
- **Primary citation:** [Reference 5](REFERENCES.md#ref-5), JASPAR 2026 database paper.

JASPAR is licensed under the **Creative Commons Attribution 4.0 International
Licence (CC BY 4.0)**. JASPAR-derived PFM/PWM, metadata and threshold files retain that
attribution requirement. Cite the JASPAR release/publication appropriate to the version
used and identify modifications such as PFM-to-PWM conversion and locally computed
background thresholds.

Suggested attribution:

> Contains data from JASPAR 2026 CORE Plants, licensed under CC BY 4.0. PWM and
> background-threshold files were generated locally by WheatGeneToolkit.

## 6. Gene Ontology mappings / GO 注释映射

- **Source:** **Ensembl Plants BioMart**, wheat gene dataset, including Gene Ontology
  associations; underlying GO terms and ontology content originate from the Gene
  Ontology Consortium.
- **Used for:** local wheat gene-to-GO mappings, GO metadata and enrichment background.
- **Relevant repository content:** `data/go_mapping/`.
- **Ensembl Plants:** https://plants.ensembl.org/biomart/martview
- **Ensembl data export policy:** https://plants.ensembl.org/info/data/export.html
- **GO licence and citation policy:** https://geneontology.org/docs/go-citation-policy/
- **Primary citations:** [References 6-9](REFERENCES.md#ref-6), covering Ensembl,
  BioMart and the Gene Ontology knowledgebase.

Ensembl describes exported data as open access and reusable, while Gene Ontology data
products are licensed under **CC BY 4.0**. Redistributed GO-derived files must preserve
attribution. For reproducibility, the exact Ensembl Plants release, BioMart dataset,
query date, selected attributes, GO release date and (where available) release DOI
should be recorded. These release details are currently not recoverable from the files
alone and should be added when known.

## 7. KEGG and BlastKOALA data / KEGG 与 BlastKOALA 数据

- **Source:** **KEGG BlastKOALA**, used for KO assignment; KEGG KO-to-pathway and
  pathway-name information used for local enrichment analysis.
- **Relevant repository content:** `data/kegg_mapping/`.
- **BlastKOALA:** https://www.kegg.jp/blastkoala/
- **KEGG copyright and licensing terms:** https://www.kegg.jp/kegg/legal.html
- **Primary citations:** [References 10-11](REFERENCES.md#ref-10), covering BlastKOALA
  and the KEGG database.

KEGG is a copyrighted database product and expressly states that it is not a public
database. Academic use of the website may be free, but academic users providing a
service based on KEGG are requested to obtain an academic service-provider licence;
non-academic use requires a commercial licence. Therefore:

1. KEGG-derived mapping tables are **not** licensed under Apache-2.0 or CC BY 4.0 by
   this project.
2. Inclusion in this repository does not grant recipients a right to redistribute or
   commercially use KEGG-derived data.
3. Public deployment and redistribution should be reviewed against the current KEGG
   terms. Obtain the appropriate KEGG licence or replace/remove locally redistributed
   KEGG-derived mapping files if the intended use is not covered.
4. Cite BlastKOALA and the KEGG publications requested by KEGG for the services and
   data used.

## Locally generated outputs / 本地生成内容

The project's original parsing, statistical analysis, plotting and web-interface code
is licensed under Apache-2.0. Purely original outputs that do not reproduce protected
third-party data may be used under the user's chosen terms. Outputs containing gene
sequences, annotations, homology tables, motif matrices, GO mappings or KEGG mappings
remain subject to the relevant upstream terms described above.

本项目原创的解析、统计、绘图和网页程序采用 Apache-2.0。未复制第三方受保护数据的原创
分析结果可由使用者自行决定使用方式；但包含基因序列、注释、同源关系、motif 矩阵、GO
映射或 KEGG 映射的输出仍须遵守相应上游条款。

## No endorsement or warranty / 不代表背书且不提供担保

References to IWGSC, URGI, NCBI, TGT, WheatOmics, JASPAR, Ensembl, the Gene Ontology
Consortium or KEGG are for attribution only and do not imply endorsement of
WheatGeneToolkit. Third-party data are provided without additional warranty by this
project.

## Complete bibliography / 完整参考文献

See [REFERENCES.md](REFERENCES.md) for the complete formatted bibliography and a
module-to-reference mapping table.

完整格式化参考文献及“项目模块—参考文献”对应表见 [REFERENCES.md](REFERENCES.md)。

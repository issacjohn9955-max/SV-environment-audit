# 鸡结构变异与环境关联 ｜ Chicken SV–Environment Association

[中文](#中文) · [English](#english)

## 中文

### 1. 仓库用途与数据状态

本仓库用于鸡结构变异（structural variation，SV）的基因型—环境关联分析，也用于在建模前审核分型个体与环境表型能否逐一对应。流程从**实测 SV 分型**出发，结合已审核的个体名单和气候坐标组，运行偏秩相关、偏冗余分析（partial RDA）和 LFMM2。三种方法均可比较由同一套 SV PCA 得到的 PC1–3、PC1–5 和 PC1–10。原始 VCF/BCF、源登记工作簿及 WorldClim 栅格保存在服务器。

| 数据版本 | 规模与用途 |
| --- | --- |
| 2026-11-23 较早登记 | 1,726 个唯一 ID，是原产地证据报告和新旧登记比对的基线 |
| 2026-11-23 最新分型库存 | 1,732 个唯一 ID，包含上一版全部个体和 6 个新增个体。共有个体中有 21 条原始品种标签不同，需区分名称写法与真实身份差异 |
| 2026-11-18 历史主队列 | 565 个体 BCF 中的 405 个历史环境个体，25 个环境坐标组。保留实际环境输入和三种方法的历史诊断，供追溯 |
| 模拟验证队列 | 112 个模拟个体、28 个坐标组、160 个输入等位列。三种方法与三个 PC 规格共九套组合已试跑 |

**当前审核进度。** 最新清单的 1,732 个主 VCF 样本头均与登记的 VCF 样本 ID 精确匹配。66 条 A/B 级品种代表点记录完整提取了五个气候变量。按现有证据，1,189 个体有待按批次确认的候选环境，405 个体暂无完整候选环境，112 个体属于排除的商业配套系，10 个体因身份冲突暂挂，16 个体因文库类型待核暂挂。VCF 头核对证明文件样本列与登记一致，个体真实身份、GT 质量和环境归属仍以相应审核证据为准。 [审核统计](04_genotype_environment_audit/review/audit_summary.json)

### 2. 仓库文件索引

| 位置 | 文件与用途 |
| --- | --- |
| 仓库根目录 | `README.md` 为本说明，[LICENSE](LICENSE) 为专有许可，`.gitignore` 忽略本地运行输出和原始基因型文件，`.gitattributes` 标明 PDF 二进制属性和 TSV 空字段的差异检查规则 |
| [01_partial_correlation](01_partial_correlation/) | [README](01_partial_correlation/README.md)、[run.py](01_partial_correlation/run.py) 分析入口、[partial_rank.py](01_partial_correlation/partial_rank.py) 偏秩相关计算内核、[historical_main405_summary.tsv](01_partial_correlation/historical_main405_summary.tsv) 历史可检验与显著位点计数 |
| [02_partial_RDA](02_partial_RDA/) | [README](02_partial_RDA/README.md)、[run.R](02_partial_RDA/run.R) 建模和置换入口、[rda_metrics_vegan.tsv](02_partial_RDA/rda_metrics_vegan.tsv) 历史整体结果、[rda_axis_tests.tsv](02_partial_RDA/rda_axis_tests.tsv) 历史分轴检验 |
| [03_LFMM2](03_LFMM2/) | [README](03_LFMM2/README.md)、[run.R](03_LFMM2/run.R) 建模与校准入口、历史 [K5](03_LFMM2/historical_main405_K5_diagnostics.tsv)、[K10](03_LFMM2/historical_main405_K10_diagnostics.tsv)、[K15](03_LFMM2/historical_main405_K15_diagnostics.tsv) 逐环境校准与发现数 |
| [04_genotype_environment_audit](04_genotype_environment_audit/) | [审核指南](04_genotype_environment_audit/README.md) 和 [audit.py](04_genotype_environment_audit/audit.py) 审核重建入口 |
| [environment/](04_genotype_environment_audit/environment/) | [origin_evidence_1726_snapshot.tsv](04_genotype_environment_audit/environment/origin_evidence_1726_snapshot.tsv) 原产地和证据等级，[candidate_climate_by_breed.tsv](04_genotype_environment_audit/environment/candidate_climate_by_breed.tsv) 候选代表点气候，[historical_main405_group_climate.csv](04_genotype_environment_audit/environment/historical_main405_group_climate.csv) 历史组层原值，[historical_main405_individual_standardized_climate.csv](04_genotype_environment_audit/environment/historical_main405_individual_standardized_climate.csv) 历史 LFMM2 个体层标准化值 |
| [review/](04_genotype_environment_audit/review/) | [audit_summary.json](04_genotype_environment_audit/review/audit_summary.json) 计数、状态与输入 SHA-256，[individual_alignment.tsv](04_genotype_environment_audit/review/individual_alignment.tsv) 全量个体对照，[review_queue.tsv](04_genotype_environment_audit/review/review_queue.tsv) 待核问题子集，[breed_label_crosswalk.tsv](04_genotype_environment_audit/review/breed_label_crosswalk.tsv) 有据可查的品种别名，[historical_main405_alignment.tsv](04_genotype_environment_audit/review/historical_main405_alignment.tsv) 历史样本头与环境匹配，[原产地审核 Markdown](04_genotype_environment_audit/review/origin_review_1726_20261123.md) 与 [原产地审核 PDF](04_genotype_environment_audit/review/鸡群分型登记对照与品种原产地核查_20261123.pdf) |
| [workflow/](04_genotype_environment_audit/workflow/) | [prepare.py](04_genotype_environment_audit/workflow/prepare.py) 共享预处理，[run.template.json](04_genotype_environment_audit/workflow/run.template.json) 正式参数模板，[PC3/5/10 原理和操作 PDF](04_genotype_environment_audit/workflow/SV_PC_LFMM2_RDA_操作与原理.pdf) |
| [example/](04_genotype_environment_audit/example/) | [run.json](04_genotype_environment_audit/example/run.json) 试跑配置，[dosage.tsv](04_genotype_environment_audit/example/dosage.tsv) 模拟 SV 剂量，[samples.tsv](04_genotype_environment_audit/example/samples.tsv) 模拟个体名单，[groups.tsv](04_genotype_environment_audit/example/groups.tsv) 模拟环境组与置换区块 |

### 3. 如何审核个体与环境

1. 从 [待审核清单](04_genotype_environment_audit/review/review_queue.tsv) 查找新增、原始标签变化、VCF 头问题、环境缺口及暂挂记录。再用 [全量个体对照表](04_genotype_environment_audit/review/individual_alignment.tsv) 查看同一批次的其他个体。具有候选环境且没有额外问题的旧个体也在全量表，仍需审查。
2. 按精确 ID 核对登记 ID、登记的 VCF 样本 ID、实际 VCF 列名和来源批次。逐字标签不同是复核线索。使用名称对照表和原产地审核报告判断别名、重复、身份冲突或文库暂挂。XB、MH 等已归为商业配套系，现阶段的环境队列按审核记录排除。
3. 检查坐标用途和证据。A 为专门整理的产区代表点，B 为既有审计沿用点，C 为待落点地名，D 为分布宽或需按地区、亚种及批次处理，E 为商业配套系，U 为身份冲突。代表点气候适用性应按样本来源批次核定。缺失字段以空白保存。
4. 审核通过后，冻结 `samples.tsv` 和 `groups.tsv`。在工作记录中保存个体裁定、理由、证据、环境版本和建议的地理置换区块。历史 405 队列的两张气候表按历史口径保留，当前队列按自己的审核结果建表。

用于关联的候选气候由 **WorldClim v2.1** 的 **1970–2000 年、5 弧分**栅格按经纬度所在像元提取。BIO3 为等温性，BIO15 为降水季节性，BIO18 为最暖季度降水量（mm），`wind` 是 12 个月风速均值（m/s），`srad` 是 12 个月太阳辐射均值换算为 W/m²。品种原产区的气候代表值用于地理环境关联。 [WorldClim 数据说明](https://www.worldclim.org/data/worldclim21.html)

### 4. 方法与参数

| 方法 | 数据单位、模型和主要结果 |
| --- | --- |
| 偏秩相关 | 每个坐标组一行，响应为组 SV 等位频率。对单一环境变量和频率分别排秩，移除 SV PC 与可选经纬度后计算残差相关。按环境和 PC 规格对可检验位点的 P 值做 BH。可用 `--frequency` 传入包含部分缺失的频率表，行键和顺序须与组设计表一致 |
| 偏 RDA | 每个坐标组一行，响应为各组均有可估频率且有变异的 SV 列。联合拟合环境变量，条件化 SV PC 与可选经纬度。用 vegan 的 reduced-model 残差置换给出整体、环境边际和轴检验。整体 P 和轴 BH q 均低于阈值后，再按标准化载荷绝对值筛选探索性位点 |
| LFMM2 | 每个个体一行，使用全队列均有实测且有变异的 SV 剂量列。联合纳入五个环境、显式 SV PC、可选经纬度和固定维数的潜因子。基于 LEA 拟合潜结构并给出分环境系数、名义 P、基因组控制 P 和按环境／PC 规格的 BH q |

VCF/BCF 输入按 ALT 拆列。逐个体剂量为 `2 × ALT 拷贝数 / 已调用倍性`，组频率为 `已调用 ALT 拷贝数之和 / 已调用倍性之和`。从同一实测矩阵拟合一次 SV PCA，各组 PC 为组内个体 PC 的平均。只有 PCA 标准化工作矩阵使用观测均值投影处理缺失。关联矩阵由实测 GT 构建，全流程无 SNP 输入。

正式模板默认 `lfmm_k=5`、`permutations=999`、`min_called_per_group=3`、`adjust_geography=true`、`rda_blocked=false`、`alpha=0.05`、`loading_z=3`、`seed=20261123`。PC3、PC5、PC10 可比较显式结构校正的敏感性，LFMM2 的 K 在比较时固定。正式偏 RDA 置换应根据环境组的交换性选择自由或组内区块置换。模拟配置采用 4 个区块及 39 次置换。

结果解读应同时考虑环境共线性、PC10 剩余自由度、LFMM2 基因组控制系数、组内共享气候及样本来源证据。RDA 载荷尾部是探索性证据。偏相关和 LFMM2 的 q 值各有自己的检验家族。 [LFMM2 方法论文](https://doi.org/10.1093/molbev/msz008) · [RDA 环境关联研究](https://doi.org/10.1111/mec.14584) · [vegan 置换文档](https://vegandevs.github.io/vegan/reference/anova.cca.html)

### 5. 安装与从示例运行

使用 Linux、Python 3.11+、NumPy、pandas、SciPy ≥1.11，以及 R 里的 LEA、vegan、permute、jsonlite。审核重建还需要 Pillow、openpyxl 和 bcftools。正式 VCF/BCF 输入也需要 bcftools。LEA 可通过 Bioconductor 安装。可用 `SV_R_LIBRARY` 指定 R 附加包目录。

```bash
python3 -m pip install numpy pandas 'scipy>=1.11' Pillow openpyxl
Rscript -e 'install.packages(c("BiocManager", "vegan", "permute", "jsonlite"), repos="https://cloud.r-project.org"); BiocManager::install("LEA")'
```

在仓库根目录执行。试跑预期得到 112 个体、28 组、149 个 LFMM2 位点和 159 个偏 RDA／默认偏相关位点。

```bash
(
  set -e
  export OPENBLAS_NUM_THREADS=1
  run=work/demo
  python3 04_genotype_environment_audit/workflow/prepare.py \
    04_genotype_environment_audit/example/run.json "$run/prepared"
  for pc in 3 5 10; do
    mkdir -p "$run/rda_$pc" "$run/lfmm_$pc"
    python3 01_partial_correlation/run.py "$pc" "$run/prepared" "$run/partial_$pc"
    Rscript 02_partial_RDA/run.R "$pc" "$run/prepared" "$run/rda_$pc"
    Rscript 03_LFMM2/run.R "$pc" "$run/prepared" "$run/lfmm_$pc"
  done
)
```

### 6. 正式运行、审核重建与输出

- 将 [正式配置模板](04_genotype_environment_audit/workflow/run.template.json) 复制到本地运行目录，填写 `kind=observed`、基因型路径和两张审核冻结表。JSON 中的相对输入路径从配置文件所在目录解析。可用 `dosage_tsv` 输入首列为 `sample`、其余列为实测二倍体 0/1/2 或空值的矩阵。单倍体或混合倍性使用 `vcf_bcf`。
- `samples.tsv` 每个个体一行，必需列为 `sample`、`group`。基因型文件必须能按精确键找到所有入组个体。`groups.tsv` 每组一行，组集合须与样本表完全一致，并包含配置中指定的有限数值环境列。启用经纬度调整需要 `lon`、`lat`，启用区块置换需要 `block`。模板要求每组至少三人和 PC10 模型正残差自由度。
- 使用示例步骤中的三个方法入口，将准备阶段的配置替换为正式 JSON，并选用新的运行目录。先看 `settings.json` 中的个体、组和位点范围，再解释方法结果。
- 重建审核时，`audit.py` 读取库存 TSV、上一版 Excel、原产地审核 Markdown、名称对照表、本地栅格、历史样本表和历史 BCF。检查逐个体 VCF 头、对齐精确 ID、重提环境并生成新版本的 `review/` 和 `environment/` 审核产物。

```bash
python3 04_genotype_environment_audit/audit.py --help
```

| 输出 | 首先查看 |
| --- | --- |
| `prepared/` | `settings.json` 样本与位点范围，`marker_catalog.tsv.gz` 位点来源及各方法纳入标记，`pca_eigenvalues.tsv` 方差比例，`group_design.tsv` 与 `individual_design.tsv` 冻结后的设计，`rda_frequency.tsv.gz` 与 `lfmm_dosage.tsv.gz` 响应矩阵 |
| `partial_*/` | `associations.tsv.gz` 的 rho、P、q、有数据组数、自由度和状态，`discoveries.tsv` 筛选结果，`testability.tsv` 可检验性计数 |
| `rda_*/` | `global.tsv` 的伪 F、整体 P 和解释率，`axes.tsv` 与 `marginal_terms.tsv`，`permutations.tsv.gz`，`site_scores.tsv`、`environment_arrows.tsv`、`locus_loadings.tsv.gz` 与 `exploratory_candidates.tsv` |
| `lfmm_*/` | `diagnostics.tsv` 的有效分母和基因组控制系数，`associations.tsv.gz` 的 beta、统计量、名义 P、校准 P 和 q，`discoveries.tsv` 筛选结果 |

RDA 和 LFMM2 的结果目录同时记录 `session_info.txt`。历史诊断来自 2026-11-18 主队列，历史 LFMM2 扫描的是潜因子 K，现行脚本比较的是显式 PC 数量。现行 RDA/LFMM2 入口源自 `workflowofSV` 的 `592540b`，偏相关复用 2026-11-18 修订版内核。使用与分享遵循 [LICENSE](LICENSE)。

## English

### 1. Purpose and data status

This repository supports genotype–environment association analysis of chicken structural variants (SVs) and pre-analysis review of genotype and environmental records. With observed SV calls, a reviewed individual list and a climate-group table, it runs partial rank correlation, partial redundancy analysis (partial RDA), and LFMM2. Each method can use the first 3, 5 or 10 components of the same SV PCA. The server holds the source VCF/BCF files, registration workbook and WorldClim rasters.

| Snapshot | Scope and use |
| --- | --- |
| Earlier registry, 2026-11-23 | 1,726 unique IDs. Baseline for the origin-evidence review and version comparison |
| Latest genotype inventory, 2026-11-23 | 1,732 unique IDs, including all previous individuals and 6 additions. Raw breed labels differ for 21 shared IDs and require interpretation |
| Historical main cohort, 2026-11-18 | 405 environmental records in a 565-individual BCF, assigned to 25 coordinate groups. Preserves actual historical climate inputs and diagnostics for all three methods |
| Synthetic validation cohort | 112 simulated individuals, 28 coordinate groups and 160 input alleles. All nine method-by-PC combinations passed the smoke run |

**Review status.** The recorded VCF sample IDs match all 1,732 primary VCF headers. All five climate variables were extracted for 66 tier A/B breed representative-point records. The current inventory contains 1,189 individuals with candidate climate awaiting batch-level review, 405 without complete candidate climate, 112 excluded commercial-cross individuals, 10 identity holds and 16 library-type holds. Header matching verifies the recorded sample column. Biological identity, genotype quality and environmental assignment require supporting evidence. [Audit counts](04_genotype_environment_audit/review/audit_summary.json)

### 2. Repository inventory

| Location | Files and purpose |
| --- | --- |
| Root | This `README.md`, the proprietary [LICENSE](LICENSE), `.gitignore` for local run outputs and source genotype files, and `.gitattributes` for PDF binary handling and TSV blank fields |
| [01_partial_correlation](01_partial_correlation/) | [README](01_partial_correlation/README.md), [run.py](01_partial_correlation/run.py) entry point, [partial_rank.py](01_partial_correlation/partial_rank.py) statistical kernel, and [historical_main405_summary.tsv](01_partial_correlation/historical_main405_summary.tsv) with historical test and discovery counts |
| [02_partial_RDA](02_partial_RDA/) | [README](02_partial_RDA/README.md), [run.R](02_partial_RDA/run.R) model and permutation entry point, [rda_metrics_vegan.tsv](02_partial_RDA/rda_metrics_vegan.tsv) historical global metrics, and [rda_axis_tests.tsv](02_partial_RDA/rda_axis_tests.tsv) historical axis tests |
| [03_LFMM2](03_LFMM2/) | [README](03_LFMM2/README.md), [run.R](03_LFMM2/run.R) fitting and calibration entry point, and historical per-environment [K5](03_LFMM2/historical_main405_K5_diagnostics.tsv), [K10](03_LFMM2/historical_main405_K10_diagnostics.tsv) and [K15](03_LFMM2/historical_main405_K15_diagnostics.tsv) diagnostics |
| [04_genotype_environment_audit](04_genotype_environment_audit/) | [Review guide](04_genotype_environment_audit/README.md) and [audit.py](04_genotype_environment_audit/audit.py) to rebuild the audit |
| [environment/](04_genotype_environment_audit/environment/) | [origin_evidence_1726_snapshot.tsv](04_genotype_environment_audit/environment/origin_evidence_1726_snapshot.tsv) for evidence classes and coordinates, [candidate_climate_by_breed.tsv](04_genotype_environment_audit/environment/candidate_climate_by_breed.tsv) for candidate point climates, [historical_main405_group_climate.csv](04_genotype_environment_audit/environment/historical_main405_group_climate.csv) for historical group-level raw values, and [historical_main405_individual_standardized_climate.csv](04_genotype_environment_audit/environment/historical_main405_individual_standardized_climate.csv) for standardized individual LFMM2 inputs |
| [review/](04_genotype_environment_audit/review/) | [audit_summary.json](04_genotype_environment_audit/review/audit_summary.json) with counts, states and source hashes, [individual_alignment.tsv](04_genotype_environment_audit/review/individual_alignment.tsv) with the full individual comparison, [review_queue.tsv](04_genotype_environment_audit/review/review_queue.tsv) with flagged cases, [breed_label_crosswalk.tsv](04_genotype_environment_audit/review/breed_label_crosswalk.tsv) with documented aliases, [historical_main405_alignment.tsv](04_genotype_environment_audit/review/historical_main405_alignment.tsv) with the historical alignment, and the origin review in [Markdown](04_genotype_environment_audit/review/origin_review_1726_20261123.md) and [PDF](04_genotype_environment_audit/review/鸡群分型登记对照与品种原产地核查_20261123.pdf) |
| [workflow/](04_genotype_environment_audit/workflow/) | [prepare.py](04_genotype_environment_audit/workflow/prepare.py) for shared preparation, [run.template.json](04_genotype_environment_audit/workflow/run.template.json) for production configuration, and the Chinese [methods and operations PDF](04_genotype_environment_audit/workflow/SV_PC_LFMM2_RDA_操作与原理.pdf) |
| [example/](04_genotype_environment_audit/example/) | [run.json](04_genotype_environment_audit/example/run.json) for the synthetic run, [dosage.tsv](04_genotype_environment_audit/example/dosage.tsv) for simulated dosages, [samples.tsv](04_genotype_environment_audit/example/samples.tsv) for individual-to-group assignment, and [groups.tsv](04_genotype_environment_audit/example/groups.tsv) for group climate and permutation blocks |

### 3. Review individuals and climate

1. Start with the [review queue](04_genotype_environment_audit/review/review_queue.tsv) for new IDs, raw label differences, header issues, missing climate and holds. Use the [full alignment](04_genotype_environment_audit/review/individual_alignment.tsv) to inspect other members of each batch. Individuals with candidate climate and no additional flags also remain subject to batch review.
2. Compare exact registry IDs, expected VCF IDs, observed VCF sample columns and source batches. Investigate raw label differences using the crosswalk and origin report. Review aliases, duplicates, disputed identities and library holds with their evidence. The current environmental cohort excludes confirmed commercial crosses, including XB and MH.
3. Review the coordinate purpose and evidence. Tier A contains curated origin representatives, B contains prior audited points, C contains locality leads awaiting coordinates, D requires finer geographic, subspecies or batch resolution, E identifies commercial crosses, and U marks identity conflicts. Assess whether each representative point applies to the source batch. Blank fields indicate missing data.
4. Freeze `samples.tsv` and `groups.tsv` after review. Record decisions, supporting evidence, the climate snapshot and any proposed permutation blocks. Use the historical 405-individual tables to reproduce the historical design. Build current model inputs from the current review.

Candidate climate values come from the containing cells of **WorldClim v2.1**, **1970–2000**, at **5 arc-minute** resolution. BIO3 measures isothermality, BIO15 precipitation seasonality, and BIO18 precipitation of the warmest quarter (mm). `wind` is mean monthly wind speed (m/s). `srad` is mean monthly solar radiation converted to W/m². These values characterize representative locations for geographic association analysis. [WorldClim data documentation](https://www.worldclim.org/data/worldclim21.html)

### 4. Methods and parameters

| Method | Unit, model and inference |
| --- | --- |
| Partial rank correlation | Coordinate-group SV allele frequency is the response. Rank the response and one environment at a time, then correlate their residuals after projecting out the selected SV PCs and optional geography. Apply BH across estimable loci for each environment and PC specification. The optional `--frequency` input can contain missing group frequencies and must match the group design in keys and order |
| Partial RDA | Coordinate groups are the observations. The response contains variable loci with estimable frequencies across all groups. Jointly fit the environments and condition on the selected SV PCs and optional geography. vegan reduced-model residual permutations provide global, marginal-term and axis tests. The standardized loading screen is applied to axes that pass the global-P and axis-BH gates |
| LFMM2 | Individuals are the observations. Fit complete, variable, observed SV dosage columns using joint environments, explicit SV PCs, optional geography and a fixed latent dimension. LEA estimates latent structure. Per-environment outputs include coefficients, nominal P, genomic-control P and BH q for each PC specification |

VCF/BCF records are split into ALT-specific columns. Individual dosage is `2 × ALT copies / called ploidy`. Group frequencies use total observed ALT copies divided by total called ploidy. PCA is fitted once to the observed SV matrix. Group PCs are averages of individual scores. Observed-mean projection is used for missing entries within the standardized PCA calculation. Association inputs use observed GT calls, and the workflow uses SV inputs only.

The production template sets `lfmm_k=5`, `permutations=999`, `min_called_per_group=3`, `adjust_geography=true`, `rda_blocked=false`, `alpha=0.05`, `loading_z=3` and `seed=20261123`. The three PC settings compare explicit structure adjustment at a fixed LFMM2 latent dimension. Choose free or blocked RDA permutations according to exchangeability among the coordinate groups. The synthetic example uses four blocks and 39 permutations.

Interpret models with environmental collinearity, residual degrees of freedom at PC10, LFMM2 calibration, within-group climate sharing and provenance in view. RDA loading tails provide exploratory locus evidence. Partial rank correlation and LFMM2 q-values refer to their respective testing families. [LFMM2 paper](https://doi.org/10.1093/molbev/msz008) · [RDA association study](https://doi.org/10.1111/mec.14584) · [vegan permutation documentation](https://vegandevs.github.io/vegan/reference/anova.cca.html)

### 5. Install and run the example

Use Linux, Python 3.11+, NumPy, pandas and SciPy ≥1.11, plus R with LEA, vegan, permute and jsonlite. Rebuilding the audit also requires Pillow, openpyxl and bcftools. Observed VCF/BCF input requires bcftools. Install LEA through Bioconductor. `SV_R_LIBRARY` optionally adds an R package library.

```bash
python3 -m pip install numpy pandas 'scipy>=1.11' Pillow openpyxl
Rscript -e 'install.packages(c("BiocManager", "vegan", "permute", "jsonlite"), repos="https://cloud.r-project.org"); BiocManager::install("LEA")'
```

Run from the repository root. The example yields 112 individuals, 28 groups, and 149 LFMM2 loci and 159 partial RDA/default partial rank loci.

```bash
(
  set -e
  export OPENBLAS_NUM_THREADS=1
  run=work/demo
  python3 04_genotype_environment_audit/workflow/prepare.py \
    04_genotype_environment_audit/example/run.json "$run/prepared"
  for pc in 3 5 10; do
    mkdir -p "$run/rda_$pc" "$run/lfmm_$pc"
    python3 01_partial_correlation/run.py "$pc" "$run/prepared" "$run/partial_$pc"
    Rscript 02_partial_RDA/run.R "$pc" "$run/prepared" "$run/rda_$pc"
    Rscript 03_LFMM2/run.R "$pc" "$run/prepared" "$run/lfmm_$pc"
  done
)
```

### 6. Production runs, audit rebuild and outputs

- Copy the [production template](04_genotype_environment_audit/workflow/run.template.json) to a local run directory. Set `kind=observed` and point to the observed genotypes and frozen cohort tables. The JSON file's directory is the base for relative input paths. The `dosage_tsv` format takes a leading `sample` column and observed diploid 0/1/2 or missing calls. Use `vcf_bcf` for haploid or mixed-ploidy calls.
- `samples.tsv` requires a unique `sample` key and `group` for every included individual. The genotype source must contain every included sample under its exact key. `groups.tsv` needs exactly the groups in the sample table and finite values in the configured climate columns. Geographic adjustment uses `lon` and `lat`. Blocked permutations use `block`. The template requires at least three individuals per group and positive residual degrees of freedom at PC10.
- Run the same preparation and model commands with the production JSON and a fresh output directory. Check the individual, group and locus counts in `settings.json` first.
- To rebuild the audit, supply the inventory TSV, prior Excel registry, origin review Markdown, breed crosswalk, local rasters, historical environment sample table and historical BCF to `audit.py`. It reads VCF headers, compares exact IDs, extracts climate and writes new `review/` and `environment/` audit files.

```bash
python3 04_genotype_environment_audit/audit.py --help
```

| Output | Key files |
| --- | --- |
| `prepared/` | `settings.json` for sample and locus scope, `marker_catalog.tsv.gz` for locus provenance and eligibility, `pca_eigenvalues.tsv` for explained variance, `group_design.tsv` and `individual_design.tsv` for the frozen designs, and `rda_frequency.tsv.gz` and `lfmm_dosage.tsv.gz` for the response matrices |
| `partial_*/` | `associations.tsv.gz` for rho, P, q, observed group counts, degrees of freedom and status, plus `discoveries.tsv` and `testability.tsv` |
| `rda_*/` | `global.tsv` for pseudo-F, global P and explained variation, `axes.tsv` and `marginal_terms.tsv`, `permutations.tsv.gz`, `site_scores.tsv`, `environment_arrows.tsv`, `locus_loadings.tsv.gz` and `exploratory_candidates.tsv` |
| `lfmm_*/` | `diagnostics.tsv` for the effective testing denominator and genomic-control factors, `associations.tsv.gz` for beta, test statistic, nominal and calibrated P and q, and `discoveries.tsv` |

Both R method directories contain `session_info.txt`. Historical diagnostics come from the 2026-11-18 main cohort. Historical LFMM2 varied the latent dimension K. The current scripts vary the number of explicit PCs. The current RDA and LFMM2 entry points derive from `workflowofSV` commit `592540b`. The partial rank kernel derives from the 2026-11-18 revision. See [LICENSE](LICENSE) for use and distribution terms.

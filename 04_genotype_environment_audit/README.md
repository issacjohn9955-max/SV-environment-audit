# 分型个体与环境表型审核

| 文件 | 审核用途 |
| --- | --- |
| [review/audit_summary.json](review/audit_summary.json) | 核对计数、来源校验和、VCF 头状态和历史队列匹配 |
| [review/individual_alignment.tsv](review/individual_alignment.tsv) | 当前清单与较早登记表逐个体对照、VCF 样本头、品种和环境证据 |
| [review/review_queue.tsv](review/review_queue.tsv) | 新增、原始品种标签差异、VCF 样本头差异和环境缺口的审核工作表 |
| [review/historical_main405_alignment.tsv](review/historical_main405_alignment.tsv) | 历史队列环境表与 565 BCF 样本头对照 |
| [environment/origin_evidence_1726_snapshot.tsv](environment/origin_evidence_1726_snapshot.tsv) | 1,726 个体登记版本的品种原产地证据 |
| [environment/candidate_climate_by_breed.tsv](environment/candidate_climate_by_breed.tsv) | 从 A/B 类行政区或产区代表点提取的五项候选气候值 |
| [environment/historical_main405_group_climate.csv](environment/historical_main405_group_climate.csv) | 旧 405 队列用于组水平模型的实值 |
| [environment/historical_main405_individual_standardized_climate.csv](environment/historical_main405_individual_standardized_climate.csv) | 旧 LFMM2 实际输入的逐个体标准化环境值 |

## 如何裁定

- 个体对照使用精确 ID，并保留新旧原始品种字符串。环境关联使用精确品种名及有文献或登记依据的[名称对照表](review/breed_label_crosswalk.tsv)。跨目录冲突和疑似文库类型问题暂挂。
- VCF 样本头检查读取每条登记的主分型文件，核对样本列。原始分型记录的完整性、基因型一致性和个体的实际身份需要独立证据。
- 环境采用 [WorldClim v2.1](https://www.worldclim.org/data/worldclim21.html) 的 1970–2000 年 5 弧分栅格。BIO3 为等温性，BIO15 为降水季节性，BIO18 为最暖季降水（mm）。wind 是 12 个月风速平均（m/s），srad 是 12 个月平均太阳辐射换算为 W/m²。缺测保留空白。
- A/B 代表点提取结果标为候选，建议先核定品种、来源批次和坐标用途。C/D/U 的待核地点和身份保留原证据。XB、MH 等 E 类商业配套系排除于这版主分析。两套历史环境提取表保留上次的队列范围。
- 审核输出是版本对照，另一个审核进程可以据此填写最终裁定。通过的个体名单和统一坐标组表才进入正式 PC3/5/10 分析。原始 VCF/BCF 在服务器本地。

`audit.py --help` 显示复建审核快照所需的输入。运行时请用新输出目录保存日期版本。原产地报告和审查来源见 [既有审核报告](review/origin_review_1726_20261123.md)。共享解析及模型说明见 [操作与原理 PDF](workflow/SV_PC_LFMM2_RDA_操作与原理.pdf)。

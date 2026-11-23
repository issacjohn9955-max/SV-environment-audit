# 鸡 SV—环境关联与输入审核

三种关联方法共用实测 SV 分型、冻结样本名单和环境组表。PC3、PC5、PC10 分别使用同一个 SV PCA 的前 3、5、10 个主成分。正式入组以逐个体审核通过的队列为准。

| 目录 | 内容 |
| --- | --- |
| [01_partial_correlation](01_partial_correlation/) | 坐标组频率逐环境偏秩相关，协变量为 SV PC 和可选经纬度 |
| [02_partial_RDA](02_partial_RDA/) | 坐标组频率的联合环境偏冗余分析、组水平置换与载荷 |
| [03_LFMM2](03_LFMM2/) | 实测个体剂量的联合环境 LFMM2，逐环境报告条件检验 |
| [04_genotype_environment_audit](04_genotype_environment_audit/) | 原产地证据、WorldClim 环境、基因型样本头核对、个体差异和审核入口 |

## 审核入口

从 [审核说明](04_genotype_environment_audit/README.md) 开始。当前分型清单、较早登记表、历史 405 队列与新坐标证据是不同版本。当前清单对应的环境逐行标明待核状态。A/B 级气候值来自品种代表点，个体归属及代表点对批次的适用性仍需审核。TSV 空白代表缺失。

## 试跑

Linux 下准备 Python 3.11+、NumPy、pandas、SciPy ≥1.11，以及 R 的 LEA、vegan、permute、jsonlite。复现环境审核还需要 Pillow、openpyxl 和 bcftools。可用 `SV_R_LIBRARY` 指定 R 附加包目录。

```bash
export OPENBLAS_NUM_THREADS=1
python3 04_genotype_environment_audit/workflow/prepare.py \
  04_genotype_environment_audit/example/run.json work/prepared
for pc in 3 5 10; do
  mkdir -p "work/rda_${pc}" "work/lfmm_${pc}"
  python3 01_partial_correlation/run.py "$pc" work/prepared "work/partial_${pc}"
  Rscript 02_partial_RDA/run.R "$pc" work/prepared "work/rda_${pc}"
  Rscript 03_LFMM2/run.R "$pc" work/prepared "work/lfmm_${pc}"
done
```

示例是 112 个模拟个体、28 个组。正式接口见 [配置模板](04_genotype_environment_audit/workflow/run.template.json)。新队列的结果应使用独立运行目录。RDA/LFMM2 代码来自 `workflowofSV` 的 `592540b` 并拆成独立入口。偏相关复用 2026-11-18 修订版本的偏秩相关内核。历史诊断表来自该次 405 个体队列，便于对照原分析口径。

[English](README.md) | [简体中文](README.zh-CN.md)

# 主图分析与绘图方法

本研究比较 Wagner 4 级与 1–3 级糖尿病足溃疡患者的肠道病毒组。分析脚本计算各分图所用的统计量，绘图脚本读取整理好的汇总表、图形配置或已有图稿。以下命令从仓库根目录运行。研究数据需另行提供，输出请写入单独的可写目录。[仓库概览](../README.zh-CN.md)

| 图号 | 分析或描述性汇总 | 绘图 |
| --- | --- | --- |
| **1b–f** | [figure1_catalogue_analysis.py](figure1_catalogue_analysis.py) · `summarize()` | [figure1_reflow.py](figure1_reflow.py) · `draw_panels()`，调用 [main_1.py](main_1.py) 中的 `bars()` |
| **2a–e** | [figure2_lifestyle_analysis.R](figure2_lifestyle_analysis.R) · `lifestyle_analysis()`、`same_genus_coverage()`、`figure2_sensitivity()` | [figure2_render.py](figure2_render.py) · `main_figure()` |
| **3a–c** | [figure3_paired_analysis.py](figure3_paired_analysis.py) · `analyze()` | [figure3_render.py](figure3_render.py) · `host_overview()` |
| **4a–e** | [figure4_host_analysis.R](figure4_host_analysis.R) · `paired_ranks()`、`prior_mag_score()`、`run_host_models()` | [figure4_render.py](figure4_render.py) · `focal_hosts()` |
| **5a–c** | [figure5_analysis.py](figure5_analysis.py) · `marker_span()`、`marker_consistency()`、`tree_neighbourhood()` | [figure5_reflow.py](figure5_reflow.py) · 命令行脚本 |

图 2 和图 4 共用 `fit_hc3()`、`hc3_result()` 进行回归推断。[hc3.R](hc3.R)

图 2–4 的绘图入口均需要已整理的 TSV 所在目录、对应图的 JSON 配置文件和输出目录。例如，图 2 使用 `python main/figure2_render.py DATA_DIR CONFIG_JSON OUTPUT_DIR`。图 3 和图 4 分别使用 `main/figure3_render.py` 和 `main/figure4_render.py`，参数顺序相同。三个脚本均导出 PDF、SVG 和 PNG。

## 图 1　病毒目录与分析人群

**1a** 为使用 PowerPoint 绘制的研究设计图。**1b–d** 分别按 CheckV 质量类别、病毒簇归属和可注释到的最深分类层级统计病毒操作分类单元（vOTU）数。**1e** 统计每个 vOTU 在多少名患者中检出，检出要求比对计数为正且覆盖广度至少为 50%。已完成测量但未检出任何 vOTU 的患者仍计入检出率分母，其病毒相对组成无定义。**1f** 分别列出两个 Wagner 组中纳入检出、组成和 ANCOM-BC2 分析的人数。

`python main/figure1_catalogue_analysis.py --help` 列出所需输入表和字段。将生成的 JSON 传给 `python main/figure1_reflow.py INPUT_JSON OUTPUT_DIR`，即可绘制 `Figure1b-f.pdf`。

## 图 2　病毒生活方式与溃疡严重程度

运行 `source("main/figure2_lifestyle_analysis.R")` 加载分析函数。**2a–c** 使用 `lifestyle_analysis(rpkm, category, meta)` 计算各生活方式类别的 RPKM 占该患者保留病毒目录总 RPKM 的百分比，再估计经年龄和性别调整的组间差值（Wagner 4 级减去 1–3 级，单位为百分点）。HC3 区间采用残差自由度。对每种目录方案，三种生活方式类别的检验共同进行一次 Benjamini–Hochberg（BH）校正。

**2d** 使用 `same_genus_coverage()` 统计符合两个条件的 vOTU 数，即最佳匹配的宏基因组组装基因组（MAG）与预测宿主同属，且比对覆盖查询序列的比例**严格大于 50%**。分母为该生活方式类别的全部 vOTU。**2e** 通过 `figure2_sensitivity(rpkm, catalog, meta)` 返回预测裂解型病毒占比的五种组间比较，依次使用原始目录、剔除预测宿主为 *Escherichia* 的 vOTU、放宽裂解型标签、调整 *E. coli* 协变量，或剔除宿主基因多于病毒基因的 vOTU。剔除目录成员后重新计算分母，放宽标签时沿用完整目录的分母。放宽标签和调整细菌协变量的两项比较报告原始 P 值，未计算 BH q 值。细菌协变量须在全队列中完成标准化。

## 图 3　病毒与细菌的丰度相关性

使用 `from main.figure3_paired_analysis import analyze` 导入分析函数。输入包括按覆盖广度过滤的病毒 RPKM 和计数、已有宿主指派、MAG 计数和属注释，以及置换和自助法各自的随机种子。序列比对和宿主预测结果需另行提供。

**3a** 对每个符合条件的宿主属，估计以其为预测宿主的病毒集合相对丰度与该细菌属相对丰度的相关性。**3b** 对满足检出门槛的单个 vOTU 作相同比较。两者均在病毒组成可计算的患者中求 Spearman 秩相关，通过置换细菌丰度对应的患者标签进行双侧检验。仅属水平计算 95% 区间，方法是按患者对病毒和细菌数据配对重抽样。**3c** 分别统计预测宿主匹配和不匹配的配对中 q < 0.05 的关联数。属水平、单个 vOTU 和完整交叉矩阵各自作 BH 校正，每个矩阵内的匹配与不匹配配对共同校正。将剔除的 vOTU 编号传入 `excluded`，即可在相同患者及细菌丰度分母下，重新计算病毒相对丰度、筛选符合条件的集合并构建剔除后的矩阵。函数返回属水平、单个 vOTU、交叉矩阵及配对汇总表。

## 图 4　重点宿主属

运行 `source("main/figure4_host_analysis.R")` 加载分析函数。**4a–b** 对 *Escherichia* 与 *Collinsella* 分别调用 `paired_ranks()`。每个属均合并两个 Wagner 组的患者，分别对病毒和细菌丰度排序，保留有效零值。图中直线描述秩之间的关系，相关系数、区间和 q 值来自图 3 的属水平分析。此处的细菌丰度按属汇总其全部 MAG。

**4c–e** 使用 `run_host_models()`。4c 以双侧 Fisher 精确检验比较两个 Wagner 组的检出率。4d 和 4e 仅纳入病毒丰度为正的患者，以病毒相对丰度的 log2 值为因变量，采用 HC3 线性回归，分别估计严重程度组间差异以及病毒丰度与细菌复合分数之间的关联。BH 校正的家族包括八项检验，纳入其中可估计的 P 值。每个阳性丰度模型要求每组至少有 5 名资料完整的阳性患者、两组合计至少 20 名，且残差自由度至少为 10。

筛选阳性患者前，先以 `prior_mag_score()` 构建细菌评分。每名患者的三个目标 MAG 均以其占该患者全部 MAG 计数的百分比表示。各百分比加上这三个 MAG 在全体患者中的最小正百分比的一半后取 log10，再按每个 MAG 在全队列中分别标准化。评分等于 *E. coli* MAG 的 z 值减去两个 *Collinsella* MAG 的平均 z 值。

将五份输入表（`rpkm.tsv`、`counts.tsv`、`metadata.tsv`、`bacteria.tsv`、`hosts.tsv`）置于单独的可写目录，再运行 `Rscript main/figure4_host_analysis.R EXTERNAL_WORKDIR`。命令将检验结果、集合成员、患者指标、诊断结果及已拟合的模型对象写入该目录。

## 图 5　目标病毒代表序列的证据

`python -m main.figure5_analysis --help` 列出三个子命令。**5a** 使用 `span` 从已有注释表汇总一段连续标志基因的坐标范围。**5b** 使用 `markers` 比较各标志基因中排名第一且达到支持度要求的分支划分，检查其中含目标序列的一侧是否一致，一致时按类别统计该侧成员。**5c** 使用 `neighbourhood`，按已有系统发育树上带枝长的路径距离选择要展示的末端。每个命令输出一份 JSON 汇总。全长 MAG 匹配、系统发育树、分支支持度及原始图稿均来自已有结果。

对基因组与系统发育树图稿重新排版时，使用结构兼容的源 SVG 和标签文件运行 `python main/figure5_reflow.py SOURCE_SVG LABEL_JSON OUTPUT_DIR`。脚本保留树的绘图路径，并导出 SVG、PDF 和 PNG。

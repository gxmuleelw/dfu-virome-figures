[English](README.md) | [简体中文](README.zh-CN.md)

# 表格分析方法

[仓库总览](../README.zh-CN.md) · [主图分析方法](../main/README.zh-CN.md)

本目录的脚本使用外部准备好的计数矩阵、注释、模型对象或系统发育树。与主图共用方法的表格，下文直接指向相应主图的分析代码。

## 主表 1 和主表 2

| 主表 | 数值来源 |
| --- | --- |
| 表 1：宿主相关病毒集合的主要检验 | [图 4 分析](../main/figure4_host_analysis.R)中的 `run_host_models()`：检出率用双侧 Fisher 精确检验；在病毒丰度为正的患者中，以病毒相对丰度的 log2 值为因变量拟合 HC3 线性回归；BH 校正的家族包括八项检验。 |
| 表 2：两个重点病毒代表序列的证据与分类归属 | 已有的基因注释、MAG 比对和标记基因树。各标记间含目标序列分组的比较，以及按树上距离最近的正式参考序列，分别由[图 5 分析](../main/figure5_analysis.py)中的 `marker_consistency()` 和 `tree_neighbourhood()` 汇总。 |

## 样本与病毒目录的汇总

| 脚本与函数 | 输入、计算及结果含义 |
| --- | --- |
| [比对质量控制](alignment_qc.py) · `summarize()` | 输入行列对应的比对计数、RPKM 和覆盖广度矩阵（行是 vOTU，列是样本，覆盖广度用比例表示），另提供每个样本质控后的读段总数及与列顺序一致的样本标签。函数分别汇总筛选前的计数，以及覆盖广度 ≥0.5 的计数和 RPKM，并核对筛选前的比对计数是否超过质控后的读段总数。仅当筛选后计数和 RPKM 的样本总和均为正，才可计算该患者的病毒相对组成。 |
| [检出率筛选](prevalence.R) · `prevalence_filter()` | 输入计数矩阵、判定检出所需的最低计数、检出率阈值及方案名称。函数逐个 vOTU 统计达到最低计数要求的样本数，与“输入样本数 × 检出率阈值”向上取整后的值比较，记录是否保留及原因。分母为本次传入矩阵的列数。 |
| [α 多样性](diversity.R) · `alpha_values()`、`diversity_values()` | 输入每个样本合计为 1 的相对丰度矩阵，以及行名、列名均对应的比对计数矩阵。函数计算 Shannon 指数（自然对数）、逆 Simpson 指数和计数大于零的 vOTU 数。`alpha_values()` 要求每个输入样本的比对计数总和为正。 |
| [宿主属预测结果整理](host_assignment.py) · `assign()`、`number15()`、`text_bool()` | 输入病毒目录的 vOTU 标识和含 `Virus`、`Host genus`、`Confidence score` 三列的 iPHoP 记录。保留评分 ≥90 的记录。若最高分仅对应一个属，即使该属多次出现，仍归入该属。若多个属并列最高分，记为 `Ambiguous`。无达标记录则记为 `Unassigned`。两个格式化函数分别处理数值评分和布尔标记。所得宿主属为计算预测结果。 |

## 模型与相关性分析

| 代码 | 输入、方法及输出 |
| --- | --- |
| [图 2 分析](../main/figure2_lifestyle_analysis.R) · `lifestyle_analysis()` | 在 R 中载入该脚本，向 `lifestyle_analysis()` 传入按覆盖广度筛选的 RPKM、以 vOTU 标识命名的类别向量和患者资料。先计算每名患者各类别丰度占比，再用 `safe_lifestyle_fit()` 比较严重程度组，并校正年龄和性别。HC3 置信区间使用残差自由度。每套目录方案内的三个类别构成一组 BH 校正检验。`figure2_sensitivity()` 汇总五种比较，`same_genus_coverage()` 根据已有的最佳 MAG 匹配，统计匹配属一致且查询序列覆盖率严格大于 50% 的比例。 |
| [图 3 分析](../main/figure3_paired_analysis.py) · `analyze()` | 相关分析用到五个函数。`unit_ranks()` 对并列值取平均秩并将秩向量归一化，`rho()` 计算 Spearman 相关系数，`permutation()` 置换细菌数据的患者标签进行双侧检验，`boot_ci()` 以配对患者为单位进行自助抽样，`bh()` 按传入的家族大小校正 P 值。完整分析调用 `analyze()`，传入按患者对齐的病毒计数与 RPKM、宿主属归属、MAG 计数及属注释，以及两个随机种子。仅属水平相关系数计算自助法置信区间。 |
| [图 4 分析](../main/figure4_host_analysis.R) · `run_host_models()` | `prior_mag_score()` 依据全队列的 MAG 计数构建细菌复合分数。`paired_ranks()` 针对一个属，合并两个严重程度组，保留有效的零值后分别计算两种丰度的秩。`run_host_models()` 读取五个 TSV 文件，用双侧 Fisher 精确检验比较组间检出情况，并在丰度为正的患者中对病毒相对丰度取 log2 后拟合 HC3 线性模型。可估计的 P 值在八项检验构成的家族内进行 BH 校正。 |
| [共线性诊断](ecology_diagnostics.R) · `vif_table()` | 输入含截距、且有多个预测变量列的已拟合模型及模型标签。函数对设计矩阵中每个非截距列，以其余预测变量列为自变量作回归，返回方差膨胀因子 `1/(1 − R²)`。若一个分类变量展开为多列，每列各有一个结果。 |

图 4 的模型在仓库根目录运行 `Rscript main/figure4_host_analysis.R EXTERNAL_WORKDIR`。在单独的可写目录中准备 `rpkm.tsv`、`counts.tsv`、`metadata.tsv`、`bacteria.tsv` 和 `hosts.tsv`。检验结果、病毒集合成员表、样本指标、诊断结果和拟合模型对象均写入该目录。图 2 和图 4 共用 [HC3 函数](../main/hc3.R)，依赖 R 包 `sandwich`。相关性分析依赖 NumPy、pandas 和 SciPy。

## 树上距离与坐标区间

[树和区间汇总脚本](tree_summary.py) 提供两个函数。`neighbourhood()` 读取一棵已有的、带分支长度的树，另需目标末端、各末端的类别，以及希望展示的近邻数和 ICTV 参考序列数。它返回入选末端的集合，并为每个末端记录至目标的距离和展示理由。距离相同时按标识符排序，此顺序仅用于排版。[图 5 分析](../main/figure5_analysis.py) 通过此函数选取局部树图中的末端。`union_bases()` 读取包含起止位置的闭区间，合并重叠或相邻区间，返回覆盖的碱基数及合并后的区间。建树和分支支持度估计均在上游完成。

[English](README.md) | [简体中文](README.zh-CN.md)

# 肠道病毒组与糖尿病足溃疡严重度

本仓库为论文 “Lower relative abundance of lytic gut phages in patients with Wagner grade 4 diabetic foot ulcers” 的代码。

本研究比较了 Wagner 4 级与 1–3 级糖尿病足溃疡患者的肠道病毒组。下表按研究问题列出每张主图的两个入口：“分析或描述性汇总”说明图中的量如何算出，“绘图”说明这些量如何呈现。仓库收录代码和说明；研究数据、模型拟合结果及原始矢量图稿另行保存。

| 主图与研究问题 | 分析或描述性汇总 | 绘图 |
| --- | --- | --- |
| **图 1b–f**：病毒目录由哪些成员构成？各项分析分别纳入了哪些患者？ | [统计目录类别、检出频数和各分析的纳入人数](main/figure1_catalogue_analysis.py) | [病毒目录与队列分图](main/figure1_reflow.py) |
| **图 2**：按预测生活方式划分的噬菌体相对丰度是否随 Wagner 分级而异？ | [计算每位患者的丰度百分比，估计调整后的组间差值，并进行敏感性比较](main/figure2_lifestyle_analysis.R) | [分组分布与置信区间](main/figure2_render.py) |
| **图 3**：病毒丰度与其预测宿主属的相关性，与其他细菌属相比有何不同？ | [构建属集合及单个病毒成员的配对，置换患者标签，仅在属集合层面计算 bootstrap 区间，再比较匹配与非匹配配对](main/figure3_paired_analysis.py) | [相关系数与配对计数](main/figure3_render.py) |
| **图 4**：预测宿主为 *Escherichia* 或 *Collinsella* 的病毒集合，其检出和丰度与溃疡严重度、细菌丰度有何关联？ | [计算配对丰度的秩，比较检出率，构建细菌复合分数，并在丰度为正的样本中拟合模型](main/figure4_host_analysis.R) | [检出与丰度分图](main/figure4_render.py) |
| **图 5**：已有基因注释和系统发育树为代表性病毒提供了哪些证据？ | [汇总标记基因坐标，比较已核定的分支划分，并在已有树上选取邻近树梢](main/figure5_analysis.py) | [基于现有矢量图稿排版](main/figure5_reflow.py) |

图 1a 为使用 PowerPoint 手绘的研究设计示意图。

## 按目录阅读

- [主图](main/README.zh-CN.md)：按分图介绍分析方法、所需输入、关键函数和绘图命令。
- [补充图](supp/README.zh-CN.md)：查找图 S1–S14 的绘图脚本及共用绘图函数。
- [表格](table_analysis/README.zh-CN.md)：主表 1、2 数值的来源，以及各项表格指标的计算方法。

## 数据与软件

原始宏基因组测序数据已存入 NCBI Sequence Read Archive，登录号 PRJNA1471479。统计分析使用 R 4.5.0 和 Python 3.12.10，软件包版本见补充方法第 11 节的关键资源表。

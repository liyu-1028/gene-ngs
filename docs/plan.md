# 全流程实施计划：合成 NGS 数据生成系统 (Full Implementation Plan)

本文档规划了系统五个阶段的实施步骤，全部在 **Ubuntu 22.04 LTS** 生态环境下运行，确保研究的可复用性与环境的纯洁一致。

---

## 🏗️ 部署环境概览

### **1. 运行环境 (Ubuntu 22.04 LTS)**
*   **任务**: 运行 Python 主逻辑、训练 CTGAN (利用 CPU/GPU)、调用 DeepSeek API、执行生信仿真与质控、进行最终评估。
*   **环境**: Python 3.10, PyTorch, SDV, Synthcity, MinerU (PDF 解析), Synggen (生信仿真), Sandy (质控模拟), Samtools, BWA.

---

## 🎯 第一阶段：数据治理 (Data Governance)
**目标**：将 24 份非结构化 PDF 转化为标准化的特征 CSV 文件。

*   **步骤 1**: 在 Ubuntu 本地部署并运行 MinerU HTTP 服务，用于解析 PDF 并返回 Markdown 文本。 (✅ 已完成)
*   **步骤 2**: 运行 `src/stage1/batch_processor.py` 发起批量解析。 (✅ 已完成)
*   **步骤 3**: 运行 `src/stage1/run_extraction.py` 通过大模型提取关键临床特征。 (✅ 已完成)

---

## 🎯 第二阶段：混合生成 (Hybrid Generation)
**目标**：基于 24 份真实样本，合成 1000 份虚拟患者数据。

*   **步骤 1**: 运行 `src/stage2/ctgan_generator.py` 训练 CTGAN 模型并生成连续特征。 (✅ 已完成)
*   **步骤 2**: 运行 `src/stage2/offspring_generator.py`（集成在 `run_generation.py` 中）实现变异特征的交叉重组。 (✅ 已完成)
*   **结果**: 生成了 `data/synthetic/raw_synthetic_data.csv`。

---

## 🎯 第三阶段：临床逻辑审计 (Clinical Auditing)
**目标**：剔除生成数据中的医学逻辑错误，确保“生物学合理性”。

*   **步骤 1**: 运行审计脚本 `src/stage3/clinical_auditor.py`：
    *   **硬规则审计**: 检查 TMB-MSI 矛盾、VAF 异常。 (✅ 已完成)
    *   **LLM 深度审计**: 调用 **DeepSeek-Chat**，利用推理能力识别基因互斥和临床矛盾。 (✅ 已完成)
*   **步骤 2**: 得到通过审计的高质量记录，保存为 `data/synthetic/audited_synthetic_data.csv`。 (✅ 已完成)

---

## 🎯 第四阶段：Read-Level 物理仿真 (Simulation)
**目标**：将表格变异注入底层测序数据，生成高仿真 FASTQ 文件。

*   **步骤 1**: 编译 **Synggen**。进入 `src/stage4/synggen_src` 执行 `make`。 (✅ 已完成)
*   **步骤 2**: 运行 `src/stage4/simulation_orchestrator.py`：
    *   **核心逻辑**:
        1. **解析变异**: 从 `audited_synthetic_data.csv` 提取 SNV/Indel/CNV。
        2. **坐标映射**: 自动对齐 GRCh38 坐标并处理染色体命名规范 (chr1 vs 1)。
        3. **指令生成**: 生成 `variants.pm` 等 Synggen 专用变异描述文件。 (✅ 已完成)
*   **步骤 3**: 执行物理注入流水线：
    *   **背景建模**: 利用 HG001 (NA12878) 标准品训练测序仪背景模型。 (✅ 已完成)
    *   **批量仿真**: 运行 `run_simulation.sh`，为 1000 个虚拟病人生成独立的 FASTQ 文件夹。 (✅ 已完成)

---

## 🎯 第五阶段：保真度与可用性评估 (Evaluation)
**目标**：定量证明合成数据集的科学价值、生物学一致性与隐私合规性，达到《Scientific Data》期刊的发表标准。

*   **步骤 1: 物理层测序保真度 (Read-Level NGS Fidelity)**
    *   **任务**: 验证仿真生成的 FASTQ/BAM 文件是否准确携带了设定的生物学变异。
    *   **操作**: 运行生信验证脚本 `src/stage5/verify_fidelity.sh`，抽取虚拟患者的数据，统计突变位点的 K-mer 匹配数，计算观测到的等位基因丰度 (VAF)。
    *   **指标**: 计算变异注入的召回率 (Recall) 和 VAF 绝对误差。
*   **步骤 2: 统计与生物学保真度 (Statistical & Tabular Fidelity)**
    *   **任务**: 评估 `audited_synthetic_data.csv` 与真实临床特征分布的一致性。
    *   **操作**: 运行 `src/stage5/statistical_fidelity.py`。
    *   **指标**: 
        *   连续变量 (TMB, 年龄): Kolmogorov-Smirnov (KS) 检验。
        *   分类变量 (基因, MSI): 总变分距离 (TVD)。
*   **步骤 3: 下游任务实用性 (Clinical Utility & Downstream Performance)**
    *   **任务**: 证明合成数据可以代替真实数据用于医学机器学习预训练。
    *   **操作**: 运行 `src/stage5/utility_tstr.py`。
    *   **指标**:
        *   TSTR (Train on Synthetic, Test on Real) 性能逼近度：比较 TSTR 模型与 TRTR (基线) 模型的 AUROC 和 F1-score，目标差距 $\Delta \le 0.05$。
        *   SHAP 特征稳定性：计算真实与合成模型核心决策特征排名的 Spearman 相关系数 ($\rho > 0.8$)。
*   **步骤 4: 隐私合规与防泄露 (Privacy Preservation)**
    *   **任务**: 从数学上证明合成数据没有泄露最初 24 位患者的隐私。
    *   **操作**: 运行 `src/stage5/privacy_assessment.py`。
    *   **指标**:
        *   精确匹配得分 (Exact Match Score) 必须为 0。
        *   最近邻距离 (DCR) 与成员推理攻击 (MIA) 防御度验证。

---

## 🖥️ Ubuntu 22.04 虚拟机服务部署手册 (环境已就绪)
...
[此处省略中间内容]
...

---

## 📈 实施进度跟踪 (Progress Tracking)

| 阶段 | 任务描述 | 状态 | 备注 |
| :--- | :--- | :--- | :--- |
| **环境** | Ubuntu 22.04 LTS 生态环境 | ✅ 已完成 | 所有 Python 依赖、API、生信工具均在本地打通 |
| **第一阶段** | 数据结构化提取 (PDF->CSV) | ✅ 已完成 | 24 份报告已提取完成 |
| **第二阶段** | 混合生成 (CTGAN + Offspring) | ✅ 已完成 | 1000 条合成数据已生成 |
| **第三阶段** | 临床逻辑审计 (Rule + LLM) | ✅ 已完成 | 成功产出通过严苛临床规则审计的洁净表格数据集 |
| **第四阶段** | 测序物理仿真流水线 | ✅ 已完成 | 物理变异注入完毕，成功产出虚拟患者 FASTQ 序列 |
| **第五阶段** | 保真度与可用性评估 | ✅ 已完成 | 物理、统计、实用性与隐私评估全部通过并输出指标 |

---

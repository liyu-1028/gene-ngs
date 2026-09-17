# Stage 0: 数据源接入、前置处理与运维脚本

> 运维脚本 (resume/flatten/verify/overnight/zen*) 使用 `PROJECT_ROOT` 环境变量定位
> 项目根目录，默认 `$HOME/gene-ngs`，可用 `export PROJECT_ROOT=/path/to/gene-ngs` 覆盖。


> v1 见 git 历史：数据源为 msk_impact_2017 (ODbL, 10,945 样本)。
> v2 切换为 msk_impact_50k_2026，并在 `src/stage0/` 维护因原 stage3/4 目录
> root 权限冻结而无法就地升级的 v2 脚本。

## 数据源

| 项 | v2 值 |
|---|---|
| studyId | `msk_impact_50k_2026` (cBioPortal DataHub) |
| 文献 | Bandlamudi et al., *Cancer Cell* 2026, PMID 41895280 (MSK-50K) |
| 队列 | 54,331 肿瘤样本 + 配对正常；IMPACT341/410/468/505 |
| 压缩包 SHA256 | `c9522b3ed7cfcf803fe2838868bdc566166c45a603de7962859416c24c9e7d38` |
| 许可 | **CC BY-NC-ND 4.0**（⚠ 与 2017 版 ODbL 不同，发布衍生数据前需确认 ND 条款） |
| 上游原始文件 | `data/msk_impact_50k_2026/`（未改动） |

相对 2017 版的关键增益：`AGE_AT_DX`（覆盖 96.65%）、`MSI_SCORE`/`MSI_TYPE`（86.95%）、
连续 `TMB_SCORE`、`FACETS_QC`（真实 Pass/Fail QC）。

## 脚本清单与执行顺序

```
1. src/stage0/convert_msk_impact_50k_to_pipeline.py      # 数据源 -> 流水线标准模式
   -> data/real/msk_impact_50k_2026/real_dataset.csv           (54,331 行 × 10 列)
   -> data/real/msk_impact_50k_2026/real_dataset_clinical_wide.csv
2. src/stage0/build_mutation_coordinate_db.py             # stage4 坐标映射库重建
   -> configs/mutation_coordinates_50k.json                    (222,197 SNV + 501 CNA, GRCh38)
   # 依赖: pyliftover (PYTHONPATH=.local_pkgs) + UCSC hg19ToHg38 chain
3. src/stage2/run_generation.py (v2)                      # CTGAN 1:1 生成 54,331 条
   -> data/synthetic/raw_synthetic_data.csv
4. src/stage0/clinical_auditor_v2.py                      # 硬规则全量 + LLM 抽样 2000
   -> data/synthetic/audited_synthetic_data.csv
5. src/stage0/run_stage4_v2.py                            # FASTQ 子集 2000 + Synggen 批处理
   -> data/synthetic/variant_specs/ + run_simulation.sh
   -> bash data/synthetic/run_simulation.sh                   (后台, 每患者 ~20 MB)
```

## 关键规模决策

| 层 | 规模 | 依据 |
|---|---|---|
| 表格层 (stage2/3/5) | 54,331（与真实 1:1） | "对应数据量级"；CTGAN 以 52,204 完整病例训练 |
| LLM 深审 (stage3) | 随机抽样 2,000 (seed=42) | 5 万条 × API 调用不可行；硬规则仍全量 |
| FASTQ 仿真 (stage4) | 随机抽样 2,000 (seed=42) | 每患者 ≈20 MB，全量需 ~1 TB > 磁盘 |

## stage3 v2 审计规则说明

- **双向 TMB-MSI 交叉校验**（50K 有真实 MSI，恢复完整规则）：
  - `TMB > 20 且 MSS` → REJECT（原单向规则）
  - `MSI-H 且 TMB < 2` → REJECT（v2 新增反向）
- **TMB 单边校验**仅当数据源无 MSI 字段时作为降级方案（2017 版曾被迫采用），
  现已不再需要。
- VAF 范围检查 (0,1] 保持。

## stage4 v2 坐标库说明

旧库 12 SNV + 7 CNA，且混装 GRCh37 (EGFR)/GRCh38 (MYC/TERT) 坐标。
新库统一 hg19→hg38 liftover (UCSC chain)，热点验证：
`EGFR:p.L858R -> chr7:55,191,822`、`BRAF:p.V600E -> chr7:140,753,336` ✓
CNA 跨度取自 GRCh38 target_regions.bed 基因外显子聚合 (501/541 命中)。

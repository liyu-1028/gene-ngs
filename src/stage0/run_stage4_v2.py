#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
run_stage4_v2.py  (Stage 4 v2, MSK-IMPACT 50K 量级)

复用原 src/stage4/simulation_orchestrator.py 的 SimulationOrchestrator 类
(root 冻结，不可改)，v2 的差异:

1. 坐标映射库切换为 configs/mutation_coordinates_50k.json
   (222,197 SNV + 501 CNA, 统一 GRCh38; 旧库仅 12+7 且混装 GRCh37/38)
2. FASTQ 物理仿真规模: 从审计通过的全量 (5 万+) 中随机抽取 FASTQ_N 个
   虚拟患者 (默认 2000, seed=42)。磁盘约束: 每患者 FASTQ ≈ 20 MB,
   全量 5.4 万患者需 ~1 TB，超出本机磁盘，故采用随机子集；
   表格层数据仍为全量，stage5 的统计/实用性/隐私评估不受影响。
3. 患者编号改用 audited 数据的原始行号，保证可追溯。

用法:
   python src/stage0/run_stage4_v2.py [--fastq-n 2000]
"""

import argparse
import os
import sys

import pandas as pd

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../..')))
from src.stage4.simulation_orchestrator import SimulationOrchestrator   # 复用原实现


class SimulationOrchestratorV2(SimulationOrchestrator):
    """增加 reads_root 参数：不同数据源路线的 FASTQ 产出隔离到不同目录"""

    def __init__(self, *a, reads_root='data/synthetic/reads', **kw):
        super().__init__(*a, **kw)
        self.reads_root = reads_root

    def generate_bash_script(self, var_base_dir, patient_list):
        script_content = f"""#!/bin/bash
# Batch simulation script with unique seeds (Dynamic Absolute Paths)
PROJECT_BASE=$(cd "$(dirname "${{BASH_SOURCE[0]}}")/../.." && pwd)
READS_ROOT="$PROJECT_BASE/{self.reads_root}"

PATIENTS=({ " ".join(patient_list) })

for ((i=0; i<${{#PATIENTS[@]}}; i++)); do
    PT_ID="${{PATIENTS[i]}}"
    SEED=$((42 + i))
    echo "Simulating $PT_ID with seed $SEED..."
    OUT_DIR="$READS_ROOT/$PT_ID"
    mkdir -p "$OUT_DIR"
    VAR_DIR="{var_base_dir}/$PT_ID"

    "$PROJECT_BASE/bin/synggen" mode=1 \\
        bed="$PROJECT_BASE/data/raw/reference/target_regions.bed" \\
        fasta="$PROJECT_BASE/data/raw/reference/GRCh38_full_analysis_set_plus_decoy_hla.fa" \\
        rdm="$PROJECT_BASE/data/raw/reference/models/wes_model/target_regions.rdm.gz" \\
        pbe="$PROJECT_BASE/data/raw/reference/models/wes_model/target_regions.pbe.gz" \\
        qm="$PROJECT_BASE/data/raw/reference/models/wes_model/target_regions.qm.gz" \\
        pm="$VAR_DIR/variants.pm" \\
        indel="$VAR_DIR/variants.indel" \\
        cna="$VAR_DIR/variants.cna" \\
        threads=4 \\
        nreads=200000 \\
        tc=0.3 \\
        out="$OUT_DIR" \\
        seed=$SEED
done
echo "Batch simulation completed."
"""
        script_path = os.path.join(self.output_dir, 'run_simulation.sh')
        with open(script_path, 'w', newline='\n') as f:
            f.write(script_content)
        print(f"Generated script: {script_path} (reads -> {self.reads_root})")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--input', default='data/synthetic/audited_synthetic_data.csv')
    ap.add_argument('--output-dir', default='data/synthetic')
    ap.add_argument('--reads-root', default='data/synthetic/reads')
    ap.add_argument('--config', default='configs/mutation_coordinates_50k.json')
    ap.add_argument('--fastq-n', type=int, default=2000,
                    help='进入 FASTQ 物理仿真的虚拟患者数 (默认 2000)')
    ap.add_argument('--seed', type=int, default=42)
    args = ap.parse_args()

    orch = SimulationOrchestratorV2(data_path=args.input,
                                    output_dir=args.output_dir,
                                    config_path=args.config,
                                    reads_root=args.reads_root)
    df = orch.load_data()
    print(f"Audited records: {len(df):,}")

    # 随机抽取进入 FASTQ 仿真的子集
    if args.fastq_n and args.fastq_n < len(df):
        sub = df.sample(n=args.fastq_n, random_state=args.seed)
        print(f"FASTQ simulation subset: {args.fastq_n:,} / {len(df):,} "
              f"(seed={args.seed}); 全量需 ~{len(df)*32/1e3:.0f} GB > 本机磁盘，采用子集")
    else:
        sub = df
        print(f"FASTQ simulation on ALL {len(sub):,} patients")

    sub = sub.reset_index().rename(columns={'index': 'src_row'})
    patient_data = orch.process_variants(sub)
    var_base_dir, pt_list = orch.write_synggen_files(patient_data)
    orch.generate_bash_script(var_base_dir, pt_list)
    print(f"Ready. {len(pt_list):,} patients queued for Synggen simulation.")
    print("Run: nohup bash data/synthetic/run_simulation.sh > data/synthetic/simulation.log 2>&1 &")


if __name__ == '__main__':
    main()

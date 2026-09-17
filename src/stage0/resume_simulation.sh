#!/bin/bash
# resume_simulation.sh — Synggen 断点续跑（通用）
# 用法: resume_simulation.sh <患者清单文件> <specs目录> <reads根目录> <日志文件>
# 种子规则与原脚本一致: seed = 42 + 患者编号
set -u
PROJECT_BASE=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
LIST_FILE="$1"; VAR_SPECS="$2"; READS_ROOT="$3"; LOG_FILE="$4"

: > "$LOG_FILE"
N_TOTAL=$(wc -l < "$LIST_FILE")
N_DONE=0

while read -r PT_ID; do
    [ -z "$PT_ID" ] && continue
    IDX=${PT_ID#virtual_pt_}
    IDX=$((10#$IDX))
    SEED=$((42 + IDX))
    N_DONE=$((N_DONE + 1))
    echo "[resume $N_DONE/$N_TOTAL] Simulating $PT_ID with seed $SEED..." | tee -a "$LOG_FILE"
    OUT_DIR="$PROJECT_BASE/$READS_ROOT/$PT_ID"
    mkdir -p "$OUT_DIR"
    "$PROJECT_BASE/bin/synggen" mode=1 \
        bed="$PROJECT_BASE/data/raw/reference/target_regions.bed" \
        fasta="$PROJECT_BASE/data/raw/reference/GRCh38_full_analysis_set_plus_decoy_hla.fa" \
        rdm="$PROJECT_BASE/data/raw/reference/models/wes_model/target_regions.rdm.gz" \
        pbe="$PROJECT_BASE/data/raw/reference/models/wes_model/target_regions.pbe.gz" \
        qm="$PROJECT_BASE/data/raw/reference/models/wes_model/target_regions.qm.gz" \
        pm="$PROJECT_BASE/$VAR_SPECS/$PT_ID/variants.pm" \
        indel="$PROJECT_BASE/$VAR_SPECS/$PT_ID/variants.indel" \
        cna="$PROJECT_BASE/$VAR_SPECS/$PT_ID/variants.cna" \
        threads=4 \
        nreads=200000 \
        tc=0.3 \
        out="$OUT_DIR" \
        seed=$SEED
done < "$LIST_FILE"
echo "Resume batch completed: $N_TOTAL patients." | tee -a "$LOG_FILE"

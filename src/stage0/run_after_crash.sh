#!/bin/bash
# run_after_crash.sh — 崩溃恢复链: 2017 全量仿真(优先/发表路线) → 50K 断点续跑
# 用法: nohup bash src/stage0/run_after_crash.sh > crash_recovery.log 2>&1 &
set -e
cd "$(dirname "$0")/../.."

echo "=== [$(date '+%F %T')] Phase 1: 2017 发表路线 FASTQ 全量仿真 (2000 患者) ==="
bash data/synthetic_2017/run_simulation.sh > data/synthetic_2017/simulation.log 2>&1
echo "=== [$(date '+%F %T')] Phase 1 完成: 2017 FASTQ ==="

echo "=== [$(date '+%F %T')] Phase 2: 50K FASTQ 断点续跑 (1183 患者) ==="
bash src/stage0/resume_simulation.sh \
    /tmp/50k_missing.txt \
    data/synthetic/variant_specs \
    data/synthetic/reads \
    data/synthetic/simulation_resume.log
echo "=== [$(date '+%F %T')] Phase 2 完成: 50K FASTQ 续跑 ==="
echo "ALL DONE."

#!/bin/bash
# run_2017_route.sh — 2017 (ODbL) 发表路线全链执行
# 用法: nohup bash src/stage0/run_2017_route.sh <50K仿真bash的PID> > data/synthetic_2017/route.log 2>&1 &
set -e
cd "$(dirname "$0")/../.."
PY=.venv/bin/python
SIM50K_PID="${1:-}"
mkdir -p data/synthetic_2017

echo "=== [$(date '+%F %T')] 2017 路线 Stage 2: CTGAN 生成 (1:1, 10,945) ==="
$PY -m src.stage2.run_generation \
    --input data/real/msk_impact_2017/real_dataset.csv \
    --output-dir data/synthetic_2017

echo "=== [$(date '+%F %T')] 2017 路线 Stage 3: 硬规则全量 + LLM 抽审 ==="
env -u all_proxy -u ALL_PROXY -u http_proxy -u https_proxy -u HTTP_PROXY -u HTTPS_PROXY \
    bash -c 'set -a && source .env && set +a && \
    .venv/bin/python src/stage0/clinical_auditor_v2.py \
        --input data/synthetic_2017/raw_synthetic_data.csv \
        --output-csv data/synthetic_2017/audited_synthetic_data.csv \
        --log-json data/synthetic_2017/audit_log.json'

echo "=== [$(date '+%F %T')] 2017 路线 Stage 4: 变异规格 + 仿真脚本编排 ==="
$PY src/stage0/run_stage4_v2.py \
    --input data/synthetic_2017/audited_synthetic_data.csv \
    --output-dir data/synthetic_2017 \
    --reads-root data/synthetic_2017/reads \
    --config configs/mutation_coordinates_2017.json \
    --fastq-n 2000

if [ -n "$SIM50K_PID" ]; then
    echo "=== [$(date '+%F %T')] 等待 50K 仿真 (PID $SIM50K_PID) 结束后再启动 2017 FASTQ ==="
    while kill -0 "$SIM50K_PID" 2>/dev/null; do sleep 120; done
fi

echo "=== [$(date '+%F %T')] 启动 2017 FASTQ 仿真 (2000 患者) ==="
bash data/synthetic_2017/run_simulation.sh
echo "=== [$(date '+%F %T')] 2017 路线全部完成 ==="

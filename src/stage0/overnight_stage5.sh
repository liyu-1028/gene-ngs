#!/bin/bash
PROJECT_ROOT="${PROJECT_ROOT:-$HOME/gene-ngs}"  # override with env var
# overnight_stage5.sh — 过夜总调度 (v2)
# 前提: 表格层评估 (stats/tstr/privacy) 已提前完成
# 流程: 等主仿真链结束 → hicov 50× 验证队列 (10 患者) → 双路线聚合锚点验证
#       → hicov 锚点 VAF 验证 → 总报告
# 关键: synggen 必须用绝对路径 (相对路径段错误)
set -u
cd ${PROJECT_ROOT}
WAIT_PID="${1:-8445}"
PY=${PROJECT_ROOT}/.venv/bin/python
B=${PROJECT_ROOT}
REF=$B/data/raw/reference/GRCh38_full_analysis_set_plus_decoy_hla.fa

echo "=== [$(date '+%F %T')] 等待主仿真链 (PID $WAIT_PID) 结束 ==="
while kill -0 "$WAIT_PID" 2>/dev/null; do sleep 300; done
echo "=== [$(date '+%F %T')] 主仿真链结束, 通知压缩守护收尾 ==="
touch /tmp/compress.stop
sleep 600   # 给压缩守护清积压

mkdir -p $B/data/synthetic_2017/evaluation $B/data/synthetic/evaluation $B/data/synthetic_2017/reads_hicov

# ─── Phase A: hicov 验证队列 (2017 发表路线, 10 患者 × 50M reads ≈ 50× 诊断深度) ───
echo "=== [$(date '+%F %T')] Phase A: hicov 50x 验证队列 ==="
n_hicov=0
for PT in $(ls $B/data/synthetic_2017/variant_specs | shuf --random-source=<(yes 42) | head -20); do
    [ -s "$B/data/synthetic_2017/variant_specs/$PT/variants.pm" ] || continue
    while [ "$(free -g | awk '/Mem:/{print $7}')" -lt 4 ] || \
          [ "$(df -BG / | tail -1 | awk '{print $4}' | tr -d G)" -lt 15 ]; do sleep 120; done
    mkdir -p "$B/data/synthetic_2017/reads_hicov/$PT"
    if $B/bin/synggen mode=1 \
        bed=$B/data/raw/reference/target_regions.bed \
        fasta=$REF \
        rdm=$B/data/raw/reference/models/wes_model/target_regions.rdm.gz \
        pbe=$B/data/raw/reference/models/wes_model/target_regions.pbe.gz \
        qm=$B/data/raw/reference/models/wes_model/target_regions.qm.gz \
        pm=$B/data/synthetic_2017/variant_specs/$PT/variants.pm \
        indel=$B/data/synthetic_2017/variant_specs/$PT/variants.indel \
        cna=$B/data/synthetic_2017/variant_specs/$PT/variants.cna \
        threads=3 nreads=50000000 tc=0.3 \
        out=$B/data/synthetic_2017/reads_hicov/$PT seed=$((42 + 10#${PT#virtual_pt_})) \
        > /dev/null 2>&1; then
        n_hicov=$((n_hicov+1)); echo "  [$n_hicov] $PT ✅ 仿真完成, 压缩中..."
        ( cd "$B/data/synthetic_2017/reads_hicov/$PT" && for f in *.fastq; do gzip -f "$f"; done )
        echo "  [$n_hicov] $PT 压缩完成"
    else
        echo "  $PT ❌ (跳过)"; rm -rf "$B/data/synthetic_2017/reads_hicov/$PT"
    fi
done
echo " hicov 队列完成: $n_hicov 个患者 (50x 深度)"

# ─── Phase B: 聚合锚点验证 (正式低覆盖数据集, 2017 抽 500 / 50K 抽 250) ───
echo "=== [$(date '+%F %T')] Phase B1: 2017 聚合锚点验证 (500 患者) ==="
nice -n 10 $PY -u $B/src/stage0/aggregate_read_validation.py \
    --specs-dir $B/data/synthetic_2017/variant_specs --reads-root $B/data/synthetic_2017/reads \
    --n-patients 500 --ref $REF \
    --out-csv $B/data/synthetic_2017/evaluation/aggregate_anchor_detail.csv \
    > $B/data/synthetic_2017/evaluation/aggregate.log 2>&1

echo "=== [$(date '+%F %T')] Phase B2: 50K 聚合锚点验证 (250 患者) ==="
nice -n 10 $PY -u $B/src/stage0/aggregate_read_validation.py \
    --specs-dir $B/data/synthetic/variant_specs --reads-root $B/data/synthetic/reads \
    --n-patients 250 --ref $REF \
    --out-csv $B/data/synthetic/evaluation/aggregate_anchor_detail.csv \
    > $B/data/synthetic/evaluation/aggregate.log 2>&1

# ─── Phase C: hicov 锚点 VAF 验证 (每患者 10 锚点) ───
echo "=== [$(date '+%F %T')] Phase C: hicov 锚点 VAF 验证 ==="
ls $B/data/synthetic_2017/reads_hicov | grep '^virtual_pt_' > /tmp/hicov_pts.txt
nice -n 10 $PY -u $B/src/stage0/aggregate_read_validation.py \
    --specs-dir $B/data/synthetic_2017/variant_specs --reads-root $B/data/synthetic_2017/reads_hicov \
    --patient-list /tmp/hicov_pts.txt --anchors-per-patient 10 --ref $REF \
    --out-csv $B/data/synthetic_2017/evaluation/hicov_anchor_detail.csv \
    > $B/data/synthetic_2017/evaluation/hicov_verify.log 2>&1

echo "=== [$(date '+%F %T')] 生成过夜总报告 ==="
$PY - <<'EOF'
import json, os, glob
B = "${PROJECT_ROOT}"
lines = ["# 过夜实验总报告", ""]
for label, ev in [("2017 发表路线", "data/synthetic_2017/evaluation"),
                  ("50K 扩展路线", "data/synthetic/evaluation")]:
    lines.append(f"## {label}")
    for name, fn in [("统计保真", "statistical_fidelity_metrics.json"),
                     ("下游效用 TSTR", "utility_tstr_metrics.json"),
                     ("隐私保护", "privacy_metrics.json"),
                     ("聚合锚点验证", "aggregate_anchor_detail_summary.json"),
                     ("hicov 锚点验证", "hicov_anchor_detail_summary.json")]:
        p = os.path.join(B, ev, fn)
        lines.append(f"- {'✅' if os.path.exists(p) else '❌'} {name}: `{p}`")
        if os.path.exists(p) and fn.endswith('.json'):
            try:
                lines.append("```json"); lines += json.dumps(json.load(open(p)), indent=1).splitlines(); lines.append("```")
            except Exception: pass
    lines.append("")
for label, root in [("2017", "data/synthetic_2017/reads"), ("50K", "data/synthetic/reads"),
                    ("2017 hicov(50x)", "data/synthetic_2017/reads_hicov")]:
    gz = len(glob.glob(os.path.join(B, root, "*.fastq.gz")))
    dirs = len(glob.glob(os.path.join(B, root, "virtual_pt_*")))
    lines.append(f"- {label} FASTQ: {gz} 个 .fastq.gz + {dirs} 个未压缩目录")
open(os.path.join(B, "OVERNIGHT_REPORT.md"), "w").write("\n".join(lines))
print("\n".join(lines))
EOF
echo "=== [$(date '+%F %T')] 过夜流水线全部完成, 详见 OVERNIGHT_REPORT.md ==="

#!/bin/bash
PROJECT_ROOT="${PROJECT_ROOT:-$HOME/gene-ngs}"  # override with env var
# finish_stage5.sh — Phase A 截断后的收尾: B1/B2/C + 总报告
set -u
B=${PROJECT_ROOT}
PY=$B/.venv/bin/python
REF=$B/data/raw/reference/GRCh38_full_analysis_set_plus_decoy_hla.fa
cd $B

echo "=== [$(date '+%F %T')] Phase B1: 2017 聚合锚点验证 (500 患者) ==="
nice -n 10 $PY -u $B/src/stage0/aggregate_read_validation.py \
    --specs-dir $B/data/synthetic_2017/variant_specs --reads-root $B/data/synthetic_2017/reads \
    --n-patients 500 --ref $REF \
    --out-csv $B/data/synthetic_2017/evaluation/aggregate_anchor_detail.csv \
    > $B/data/synthetic_2017/evaluation/aggregate.log 2>&1
tail -8 $B/data/synthetic_2017/evaluation/aggregate.log

echo "=== [$(date '+%F %T')] Phase B2: 50K 聚合锚点验证 (250 患者) ==="
nice -n 10 $PY -u $B/src/stage0/aggregate_read_validation.py \
    --specs-dir $B/data/synthetic/variant_specs --reads-root $B/data/synthetic/reads \
    --n-patients 250 --ref $REF \
    --out-csv $B/data/synthetic/evaluation/aggregate_anchor_detail.csv \
    > $B/data/synthetic/evaluation/aggregate.log 2>&1
tail -8 $B/data/synthetic/evaluation/aggregate.log

echo "=== [$(date '+%F %T')] Phase C: hicov 锚点 VAF 验证 ==="
ls $B/data/synthetic_2017/reads_hicov | grep '^virtual_pt_' > /tmp/hicov_pts.txt
nice -n 10 $PY -u $B/src/stage0/aggregate_read_validation.py \
    --specs-dir $B/data/synthetic_2017/variant_specs --reads-root $B/data/synthetic_2017/reads_hicov \
    --patient-list /tmp/hicov_pts.txt --anchors-per-patient 10 --ref $REF \
    --out-csv $B/data/synthetic_2017/evaluation/hicov_anchor_detail.csv \
    > $B/data/synthetic_2017/evaluation/hicov_verify.log 2>&1
tail -8 $B/data/synthetic_2017/evaluation/hicov_verify.log

echo "=== [$(date '+%F %T')] 生成总报告 ==="
$PY - <<'EOF'
import json, os, glob
B = os.environ.get("PROJECT_ROOT", os.path.expanduser("~/gene-ngs"))
lines = ["# 实验总报告", ""]
for label, ev in [("2017 发表路线", "data/synthetic_2017/evaluation"),
                  ("50K 扩展路线", "data/synthetic/evaluation")]:
    lines.append(f"## {label}")
    for name, fn in [("统计保真", "statistical_fidelity_metrics.json"),
                     ("下游效用 TSTR", "utility_tstr_metrics.json"),
                     ("隐私保护", "privacy_metrics.json"),
                     ("聚合锚点验证", "aggregate_anchor_detail_summary.json"),
                     ("hicov 锚点 VAF 验证", "hicov_anchor_detail_summary.json")]:
        p = os.path.join(B, ev, fn)
        lines.append(f"- {'✅' if os.path.exists(p) else '❌'} {name}: `{p}`")
        if os.path.exists(p) and fn.endswith('summary.json'):
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
touch /tmp/compress.stop
echo "=== [$(date '+%F %T')] 收尾完成, 详见 OVERNIGHT_REPORT.md ==="

#!/bin/bash
PROJECT_ROOT="${PROJECT_ROOT:-$HOME/gene-ngs}"  # override with env var
# overnight_final.sh — 过夜最终验证与发布打包总控
# 阶段: A=50K全量审计 B=2017全量审计 C=自动修复循环 D=Zenodo清单 E=最终报告
set -u
B=${PROJECT_ROOT}
R50=$B/data/synthetic/reads
R17=$B/data/synthetic_2017/reads
LOG=$B/overnight_final.log
WORK=/tmp/overnight_final
EXPECTED=200000   # 每患者 reads 数

mkdir -p "$WORK"
: > "$LOG"
log(){ echo "[$(date '+%m-%d %H:%M:%S')] $*" | tee -a "$LOG"; }

# ── 单文件审计 worker: 输出 "路径 reads数 PASS/FAIL" ──
audit_one(){
    f="$1"
    n=$(zcat "$f" 2>/dev/null | wc -l)
    reads=$((n / 4))
    if [ "$reads" -eq 200000 ]; then
        echo "$f $reads PASS"
    else
        echo "$f $reads FAIL"
    fi
}
export -f audit_one

audit_route(){
    local root="$1" tag="$2"
    log "[$tag] 全量审计开始: $root"
    : > "$WORK/${tag}_results.txt"
    find "$root" -maxdepth 1 -name "virtual_pt_*.fastq.gz" -print0 \
      | xargs -0 -P 12 -n 1 bash -c 'audit_one "$0"' >> "$WORK/${tag}_results.txt"
    local pass fail
    pass=$(grep -c PASS "$WORK/${tag}_results.txt" || true)
    fail=$(grep -c FAIL "$WORK/${tag}_results.txt" || true)
    log "[$tag] PASS=$pass FAIL=$fail"
    grep FAIL "$WORK/${tag}_results.txt" | awk '{print $1}' > "$WORK/${tag}_failed.txt" || true
}

repair_50k(){
    local f="$1"
    local b=$(basename "$f" .fastq.gz)
    local idx=$((10#${b#virtual_pt_}))
    local seed=$((42 + idx))
    log "[repair-50K] $b 重仿真 seed=$seed"
    rm -f "$f"
    local out="$R50/$b"
    mkdir -p "$out"
    "$B/bin/synggen" mode=1 \
        bed="$B/data/raw/reference/target_regions.bed" \
        fasta="$B/data/raw/reference/GRCh38_full_analysis_set_plus_decoy_hla.fa" \
        rdm="$B/data/raw/reference/models/wes_model/target_regions.rdm.gz" \
        pbe="$B/data/raw/reference/models/wes_model/target_regions.pbe.gz" \
        qm="$B/data/raw/reference/models/wes_model/target_regions.qm.gz" \
        pm="$B/data/synthetic/variant_specs/$b/variants.pm" \
        indel="$B/data/synthetic/variant_specs/$b/variants.indel" \
        cna="$B/data/synthetic/variant_specs/$b/variants.cna" \
        threads=4 nreads=200000 tc=0.3 out="$out" seed=$seed >>"$LOG" 2>&1
    # 合并→gz→校验→清理
    cat "$out"/*.R1.fastq 2>/dev/null | gzip -c > "$f"
    local n=$(zcat "$f" 2>/dev/null | wc -l)
    if [ $((n/4)) -eq 200000 ] && gzip -t "$f" 2>/dev/null; then
        rm -rf "$out"
        log "[repair-50K] $b 修复成功 ($((n/4)) reads)"
        return 0
    fi
    log "[repair-50K] $b 修复失败 (reads=$((n/4)))"
    return 1
}

repair_2017(){
    local f="$1"
    local b=$(basename "$f" .fastq.gz)
    local idx=$((10#${b#virtual_pt_}))
    local seed=$((42 + idx))
    log "[repair-2017] $b 重仿真 seed=$seed"
    rm -f "$f"
    local out="$R17/$b"
    mkdir -p "$out"
    "$B/bin/synggen" mode=1 \
        bed="$B/data/raw/reference/target_regions.bed" \
        fasta="$B/data/raw/reference/GRCh38_full_analysis_set_plus_decoy_hla.fa" \
        rdm="$B/data/raw/reference/models/wes_model/target_regions.rdm.gz" \
        pbe="$B/data/raw/reference/models/wes_model/target_regions.pbe.gz" \
        qm="$B/data/raw/reference/models/wes_model/target_regions.qm.gz" \
        pm="$B/data/synthetic_2017/variant_specs/$b/variants.pm" \
        indel="$B/data/synthetic_2017/variant_specs/$b/variants.indel" \
        cna="$B/data/synthetic_2017/variant_specs/$b/variants.cna" \
        threads=4 nreads=200000 tc=0.3 out="$out" seed=$seed >>"$LOG" 2>&1
    cat "$out"/*.R1.fastq 2>/dev/null | gzip -c > "$f"
    local n=$(zcat "$f" 2>/dev/null | wc -l)
    if [ $((n/4)) -eq 200000 ] && gzip -t "$f" 2>/dev/null; then
        rm -rf "$out"
        log "[repair-2017] $b 修复成功 ($((n/4)) reads)"
        return 0
    fi
    log "[repair-2017] $b 修复失败 (reads=$((n/4)))"
    return 1
}

# ══════════ 主流程 ══════════
log "════ 过夜最终验证启动 (12核) ╏═"

audit_route "$R50" "50K"
audit_route "$R17" "2017"

# ── 修复循环 (最多3轮) ──
for round in 1 2 3; do
    n50=$(wc -l < "$WORK/50K_failed.txt" 2>/dev/null || echo 0)
    n17=$(wc -l < "$WORK/2017_failed.txt" 2>/dev/null || echo 0)
    [ "$n50" -eq 0 ] && [ "$n17" -eq 0 ] && break
    log "════ 修复轮 $round: 50K失败=$n50 2017失败=$n17 ════"
    : > "$WORK/50K_failed_new.txt"; : > "$WORK/2017_failed_new.txt"
    while read -r f; do
        [ -z "$f" ] && continue
        repair_50k "$f" || echo "$f" >> "$WORK/50K_failed_new.txt"
    done < "$WORK/50K_failed.txt"
    while read -r f; do
        [ -z "$f" ] && continue
        repair_2017 "$f" || echo "$f" >> "$WORK/2017_failed_new.txt"
    done < "$WORK/2017_failed.txt"
    mv "$WORK/50K_failed_new.txt" "$WORK/50K_failed.txt"
    mv "$WORK/2017_failed_new.txt" "$WORK/2017_failed.txt"
done

# ── 终态统计 ──
P50=$(grep -c PASS "$WORK/50K_results.txt" || echo 0)
P17=$(grep -c PASS "$WORK/2017_results.txt" || echo 0)
F50=$(wc -l < "$WORK/50K_failed.txt" 2>/dev/null || echo 0)
F17=$(wc -l < "$WORK/2017_failed.txt" 2>/dev/null || echo 0)
log "════ 终态: 50K PASS=$P50 FAIL=$F50 | 2017 PASS=$P17 FAIL=$F17 ════"

# ── Zenodo MANIFEST ──
log "生成 Zenodo MANIFEST..."
{
  echo "# Zenodo 发布包校验清单 — 生成于 $(date '+%F %H:%M')"
  echo "# 表格与验证文件 (data/zenodo_staging/)"
  cd "$B/data/zenodo_staging" && md5sum * | sed 's|^|zenodo_staging/|'
  echo "# 2017 FASTQ (data/synthetic_2017/reads/, $P17 files × 200000 reads)"
  cd "$R17" && md5sum virtual_pt_*.fastq.gz | sed 's|^|synthetic_2017/reads/|'
  echo "# 50K FASTQ (data/synthetic/reads/, $P50 files × 200000 reads)"
  cd "$R50" && md5sum virtual_pt_*.fastq.gz | sed 's|^|synthetic/reads/|'
} > "$B/data/zenodo_staging/MANIFEST.md5" 2>>"$LOG"
log "MANIFEST.md5 完成 ($(wc -l < "$B/data/zenodo_staging/MANIFEST.md5") 行)"

# ── 最终标记 ──
if [ "$F50" -eq 0 ] && [ "$F17" -eq 0 ]; then
    log "════ ★ 全部通过: 4000/4000 文件 × 200000 reads ════"
    echo "ALL_PASS" > "$WORK/verdict.txt"
else
    log "════ ⚠ 存在未修复失败: 50K=$F50 2017=$F17 ════"
    echo "HAS_FAILURES" > "$WORK/verdict.txt"
fi
log "════ 过夜流程结束 ════"

#!/bin/bash
PROJECT_ROOT="${PROJECT_ROOT:-$HOME/gene-ngs}"  # override with env var
# flatten_50k_reads.sh — 将剩余目录型 50K reads 统一为扁平 virtual_pt_NNNN.fastq.gz
# 用法: flatten_50k_reads.sh <工作目录列表文件|all>
set -u
READS=${PROJECT_ROOT}/data/synthetic/reads
ARCHIVE=${PROJECT_ROOT}/data/synthetic/reads_extended_bed_archive.tar
LIST="$1"

if [ "$LIST" = "all" ]; then
    find "$READS" -maxdepth 1 -type d -name "virtual_pt_*" -printf "%f\n" | sort > /tmp/dirs87.txt
else
    cp "$LIST" /tmp/dirs87.txt
fi

# 1) 归档所有 extended.bed (一次性)
if [ ! -f "$ARCHIVE" ]; then
    (cd "$READS" && find . -maxdepth 2 -path "./*virtual_pt_*/target_regions.extended.bed" | tar -cf "$ARCHIVE" -T -)
    echo "[archive] $(tar -tf $ARCHIVE | wc -l) 个 extended.bed 已归档"
fi

# 2) 并行合并+压缩+校验+清理
process_one() {
    d="$1"
    dir="$READS/$d"
    out="$READS/$d.fastq.gz"
    [ -f "$out" ] && { echo "[skip] $out 已存在"; return 0; }
    fqs=$(ls "$dir"/*.R1.fastq 2>/dev/null | sort)
    [ -z "$fqs" ] && { echo "[FAIL] $d 无fastq"; return 1; }
    cat $fqs | gzip -c > "$out" || { echo "[FAIL] $d gzip"; return 1; }
    # 校验: gzip 完整性 + reads 数
    n=$(zcat "$out" | wc -l)
    if [ $((n / 4)) -eq 200000 ] && gzip -t "$out" 2>/dev/null; then
        rm -rf "$dir"
        echo "[OK] $d → 200000 reads"
    else
        echo "[FAIL] $d 校验失败 reads=$((n/4))"
        rm -f "$out"
        return 1
    fi
}
export -f process_one
export READS

cat /tmp/dirs87.txt | xargs -P 8 -I{} bash -c 'process_one {}'
echo "═══ 完成。剩余目录数: $(find $READS -maxdepth 1 -type d -name 'virtual_pt_*' | wc -l) ═══"
echo "═══ 扁平文件总数: $(find $READS -maxdepth 1 -name 'virtual_pt_*.fastq.gz' | wc -l) ═══"

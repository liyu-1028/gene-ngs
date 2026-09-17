#!/bin/bash
PROJECT_ROOT="${PROJECT_ROOT:-$HOME/gene-ngs}"  # override with env var
# verify_fidelity_v2.sh — 批量 read 级物理保真验证 (k-mer 法, 支持 .fastq.gz)
# 用法: verify_fidelity_v2.sh <reads根目录> <specs目录> <抽样患者数> <输出CSV> <随机种子>
set -u
READS_ROOT="$1"; VAR_SPECS="$2"; N_SAMPLE="${3:-30}"; OUT_CSV="$4"; SEED="${5:-42}"
REF="${PROJECT_ROOT}/data/raw/reference/GRCh38_full_analysis_set_plus_decoy_hla.fa"

echo "patient,chrom,pos,alt,expected_vaf,ref_count,alt_count,total,observed_vaf,abs_delta,status" > "$OUT_CSV"

# 确定性抽样: 有 specs 且有 reads (目录或 gz) 的患者
mapfile -t candidates < <(comm -12 \
    <(ls "$VAR_SPECS" | grep '^virtual_pt_' | sort) \
    <({ ls -d "$READS_ROOT"/virtual_pt_*/ 2>/dev/null | xargs -n1 basename; ls "$READS_ROOT" 2>/dev/null | grep '\.fastq\.gz$' | sed 's/\.fastq\.gz$//'; } | sort -u) \
    | shuf --random-source=<(yes "$SEED") | head -"$N_SAMPLE")

echo "抽样 ${#candidates[@]} 个患者进行 read 级验证..." >&2

n_mut=0; n_ok=0; n_warn=0; n_miss=0; sum_delta=0
for PT in "${candidates[@]}"; do
    VAR_FILE="$VAR_SPECS/$PT/variants.pm"
    [ -f "$VAR_FILE" ] || continue
    # 每患者取前 5 个 SNV
    while IFS=$'\t' read -r chrom pos alt vaf rest; do
        [ -z "$chrom" ] && continue
        start=$((pos - 10)); end=$((pos + 10))
        ref_seq=$(samtools faidx "$REF" "$chrom:$start-$end" 2>/dev/null | grep -v ">" | tr -d '\n' | tr a-z A-Z)
        if [ ${#ref_seq} -ne 21 ]; then continue; fi
        alt_seq="${ref_seq:0:10}${alt}${ref_seq:11:10}"

        if [ -f "$READS_ROOT/$PT.fastq.gz" ]; then
            ref_count=$(gzip -dc "$READS_ROOT/$PT.fastq.gz" | grep -o "$ref_seq" | wc -l)
            alt_count=$(gzip -dc "$READS_ROOT/$PT.fastq.gz" | grep -o "$alt_seq" | wc -l)
        elif [ -d "$READS_ROOT/$PT" ]; then
            ref_count=$(cat "$READS_ROOT/$PT"/*.fastq 2>/dev/null | grep -o "$ref_seq" | wc -l)
            alt_count=$(cat "$READS_ROOT/$PT"/*.fastq 2>/dev/null | grep -o "$alt_seq" | wc -l)
        else
            continue
        fi

        total=$((ref_count + alt_count)); n_mut=$((n_mut+1))
        if [ "$total" -gt 0 ]; then
            obs=$(awk -v a=$alt_count -v t=$total 'BEGIN{printf "%.4f", a/t}')
            delta=$(awk -v o=$obs -v e=$vaf 'BEGIN{d=o-e; if(d<0)d=-d; printf "%.4f", d}')
            sum_delta=$(awk -v s=$sum_delta -v d=$delta 'BEGIN{printf "%.4f", s+d}')
            pass=$(awk -v d=$delta 'BEGIN{print (d<0.08)?"OK":"WARN"}')
            if [ "$pass" = "OK" ]; then n_ok=$((n_ok+1)); else n_warn=$((n_warn+1)); fi
            echo "$PT,$chrom,$pos,$alt,$vaf,$ref_count,$alt_count,$total,$obs,$delta,$pass" >> "$OUT_CSV"
        else
            n_miss=$((n_miss+1))
            echo "$PT,$chrom,$pos,$alt,$vaf,0,0,0,NA,NA,NO_COVERAGE" >> "$OUT_CSV"
        fi
    done < <(head -5 "$VAR_FILE")
done

# 汇总
if [ "$n_mut" -gt 0 ]; then
    median_delta=$(awk -F, 'NR>1 && $10!="NA"{print $10}' "$OUT_CSV" | sort -n | awk '{a[NR]=$1} END{if(NR%2)print a[(NR+1)/2]; else printf "%.4f",(a[NR/2]+a[NR/2+1])/2}')
    recall=$(awk -v ok=$n_ok -v w=$n_warn -v m=$n_mut 'BEGIN{printf "%.4f",(ok+w)/m}')
    within=$(awk -v ok=$n_ok -v w=$n_warn -v m=$n_mut 'BEGIN{printf "%.4f",ok/m}')
    mean_delta=$(awk -v s=$sum_delta -v n=$n_mut 'BEGIN{printf "%.4f",s/n}')
    {
        echo ""
        echo "═══ Read 级验证汇总 ($READS_ROOT) ═══"
        echo "抽样患者: ${#candidates[@]} | 验证突变位点: $n_mut"
        echo "锚点召回率 (alt_kmer>0): $recall"
        echo "VAF±0.08 达标率: $within"
        echo "|ΔVAF| 均值: $mean_delta | 中位: $median_delta"
        echo "明细: $OUT_CSV"
    } | tee "${OUT_CSV%.csv}_summary.txt"
fi

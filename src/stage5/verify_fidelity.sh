#!/bin/bash

# ==============================================================================
# Stage 5: Physical Fidelity Verification (Fast K-mer Method)
# This script directly verifies if the simulated FASTQ files physically contain
# the injected mutations at the specified Variant Allele Frequencies (VAF).
# It bypasses the 2-hour 'bwa index' step by searching for sequence contexts.
# ==============================================================================

PT_ID=${1:-"virtual_pt_0000"}
# Dynamically resolve workspace root (two levels up from src/stage5/)
SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
WORKSPACE_DIR="$( cd "$SCRIPT_DIR/../.." && pwd )"
REF="$WORKSPACE_DIR/data/raw/reference/GRCh38_full_analysis_set_plus_decoy_hla.fa"
READS_DIR="$WORKSPACE_DIR/data/synthetic/reads/$PT_ID"
VAR_FILE="$WORKSPACE_DIR/data/synthetic/variant_specs/$PT_ID/variants.pm"

echo "=========================================================="
echo " Stage 5: Physical Read-Level Fidelity Verification"
echo " Patient: $PT_ID"
echo "=========================================================="

if [ ! -f "$VAR_FILE" ]; then
    echo "Error: Variant file not found at $VAR_FILE"
    exit 1
fi

if [ ! -d "$READS_DIR" ]; then
    echo "Error: Reads directory not found at $READS_DIR"
    exit 1
fi

echo "Injected Point Mutations from variants.pm:"
cat $VAR_FILE
echo "----------------------------------------------------------"

# Iterate over each point mutation
while IFS=$'\t' read -r chrom pos alt vaf clonality cnAllele; do
    echo ">> Checking Locus: $chrom:$pos (Alt: $alt, Expected VAF: $vaf)"
    
    # We will extract a 21-bp context (10bp before + mutated base + 10bp after)
    start=$((pos - 10))
    end=$((pos + 10))
    
    # 1. Extract reference context using samtools
    ref_seq=$(samtools faidx $REF $chrom:$start-$end | grep -v ">" | tr -d '\n' | tr a-z A-Z)
    
    if [ ${#ref_seq} -ne 21 ]; then
        echo "   [!] Could not fetch reference sequence properly. Skipping."
        continue
    fi
    
    # 2. Build the mutated sequence context
    ref_base=${ref_seq:10:1}
    alt_seq="${ref_seq:0:10}${alt}${ref_seq:11:10}"
    
    echo "   Ref Context: $ref_seq (Center Base: $ref_base)"
    echo "   Alt Context: $alt_seq (Center Base: $alt)"
    
    # 3. Count exact K-mer matches in the FASTQ files (supporting both compressed and uncompressed)
    if ls $READS_DIR/*.gz >/dev/null 2>&1; then
        ref_count=$(gzip -dc $READS_DIR/*.gz | grep -o "$ref_seq" | wc -l)
        alt_count=$(gzip -dc $READS_DIR/*.gz | grep -o "$alt_seq" | wc -l)
    elif ls $READS_DIR/*.fastq >/dev/null 2>&1; then
        ref_count=$(cat $READS_DIR/*.fastq | grep -o "$ref_seq" | wc -l)
        alt_count=$(cat $READS_DIR/*.fastq | grep -o "$alt_seq" | wc -l)
    elif ls $READS_DIR/*.fq >/dev/null 2>&1; then
        ref_count=$(cat $READS_DIR/*.fq | grep -o "$ref_seq" | wc -l)
        alt_count=$(cat $READS_DIR/*.fq | grep -o "$alt_seq" | wc -l)
    else
        echo "   [!] No FASTQ or gzipped FASTQ files found in $READS_DIR"
        ref_count=0
        alt_count=0
    fi
    
    # 4. Calculate observed VAF
    total=$((ref_count + alt_count))
    if [ "$total" -gt 0 ]; then
        obs_vaf=$(echo "scale=4; $alt_count / $total" | bc)
        echo "   Result: Ref Reads = $ref_count | Alt Reads = $alt_count"
        echo "   => Observed VAF : $obs_vaf (Expected: $vaf)"
        
        # Simple check for success (tolerance +/- 0.08)
        diff=$(echo "$obs_vaf - $vaf" | bc | awk '{print sqrt($1*$1)}')
        if (( $(echo "$diff < 0.08" | bc -l) )); then
            echo "   [SUCCESS] Mutation injected physically at high fidelity! [OK]"
        else
            echo "   [WARNING] Mutation observed, but VAF differs slightly. [WARN]"
        fi
    else
        echo "   Result: No exact 21-bp context matches found. (Coverage may be low or insert size varied)"
    fi
    echo "----------------------------------------------------------"
done < "$VAR_FILE"

echo ""
echo "=========================================================="
echo " Optional: Full BWA-MEM Alignment Guide"
echo "=========================================================="
echo "To perform a standard bioinformatics validation, run:"
echo "1. bwa index $REF  # Takes ~2 hours"
echo "2. cat $READS_DIR/*.fastq > $READS_DIR/merged.fq"
echo "3. bwa mem -t 4 $REF $READS_DIR/merged.fq | samtools sort -o $READS_DIR/mapped.bam"
echo "4. samtools index $READS_DIR/mapped.bam"
echo "5. samtools mpileup -r <chrom>:<pos>-<pos> $READS_DIR/mapped.bam"

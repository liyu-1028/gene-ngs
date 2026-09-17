#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
build_mutation_coordinate_db.py  (Stage 0 附属: stage4 坐标映射库重建)

为 MSK-IMPACT 50K 数据源重建 stage4/Synggen 的变异坐标映射库
(configs/mutation_coordinates_50k.json)，替代旧版仅 12 SNV + 7 CNA 的小库。

构建逻辑:
  SNV 库:
    - 源: 50K data_mutations.txt (GRCh37)
    - 仅 Variant_Type=='SNP' 且 alt 为单碱基 (A/C/G/T)
    - key = "Hugo_Symbol:HGVSp_Short" (缺则 HGVSc)
    - 排除 key 含 del/ins/dup 子串的 (路由到 stage4 的 indel 分支)
    - 同 key 多坐标取众数
    - 坐标 hg19 -> hg38 liftover (UCSC chain)，失败的 key 丢弃
  CNA 库:
    - 源: GRCh38 target_regions.bed (第4列基因名)，对 data_cna.txt 的 541 个
      panel 基因取外显子跨度和 (min start, max end)

用法:
    python src/stage0/build_mutation_coordinate_db.py
"""

import argparse
import json
import os
import sys
from collections import Counter

import pandas as pd
from pyliftover import LiftOver

BASE = os.path.abspath(os.path.join(os.path.dirname(__file__), '../..'))
BED = os.path.join(BASE, 'data/raw/reference/target_regions.bed')
CHAIN = os.path.join(BASE, 'data/raw/reference/liftover/hg19ToHg38.over.chain.gz')

REPORTABLE_VC = {
    'Missense_Mutation', 'Nonsense_Mutation', 'Frame_Shift_Del', 'Frame_Shift_Ins',
    'In_Frame_Del', 'In_Frame_Ins', 'Splice_Site', 'Translation_Start_Site',
    'Nonstop_Mutation',
}


def build_snv_db(liftover: LiftOver, maf_path: str, skip_header_rows: int) -> dict:
    print(f"[1/3] 读取 MAF ({maf_path}) 并筛选可注入 SNP ...")
    mut = pd.read_csv(maf_path, sep='\t', skiprows=skip_header_rows, low_memory=False,
                      usecols=['Hugo_Symbol', 'Chromosome', 'Start_Position',
                               'Variant_Classification', 'Variant_Type',
                               'Reference_Allele', 'Tumor_Seq_Allele2',
                               'HGVSp_Short', 'HGVSc'])
    mut = mut[mut['Variant_Classification'].isin(REPORTABLE_VC)]
    mut = mut[mut['Variant_Type'] == 'SNP']
    mut = mut[mut['Tumor_Seq_Allele2'].isin(['A', 'C', 'G', 'T'])]

    # key 构造: 优先 HGVSp_Short，缺则 HGVSc；排除路由歧义 key
    def make_key(r):
        for f in ('HGVSp_Short', 'HGVSc'):
            v = r[f]
            if isinstance(v, str) and v and v != '.':
                k = f"{r['Hugo_Symbol']}:{v}"
                if not any(w in k for w in ('del', 'ins', 'dup')):
                    return k
        return None

    mut = mut[mut['Hugo_Symbol'].notna()]
    mut['key'] = mut.apply(make_key, axis=1)
    mut = mut[mut['key'].notna()]
    print(f"      候选 key: {mut['key'].nunique():,}")

    # 同 key 取众数坐标
    coord_counter = {}
    for (k, chrom, pos, alt), n in (
            mut.groupby(['key', 'Chromosome', 'Start_Position', 'Tumor_Seq_Allele2'])
            .size().items()):
        coord_counter.setdefault(k, Counter())[(str(chrom), int(float(pos)), alt)] += n

    print("[2/3] hg19 -> hg38 liftover ...")
    snv_db, failed = {}, 0
    chrom_fix = {str(i): f'chr{i}' for i in list(range(1, 23)) + ['X', 'Y']}
    for k, counter in coord_counter.items():
        (chrom, pos, alt), _ = counter.most_common(1)[0]
        res = liftover.convert_coordinate(chrom_fix.get(chrom, chrom), pos - 1)
        if res:
            out_chrom, out_pos0, _, _ = res[0]
            snv_db[k] = [out_chrom, out_pos0 + 1, alt]
        else:
            failed += 1
    print(f"      SNV 库: {len(snv_db):,} | liftover 失败丢弃: {failed:,}")
    return snv_db


def build_cna_db(cna_path: str) -> dict:
    print("[3/3] 从 GRCh38 BED 构建 CNA 基因跨度库 ...")
    panel_genes = set(pd.read_csv(cna_path, sep='\t', index_col=0).index)
    spans, chrom_of = {}, {}
    with open(BED) as f:
        for line in f:
            p = line.rstrip('\n').split('\t')
            if len(p) < 4:
                continue
            gene = p[3]
            if gene not in panel_genes:
                continue
            s, e = int(p[1]), int(p[2])
            lo, hi = spans.get(gene, (s, e))
            spans[gene] = (min(lo, s), max(hi, e))
            chrom_of.setdefault(gene, p[0])
    cna_db = {g: [chrom_of[g], lo, hi] for g, (lo, hi) in spans.items()}
    print(f"      CNA 库: {len(cna_db):,} / {len(panel_genes)} 个 panel 基因命中")
    return cna_db


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--maf', default=os.path.join(BASE, 'data/msk_impact_50k_2026/data_mutations.txt'))
    ap.add_argument('--cna', default=os.path.join(BASE, 'data/msk_impact_50k_2026/data_cna.txt'))
    ap.add_argument('--out', default=os.path.join(BASE, 'configs/mutation_coordinates_50k.json'))
    ap.add_argument('--maf-skiprows', type=int, default=0,
                    help='MAF 表头前的注释行数: 50K 版=0, 2017 版=2')
    args = ap.parse_args()

    liftover = LiftOver(CHAIN)
    snv_db = build_snv_db(liftover, args.maf, args.maf_skiprows)
    cna_db = build_cna_db(args.cna)
    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    with open(args.out, 'w') as f:
        json.dump({'SNV': snv_db, 'CNA': cna_db}, f)
    size_mb = os.path.getsize(args.out) / 1e6
    print(f"\n✓ 写出 {args.out} ({size_mb:.1f} MB, SNV {len(snv_db):,} + CNA {len(cna_db):,})")
    # 抽样验证
    for k in ['EGFR:p.L858R', 'KRAS:p.G12C', 'BRAF:p.V600E', 'TP53:p.R248Q']:
        print(f"  {k} -> {snv_db.get(k, 'MISS')}")


if __name__ == '__main__':
    sys.exit(main())

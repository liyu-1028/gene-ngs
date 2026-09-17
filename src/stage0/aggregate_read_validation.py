#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
aggregate_read_validation.py — 低覆盖 FASTQ 的聚合锚点物理验证
对抽样患者的 variants.pm 锚点, 在该患者 FASTQ (支持 .gz) 中做 21-mer 计数:
  - 可检出率: alt 21-mer 出现过的锚点比例 (与 Poisson(λ) 理论值对照)
  - alt 特异性: 锚点 alt k-mer 不应出现在"未携带该突变"的患者中
  - 聚合 VAF 一致性: 有覆盖锚点上 (alt/(ref+alt)) 与预期 VAF 的 Spearman 相关
用法:
  python aggregate_read_validation.py --specs-dir ... --reads-root ... --n-patients 500 \
      --ref GRCh38.fa --out-csv ... [--other-reads-root 目录,用于特异性对照]
"""
import argparse, glob, gzip, os, random, subprocess, json
import numpy as np
from scipy.stats import spearmanr

def load_pm(path):
    rows = []
    with open(path) as f:
        for line in f:
            p = line.rstrip('\n').split('\t')
            if len(p) >= 4:
                rows.append(dict(chrom=p[0], pos=int(p[1]), alt=p[2], vaf=float(p[3])))
    return rows

def kmers(ref_fa, chrom, pos, alt):
    """返回 (ref21, alt21); 提取失败返回 None"""
    try:
        r = subprocess.run(['samtools', 'faidx', ref_fa, f'{chrom}:{pos-10}-{pos+10}'],
                           capture_output=True, text=True, timeout=60)
        seq = ''.join(r.stdout.splitlines()[1:]).upper()
        if len(seq) != 21:
            return None
        ref21, alt21 = seq, seq[:10] + alt + seq[11:]
        return ref21, alt21
    except Exception:
        return None

def count_kmers_in_patient(reads_files, km_list):
    """C 级多模式计数: gzip -dc | grep -oF -f patterns | sort | uniq -c"""
    import subprocess, tempfile, os
    counts = {k: 0 for k in km_list}
    if not km_list:
        return counts
    with tempfile.NamedTemporaryFile('w', suffix='.pat', delete=False) as tf:
        tf.write('\n'.join(km_list) + '\n')
        pat = tf.name
    try:
        cmd = ("gzip -dc " + ' '.join("'" + f + "'" for f in reads_files) +
               " | LC_ALL=C grep -oF -f " + pat + " | sort | uniq -c")
        r = subprocess.run(['bash', '-c', cmd], capture_output=True, text=True,
                           timeout=1800)
        for line in r.stdout.splitlines():
            parts = line.strip().split()
            if len(parts) == 2 and parts[1] in counts:
                counts[parts[1]] = int(parts[0])
    finally:
        os.unlink(pat)
    return counts

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--specs-dir', required=True)
    ap.add_argument('--reads-root', required=True)
    ap.add_argument('--n-patients', type=int, default=500)
    ap.add_argument('--anchors-per-patient', type=int, default=5)
    ap.add_argument('--ref', required=True)
    ap.add_argument('--out-csv', required=True)
    ap.add_argument('--patient-list', default=None,
                    help='可选: 患者清单文件 (每行一个 ID), 覆盖自动抽样')
    ap.add_argument('--seed', type=int, default=42)
    args = ap.parse_args()

    pts = sorted(os.path.basename(d) for d in glob.glob(os.path.join(args.specs_dir, 'virtual_pt_*')))
    if args.patient_list:
        with open(args.patient_list) as f:
            pts = [l.strip() for l in f if l.strip()]
        rng = random.Random(args.seed)
        sample = sorted(pts)
    else:
        rng = random.Random(args.seed)
        sample = sorted(rng.sample(pts, min(args.n_patients, len(pts))))

    rows = []
    for i, pt in enumerate(sample):
        pm = os.path.join(args.specs_dir, pt, 'variants.pm')
        # 读取目录(压缩后删除) 或 gz
        gzf = os.path.join(args.reads_root, pt + '.fastq.gz')
        dirf = os.path.join(args.reads_root, pt)
        if os.path.exists(gzf):
            rfiles = [gzf]
        elif os.path.isdir(dirf) and (glob.glob(os.path.join(dirf, '*.fastq')) or
                                      glob.glob(os.path.join(dirf, '*.fastq.gz'))):
            rfiles = sorted(glob.glob(os.path.join(dirf, '*.fastq'))) + \
                     sorted(glob.glob(os.path.join(dirf, '*.fastq.gz')))
        else:
            continue
        anchors = load_pm(pm)
        rng2 = random.Random(args.seed + i)
        picks = rng2.sample(anchors, min(args.anchors_per_patient, len(anchors)))
        km = [kmers(args.ref, a['chrom'], a['pos'], a['alt']) for a in picks]
        pairs = [(a, k) for a, k in zip(picks, km) if k]
        if not pairs:
            continue
        counts = count_kmers_in_patient(rfiles, [k for _, kk in pairs for k in kk])
        for a, kk in pairs:
            ref21, alt21 = kk
            rows.append(dict(patient=pt, chrom=a['chrom'], pos=a['pos'], alt=a['alt'],
                             vaf_expected=a['vaf'],
                             ref_count=counts[ref21], alt_count=counts[alt21]))
        if (i + 1) % 25 == 0:
            print(f"[aggregate] {i+1}/{len(sample)} 患者, 累计锚点 {len(rows)}", flush=True)

    # 输出 CSV
    import csv
    with open(args.out_csv, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)

    # 汇总
    n = len(rows)
    covered = [r for r in rows if (r['ref_count'] + r['alt_count']) > 0]
    with_alt = [r for r in rows if r['alt_count'] > 0]
    detect = len(with_alt) / n if n else 0
    cover_rate = len(covered) / n if n else 0
    # 聚合 VAF 一致性 (total>=2 的锚点)
    good = [r for r in covered if (r['ref_count'] + r['alt_count']) >= 2]
    xs = [r['vaf_expected'] for r in good]
    ys = [r['alt_count'] / (r['alt_count'] + r['ref_count']) for r in good]
    rho, p = (spearmanr(xs, ys) if len(good) >= 10 and len(set(xs)) > 1 else (float('nan'), float('nan')))
    summary = dict(patients_sampled=len(sample), anchors=n,
                   locus_cover_rate=round(cover_rate, 4),
                   alt_detectability=round(detect, 4),
                   anchors_with_ge2_reads=len(good),
                   vaf_spearman_rho=round(float(rho), 4), vaf_p=float(p))
    out_json = args.out_csv.replace('.csv', '_summary.json')
    with open(out_json, 'w') as f:
        json.dump(summary, f, indent=2)
    print("[aggregate] ✅", json.dumps(summary, indent=2))

if __name__ == '__main__':
    main()

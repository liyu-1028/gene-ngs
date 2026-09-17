#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
run_hicov_minisim_validation.py — 位点靶向高覆盖 mini 仿真 VAF 验证
对 hicov 队列每名患者: 合成 mini 模型 → synggen 高覆盖仿真 → 锚点 VAF 测量 → 与注入值对比
预期: 观测 VAF ≈ tc × spec VAF (synggen 肿瘤稀释); 覆盖率应 ~100%
输出: data/synthetic_2017/evaluation/hicov_minisim/{per_anchor.tsv, summary.json, minisim_scatter.png}
"""
import os
import gzip, glob, json, os, shutil, subprocess, sys

B = os.environ.get("PROJECT_ROOT", os.path.expanduser("~/gene-ngs"))
FASTA = f"{B}/data/raw/reference/GRCh38_full_analysis_set_plus_decoy_hla.fa"
QM = f"{B}/data/raw/reference/models/wes_model/target_regions.qm.gz"
HICOV = f"{B}/data/synthetic_2017/reads_hicov"
SPECS = f"{B}/data/synthetic_2017/variant_specs"
OUT = f"{B}/data/synthetic_2017/evaluation/hicov_minisim"
DEPTH = 150
TC = 0.3
K = 21
HALF = K // 2
COMP = str.maketrans("ACGT", "TGCA")
rc = lambda s: s.translate(COMP)[::-1]

os.makedirs(OUT, exist_ok=True)
WORK = "/tmp/hicov_minisim"
shutil.rmtree(WORK, ignore_errors=True)
os.makedirs(WORK)

def faidx(region):
    o = subprocess.run(["samtools", "faidx", FASTA, region],
                       capture_output=True, text=True).stdout
    return "".join(o.split("\n")[1:]).upper()

rows = []
patients = sorted(os.path.basename(p) for p in glob.glob(f"{HICOV}/virtual_pt_*"))
print(f"患者数: {len(patients)}", file=sys.stderr)

for pt in patients:
    ext_bed = f"{HICOV}/{pt}/target_regions.extended.bed"
    vdir = f"{SPECS}/{pt}"
    if not os.path.exists(ext_bed) or not os.path.exists(f"{vdir}/variants.pm"):
        print(f"  {pt}: 缺文件, 跳过", file=sys.stderr)
        continue

    # ── 构建模型 ──
    r = subprocess.run(["python3", f"{B}/src/stage0/build_mini_model.py",
                        ext_bed, f"{vdir}/variants.pm", f"{WORK}/{pt}", str(DEPTH)],
                       capture_output=True, text=True)
    if r.returncode != 0:
        print(f"  {pt}: builder失败 {r.stderr[:200]}", file=sys.stderr)
        continue
    nreads = int([l for l in r.stdout.split("\n") if "NREADS" in l][0].split("=")[1])

    # ── 仿真 ──
    outd = f"{WORK}/{pt}/out"
    cmd = [f"{B}/bin/synggen", "mode=1", f"bed={WORK}/{pt}/mini.bed", f"fasta={FASTA}",
           f"rdm={WORK}/{pt}/mini.rdm.gz", f"pbe={WORK}/{pt}/mini.pbe.gz", f"qm={QM}",
           f"pm={vdir}/variants.pm", f"indel={vdir}/variants.indel", f"cna={vdir}/variants.cna",
           "threads=2", f"nreads={nreads}", "tc=0.3", f"out={outd}", "seed=42"]
    r = subprocess.run(cmd, capture_output=True, text=True)
    fqs = sorted(glob.glob(f"{outd}/*.fastq"))
    if not fqs:
        print(f"  {pt}: synggen失败 {r.stderr[:200]}", file=sys.stderr)
        continue

    # ── 区域 (最终 ±100) ──
    regions = []
    with open(ext_bed) as f:
        for line in f:
            c, s, e = line.split()[:3]
            regions.append((c, int(s) - 100, int(e) + 100))

    # ── 锚点: ALT≠REF 且落在区域内 ──
    anchors = []
    with open(f"{vdir}/variants.pm") as f:
        for line in f:
            p = line.split()
            if len(p) < 4:
                continue
            c, pos, alt, vaf = p[0], int(p[1]), p[2], float(p[3])
            if not any(mc == c and ms <= pos < me for mc, ms, me in regions):
                continue
            ref = faidx(f"{c}:{pos - HALF}-{pos + HALF}")
            refb = faidx(f"{c}:{pos}-{pos}")
            if refb == alt:
                continue  # 无效突变
            anchors.append((c, pos, refb, alt, vaf, ref,
                            ref[:HALF] + alt + ref[HALF + 1:]))

    # ── k-mer 计数 (grep -oF -f, 双端序列) ──
    pats, keys = [], []
    for i, (c, pos, refb, alt, vaf, refk, altk) in enumerate(anchors):
        for tag, km in (("R", refk), ("A", altk), ("Rrc", rc(refk)), ("Arc", rc(altk))):
            pats.append(km)
            keys.append((i, tag))
    pf = f"{WORK}/{pt}/pats.txt"
    with open(pf, "w") as f:
        f.write("\n".join(pats) + "\n")
    counts = {}
    for j, km in enumerate(pats):
        counts[j] = 0
    for fq in fqs:
        r = subprocess.run(["grep", "-oF", "-f", pf, fq], capture_output=True, text=True)
        if r.stdout:
            from collections import Counter
            hit = Counter(r.stdout.split("\n")[:-1])
            for kmer, n in hit.items():
                counts[pats.index(kmer)] += n
    ref_n = {i: 0 for i in range(len(anchors))}
    alt_n = {i: 0 for i in range(len(anchors))}
    for j, (i, tag) in enumerate(keys):
        if tag in ("R", "Rrc"):
            ref_n[i] += counts[j]
        else:
            alt_n[i] += counts[j]

    for i, (c, pos, refb, alt, vaf, refk, altk) in enumerate(anchors):
        rn, an = ref_n[i], alt_n[i]
        tot = rn + an
        obs = an / tot if tot else None
        rows.append(dict(patient=pt, chrom=c, pos=pos, ref=refb, alt=alt,
                         spec_vaf=vaf, tc_vaf=round(TC * vaf, 4),
                         ref_kmer=refk, alt_kmer=altk,
                         ref_count=rn, alt_count=an, depth=tot,
                         obs_vaf=round(obs, 4) if obs is not None else None,
                         detected=(tot > 0 and an > 0)))
    print(f"  {pt}: {len(anchors)}锚点 nreads={nreads} (仿真+计数完成)", file=sys.stderr)

# ── 汇总 ──
with open(f"{OUT}/per_anchor.tsv", "w") as f:
    f.write("patient\tchrom\tpos\tref\talt\tspec_vaf\ttc_vaf\tref_count\talt_count\tdepth\tobs_vaf\tdetected\n")
    for r0 in rows:
        f.write(f"{r0['patient']}\t{r0['chrom']}\t{r0['pos']}\t{r0['ref']}\t{r0['alt']}\t"
                f"{r0['spec_vaf']}\t{r0['tc_vaf']}\t{r0['ref_count']}\t{r0['alt_count']}\t"
                f"{r0['depth']}\t{r0['obs_vaf']}\t{r0['detected']}\n")

import numpy as np
tcv = np.array([r0["tc_vaf"] for r0 in rows])
obv = np.array([r0["obs_vaf"] if r0["obs_vaf"] is not None else np.nan for r0 in rows])
dep = np.array([r0["depth"] for r0 in rows])
ok = ~np.isnan(obv)
from scipy.stats import pearsonr, spearmanr
pr, pp = pearsonr(tcv[ok], obv[ok])
sr, sp = spearmanr(tcv[ok], obv[ok])
det_rate = np.mean([r0["detected"] for r0 in rows])
cov_rate = np.mean(dep > 0)
summary = dict(
    n_patients=len(patients), n_anchors=len(rows),
    anchors_with_depth=int(np.sum(dep > 0)), coverage_rate=round(float(cov_rate), 4),
    alt_detected_rate=round(float(det_rate), 4),
    median_depth=float(np.median(dep)),
    pearson_r=round(float(pr), 4), pearson_p=float(pp),
    spearman_rho=round(float(sr), 4), spearman_p=float(sp),
    mean_abs_dev_vs_tc=round(float(np.nanmean(np.abs(obv[ok] - tcv[ok]))), 4),
    tc=TC, target_depth=DEPTH,
)
with open(f"{OUT}/summary.json", "w") as f:
    json.dump(summary, f, indent=2)

# ── 散点图 ──
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
fig, ax = plt.subplots(figsize=(6, 6))
ax.scatter(tcv[ok], obv[ok], alpha=0.7, s=42, edgecolor="k", linewidth=0.5)
lim = max(0.35, float(np.nanmax(obv[ok])) * 1.05)
ax.plot([0, lim], [0, lim], "k--", lw=1, label="y = x")
xs = np.linspace(0, lim, 50)
okm = ~np.isnan(obv)
z = np.polyfit(tcv[okm], obv[okm], 1)
ax.plot(xs, np.polyval(z, xs), "r-", lw=1.5,
        label=f"fit: y={z[0]:.2f}x{z[1]:+.3f}")
ax.set_xlabel(f"Injected VAF (tc-adjusted = {TC}×spec)")
ax.set_ylabel("Observed VAF (k-mer, mini-sim ~150×)")
ax.set_title(f"HiCov mini-sim VAF recovery\n{len(rows)} anchors, {len(patients)} patients; "
             f"Pearson r={pr:.3f} (p={pp:.1e})")
ax.legend()
ax.set_xlim(0, lim)
ax.set_ylim(0, lim)
fig.tight_layout()
fig.savefig(f"{OUT}/minisim_scatter.png", dpi=150)

print(json.dumps(summary, indent=2))
print(f"\n输出: {OUT}/per_anchor.tsv, summary.json, minisim_scatter.png")

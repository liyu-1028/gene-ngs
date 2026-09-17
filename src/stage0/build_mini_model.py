#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
build_mini_model.py — 为位点靶向高覆盖仿真合成 synggen 兼容模型
输入: 患者 extended.bed + variants.pm + rlen(100)
输出: mini.bed(原始区域, synggen 自行 ±rlen 扩展), mini.rdm.gz, mini.pbe.gz
原理: 源码 model_rdm.c 校验 (1)RDM行数==最终区域数 (2)每行列数==区域长度;
      generate_reads.c 按累计深度二分分配 reads → 合成均匀累计深度即可。
用法: python build_mini_model.py <extended.bed> <variants.pm> <outdir> [target_depth]
"""
import gzip, os, sys

def main():
    ext_bed, pm, outdir = sys.argv[1], sys.argv[2], sys.argv[3]
    depth = int(sys.argv[4]) if len(sys.argv) > 4 else 150
    rlen = 100
    os.makedirs(outdir, exist_ok=True)

    regions = []
    with open(ext_bed) as f:
        for line in f:
            c, s, e = line.split()[:3]
            regions.append((c, int(s), int(e)))

    muts = []
    with open(pm) as f:
        for line in f:
            p = line.split()
            if len(p) >= 3:
                muts.append((p[0], int(p[1])))

    # 1) 含突变的 extended 区域 (保持文件顺序)
    keep = []
    for c, s, e in regions:
        if any(mc == c and s <= mp < e for mc, mp in muts):
            keep.append((c, s, e))
    if not keep:
        sys.exit("no regions with mutations")

    # 2) 最终区域 = ±rlen (检查无合并)
    final = [(c, s - rlen, e + rlen) for c, s, e in keep]
    for a, b in zip(final, final[1:]):
        assert a[0] != b[0] or a[2] < b[1], "regions would merge: %s %s" % (a, b)

    # 3) mini.bed (输入, 未扩展)
    with open(os.path.join(outdir, "mini.bed"), "w") as f:
        for c, s, e in keep:
            f.write(f"{c}\t{s}\t{e}\n")

    # 4) mini.rdm: 逐碱基原始深度 (synggen 内部 computeRDMCumulative 会自行累计)
    #    格式: aN(=区域总深度);d;d;...;d   → 均匀深度, 区域份额 ∝ L×d
    #    注意: 累计和必须 < 2^31 (int32), d=20 时总深 ~13万, 安全
    L_total = sum(e - s for _, s, e in final)
    nreads = depth * L_total // rlen
    D = 20  # 每碱基深度权重
    with gzip.open(os.path.join(outdir, "mini.rdm.gz"), "wt") as f:
        for c, s, e in final:
            L = e - s
            f.write(";".join([str(L * D)] + [str(D)] * L) + "\n")

    # 5) mini.pbe: NA (无系统性配对碱基错误)
    with gzip.open(os.path.join(outdir, "mini.pbe.gz"), "wt") as f:
        for _ in final:
            f.write("NA\n")

    print(f"regions={len(final)} L_total={L_total}bp 建议nreads={nreads} (~{depth}x)")
    print(f"OUTPUT: {outdir}/mini.bed,{outdir}/mini.rdm.gz,{outdir}/mini.pbe.gz,NREADS={nreads}")

if __name__ == "__main__":
    main()

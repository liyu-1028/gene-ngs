#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
convert_msk_impact_to_pipeline.py  (Stage 0: 真实数据源接入)

将 cBioPortal DataHub 的 MSK-IMPACT 2017 (Zehir et al., Nat Med 2017)
原始发布文件转换为项目流水线标准模式
（与 data/processed/clinical_features_dataset.csv 完全同列）：

    somatic_variants                  -> [{'gene', 'alteration', 'vaf'}, ...]
    source_file                       -> SAMPLE_ID
    patient_info.age                  -> 公开版未提供（NaN，见 LIMITATIONS）
    patient_info.gender               -> SEX (Female->女 / Male->男)
    patient_info.pathological_diagnosis -> CANCER_TYPE_DETAILED（英文原文）
    patient_info.sample_type          -> "{SAMPLE_TYPE} / {SPECIMEN_PRESERVATION_TYPE}"
    biomarkers.tmb                    -> TMB_NONSYNONYMOUS (mut/Mb，官方值)
    biomarkers.msi                    -> 公开版未提供（空，见 LIMITATIONS）
    qc_metrics.overall_quality        -> PASS（所有已发布样本均通过 MSK 临床 QC）
    qc_metrics                        -> 空

somatic_variants 组装规则（临床报告级别）：
  1. SNV/InDel  : data_mutations.txt，仅保留可报告后果类别
                  (Missense/Nonsense/Frame_Shift_*/In_Frame_*/Splice_Site/
                   Translation_Start_Site/Nonstop)；VAF = t_alt/(t_ref+t_alt)；
                  alteration 优先 HGVSp_Short，其次 HGVSc，退化为 基因组坐标表示
  2. 融合(SV)   : data_sv.txt，alteration = Event_Info（如 "ALK-DDX11L2 Fusion"），
                  VAF = Tumor_Variant_Count/Tumor_Read_Count（可算则算）
  3. CNA        : data_cna.txt (GISTIC 离散值)，2 -> Amplification，-2 -> Deep Deletion，
                  vaf=None（与原 24 例真实数据中 Amplification 条目一致）

用法:
    python scripts/convert_msk_impact_to_pipeline.py \
        [--src data/msk_impact_2017] [--out data/real/msk_impact_2017]
"""

import argparse
import os
import sys
import json
import pandas as pd

# ---- 可报告的突变后果类别（临床报告级别，排除 Silent/UTR/Intron 等） ----
REPORTABLE_VC = {
    'Missense_Mutation',
    'Nonsense_Mutation',
    'Frame_Shift_Del',
    'Frame_Shift_Ins',
    'In_Frame_Del',
    'In_Frame_Ins',
    'Splice_Site',
    'Translation_Start_Site',
    'Nonstop_Mutation',
}


def load_tables(src: str):
    """读取 cBioPortal 标准发布文件（跳过 # 元信息头）。"""
    pat = pd.read_csv(os.path.join(src, 'data_clinical_patient.txt'),
                      sep='\t', skiprows=4)
    sam = pd.read_csv(os.path.join(src, 'data_clinical_sample.txt'),
                      sep='\t', skiprows=4)
    mut = pd.read_csv(os.path.join(src, 'data_mutations.txt'),
                      sep='\t', skiprows=2, low_memory=False)
    sv = pd.read_csv(os.path.join(src, 'data_sv.txt'), sep='\t')
    cna = pd.read_csv(os.path.join(src, 'data_cna.txt'), sep='\t', index_col=0)
    return pat, sam, mut, sv, cna


def build_snv(mut: pd.DataFrame) -> dict:
    """sample -> list of variant dicts（仅可报告类别）"""
    out = {}
    kept = mut[mut['Variant_Classification'].isin(REPORTABLE_VC)].copy()

    # 计数必须非负且分母>0（源数据存在 1 条 t_ref_count=-2 的瑕疵记录）
    ok = ((pd.to_numeric(kept['t_ref_count'], errors='coerce') >= 0) &
          (pd.to_numeric(kept['t_alt_count'], errors='coerce') >= 0))
    denom = (kept['t_ref_count'] + kept['t_alt_count']).where(ok)
    kept['VAF'] = (kept['t_alt_count'] / denom).astype(float)

    def alter(r):
        if isinstance(r['HGVSp_Short'], str) and r['HGVSp_Short']:
            return r['HGVSp_Short']
        if isinstance(r['HGVSc'], str) and r['HGVSc']:
            return r['HGVSc']
        return (f"g.{r['Chromosome']}:{int(r['Start_Position'])}"
                f"{r['Reference_Allele']}>{r['Tumor_Seq_Allele2']}")

    for sid, g in kept.groupby('Tumor_Sample_Barcode', sort=False):
        vs = [{'gene': r['Hugo_Symbol'],
               'alteration': alter(r),
               'vaf': None if pd.isna(r['VAF']) else round(float(r['VAF']), 4)}
              for _, r in g.iterrows()]
        out[sid] = vs
    return out


def build_fusion(sv: pd.DataFrame) -> dict:
    """sample -> list of fusion variant dicts"""
    out = {}
    for sid, g in sv.groupby('Sample_Id', sort=False):
        vs = []
        for _, r in g.iterrows():
            rc = pd.to_numeric(r['Tumor_Read_Count'], errors='coerce')
            vc = pd.to_numeric(r['Tumor_Variant_Count'], errors='coerce')
            # 仅当支持读数落在 (0, read_count) 开区间内才给出 VAF；
            # vc==0（无支持读数）或 vc>=rc（源数据计数口径不一致）置 None
            if pd.notna(rc) and pd.notna(vc) and 0 < vc < rc:
                vaf = round(float(vc / rc), 4)
                if vaf == 0:          # 极小 VAF (如 1/505392) 四舍五入为 0，无报告意义
                    vaf = None
            else:
                vaf = None
            ev = r['Event_Info']
            if not isinstance(ev, str) or not ev:
                # Event_Info 缺失时退化为 "Site1-Site2 Fusion"
                ev = f"{r['Site1_Hugo_Symbol']}-{r['Site2_Hugo_Symbol']} Fusion"
            vs.append({'gene': r['Site1_Hugo_Symbol'],
                       'alteration': ev,
                       'vaf': vaf})
        out[sid] = vs
    return out


def build_cna(cna: pd.DataFrame) -> dict:
    """sample -> list of CNA variant dicts（仅 2=Amplification / -2=Deep Deletion）"""
    out = {}
    for sid in cna.columns:
        col = cna[sid]
        amps = col.index[col == 2].tolist()
        dels = col.index[col == -2].tolist()
        vs = [{'gene': g, 'alteration': 'Amplification', 'vaf': None} for g in amps]
        vs += [{'gene': g, 'alteration': 'Deep Deletion', 'vaf': None} for g in dels]
        if vs:
            out[sid] = vs
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--src', default='data/msk_impact_2017')
    ap.add_argument('--out', default='data/real/msk_impact_2017')
    args = ap.parse_args()

    src, out = args.src, args.out
    os.makedirs(out, exist_ok=True)

    print(f"[1/5] 读取 cBioPortal 发布文件: {src}")
    pat, sam, mut, sv, cna = load_tables(src)
    print(f"      患者 {pat.PATIENT_ID.nunique():,} | 样本 {sam.SAMPLE_ID.nunique():,} | "
          f"突变行 {len(mut):,} | 融合 {len(sv):,} | CNA 基因 {cna.shape[0]}")

    print("[2/5] 组装三类变异 ...")
    snv = build_snv(mut)
    fus = build_fusion(sv)
    cnv = build_cna(cna)
    print(f"      SNV/InDel 样本 {len(snv):,} | 融合样本 {len(fus):,} | CNA 样本 {len(cnv):,}")
    n_snv_tot = sum(len(v) for v in snv.values())
    n_fus_tot = sum(len(v) for v in fus.values())
    n_cnv_tot = sum(len(v) for v in cnv.values())
    n_var_tot = n_snv_tot + n_fus_tot + n_cnv_tot
    vafs = [x['vaf'] for vs in snv.values() for x in vs if x['vaf'] is not None]
    vaf_series = pd.Series(
        [len(snv.get(s, [])) + len(fus.get(s, [])) + len(cnv.get(s, []))
         for s in sam['SAMPLE_ID']])

    print("[3/5] 合并临床表 (patient x sample) ...")
    df = sam.merge(pat, on='PATIENT_ID', how='left')

    def gender(x):
        return {'Female': '女', 'Male': '男'}.get(x, x)

    def sample_type(r):
        return f"{r['SAMPLE_TYPE']} / {r['SPECIMEN_PRESERVATION_TYPE']}"

    print("[4/5] 生成流水线模式数据集 ...")
    df['somatic_variants'] = df['SAMPLE_ID'].map(
        lambda s: repr(snv.get(s, []) + fus.get(s, []) + cnv.get(s, [])))
    df['source_file'] = df['SAMPLE_ID']
    df['patient_info.age'] = pd.NA                       # 公开版未提供
    df['patient_info.gender'] = df['SEX'].map(gender)
    df['patient_info.pathological_diagnosis'] = df['CANCER_TYPE_DETAILED']
    df['patient_info.sample_type'] = df.apply(sample_type, axis=1)
    df['biomarkers.tmb'] = df['TMB_NONSYNONYMOUS']
    df['biomarkers.msi'] = ''                            # 公开版未提供
    df['qc_metrics.overall_quality'] = 'PASS'            # 均已通过 MSK 临床 QC
    df['qc_metrics'] = ''

    cols = ['somatic_variants', 'source_file', 'patient_info.age',
            'patient_info.gender', 'patient_info.pathological_diagnosis',
            'patient_info.sample_type', 'biomarkers.tmb', 'biomarkers.msi',
            'qc_metrics.overall_quality', 'qc_metrics']
    real = df[cols].copy()

    out_csv = os.path.join(out, 'real_dataset.csv')
    real.to_csv(out_csv, index=False, encoding='utf-8')
    print(f"      -> {out_csv}  ({len(real):,} 行 × {len(real.columns)} 列)")

    # ---- 另存一份带完整临床信息的研究级宽表（供生存分析/TSTR 等扩展使用） ----
    clin_cols = ['PATIENT_ID', 'SAMPLE_ID', 'SEX', 'VITAL_STATUS', 'SMOKING_HISTORY',
                 'OS_MONTHS', 'OS_STATUS', 'SAMPLE_TYPE', 'SPECIMEN_PRESERVATION_TYPE',
                 'DNA_INPUT', 'SAMPLE_COVERAGE', 'TUMOR_PURITY', 'MATCHED_STATUS',
                 'PRIMARY_SITE', 'METASTATIC_SITE', 'ONCOTREE_CODE',
                 'CANCER_TYPE', 'CANCER_TYPE_DETAILED', 'TMB_NONSYNONYMOUS']
    clin = df[[c for c in clin_cols if c in df.columns]].copy()
    clin['n_snv'] = clin['SAMPLE_ID'].map(lambda s: len(snv.get(s, [])))
    clin['n_fusion'] = clin['SAMPLE_ID'].map(lambda s: len(fus.get(s, [])))
    clin['n_cna'] = clin['SAMPLE_ID'].map(lambda s: len(cnv.get(s, [])))
    out_clin = os.path.join(out, 'real_dataset_clinical_wide.csv')
    clin.to_csv(out_clin, index=False, encoding='utf-8')
    print(f"      -> {out_clin}  ({len(clin):,} 行，含全部临床列)")

    # ---- 质量汇总 ----
    # 注: SNV 存在 3 条 t_ref_count==0 的源记录，VAF=1.0 属实 (FFPE 伪影/低纯度)，保留
    report = {
        'source': 'cBioPortal DataHub msk_impact_2017 (ODbL; Zehir et al., Nat Med 2017, PMID 28481359)',
        'rows_samples': int(len(real)),
        'rows_patients': int(pat.PATIENT_ID.nunique()),
        'variants_total': int(n_var_tot),
        'variants_snv_indel': int(n_snv_tot),
        'variants_fusion': int(n_fus_tot),
        'variants_cna': int(n_cnv_tot),
        'variants_per_sample_median': float(vaf_series.median()),
        'vaf_computable_snv_pct': round(100 * len(vafs) / max(1, n_snv_tot), 2),
        'vaf_median': round(float(pd.Series(vafs).median()), 4),
        'limitations': [
            'patient_info.age: 2017 公开版未包含年龄字段，置 NaN',
            'biomarkers.msi: 公开版无 MSI 状态，置空（stage3 审计时 TMB-MSI 规则需相应放宽）',
            'somatic_variants 仅含临床可报告类别，Silent/UTR/Intron 未纳入（TMB 列仍用官方 TMB_NONSYNONYMOUS 全口径）',
            '样本为全瘤组织 (无血液对照)，sample_type 采用英文 "{SAMPLE_TYPE} / {SPECIMEN_PRESERVATION_TYPE}" 口径',
            'SNV 中极少数源记录 t_ref_count<=0 → VAF 为 1.0 或置空 (FFPE 伪影/源数据瑕疵)，如实处理；融合 VAF 仅在 0<vc<rc 时给出',
        ],
    }
    with open(os.path.join(out, 'CONVERSION_REPORT.json'), 'w', encoding='utf-8') as f:
        json.dump(report, f, ensure_ascii=False, indent=2)
    print("[5/5] 转换汇总:")
    print(json.dumps({k: v for k, v in report.items() if k != 'limitations'},
                     ensure_ascii=False, indent=2))
    for lim in report['limitations']:
        print(f"  ⚠ {lim}")
    return 0


if __name__ == '__main__':
    sys.exit(main())

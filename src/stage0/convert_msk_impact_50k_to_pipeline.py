#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
convert_msk_impact_50k_to_pipeline.py  (Stage 0: 真实数据源接入 v2)

将 cBioPortal DataHub 的 MSK-IMPACT 50K (Bandlamudi et al., Cancer Cell 2026,
studyId=msk_impact_50k_2026, 54,331 肿瘤样本 + 配对正常) 原始发布文件转换为
项目流水线标准模式（与 data/processed/clinical_features_dataset.csv 同列）：

    somatic_variants                  -> [{'gene', 'alteration', 'vaf'}, ...]
    source_file                       -> SAMPLE_ID
    patient_info.age                  -> AGE_AT_DX ✓ (2017 版缺失，此版已提供)
    patient_info.gender               -> SEX (Female->女 / Male->男)
    patient_info.pathological_diagnosis -> CANCER_TYPE_DETAILED
    patient_info.sample_type          -> SAMPLE_TYPE (Primary/Metastasis/...)
    biomarkers.tmb                    -> TMB_SCORE (连续值, mut/Mb, 官方口径)
    biomarkers.msi                    -> MSI_TYPE (Stable->MSS, Instable->MSI-H,
                                         Indeterminate/Do not report->'') ✓
    qc_metrics.overall_quality        -> FACETS_QC (Pass->PASS, Fail->FAIL) ✓
    qc_metrics                        -> 空

somatic_variants 组装规则（临床报告级别，与 2017 版脚本一致）：
  1. SNV/InDel  : data_mutations.txt（此版无 '#' 注释行, skiprows=0），
                  仅保留可报告后果类别；VAF = t_alt/(t_ref+t_alt)，
                  计数非负且分母>0 才计算
  2. 融合(SV)   : 此版无 Tumor_Read_Count，VAF = Tumor_Variant_Count /
                  (Split_Read_Support + Paired_End_Read_Support)，
                  仅当 0 < vaf < 1 时给出
  3. CNA        : data_cna.txt (GISTIC 离散值)，2 -> Amplification，-2 -> Deep Deletion

用法:
    python src/stage0/convert_msk_impact_50k_to_pipeline.py \
        [--src data/msk_impact_50k_2026] [--out data/real/msk_impact_50k_2026]
"""

import argparse
import os
import sys
import json
import pandas as pd

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

MSI_MAP = {'Stable': 'MSS', 'Instable': 'MSI-H'}   # 其余 (Indeterminate/Do not report/'') -> ''
QC_MAP = {'Pass': 'PASS', 'Fail': 'FAIL'}          # 其余 ('') -> 'UNKNOWN'


def load_tables(src: str):
    pat = pd.read_csv(os.path.join(src, 'data_clinical_patient.txt'),
                      sep='\t', skiprows=4, dtype=str)
    sam = pd.read_csv(os.path.join(src, 'data_clinical_sample.txt'),
                      sep='\t', skiprows=4, dtype=str)
    # 50K 版 MAF 无 '#' 注释行
    mut = pd.read_csv(os.path.join(src, 'data_mutations.txt'),
                      sep='\t', skiprows=0, low_memory=False)
    sv = pd.read_csv(os.path.join(src, 'data_sv.txt'), sep='\t')
    cna = pd.read_csv(os.path.join(src, 'data_cna.txt'), sep='\t', index_col=0)
    return pat, sam, mut, sv, cna


def build_snv(mut: pd.DataFrame) -> dict:
    out = {}
    kept = mut[mut['Variant_Classification'].isin(REPORTABLE_VC)].copy()

    ok = ((pd.to_numeric(kept['t_ref_count'], errors='coerce') >= 0) &
          (pd.to_numeric(kept['t_alt_count'], errors='coerce') >= 0))
    denom = (pd.to_numeric(kept['t_ref_count'], errors='coerce') +
             pd.to_numeric(kept['t_alt_count'], errors='coerce')).where(ok & (kept['t_ref_count'].notna()))
    denom = denom.where(denom > 0)
    kept['VAF'] = (pd.to_numeric(kept['t_alt_count'], errors='coerce') / denom)

    def alter(r):
        if isinstance(r['HGVSp_Short'], str) and r['HGVSp_Short']:
            return r['HGVSp_Short']
        if isinstance(r['HGVSc'], str) and r['HGVSc']:
            return r['HGVSc']
        return (f"g.{r['Chromosome']}:{int(float(r['Start_Position']))}"
                f"{r['Reference_Allele']}>{r['Tumor_Seq_Allele2']}")

    for sid, g in kept.groupby('Tumor_Sample_Barcode', sort=False):
        vs = [{'gene': r['Hugo_Symbol'],
               'alteration': alter(r),
               'vaf': None if pd.isna(r['VAF']) else round(float(r['VAF']), 4)}
              for _, r in g.iterrows()]
        out[sid] = vs
    return out


def build_fusion(sv: pd.DataFrame) -> dict:
    out = {}
    split = pd.to_numeric(sv['Split_Read_Support'], errors='coerce')
    paired = pd.to_numeric(sv['Paired_End_Read_Support'], errors='coerce')
    vc = pd.to_numeric(sv['Tumor_Variant_Count'], errors='coerce')
    denom = (split.fillna(0) + paired.fillna(0))
    vaf = (vc / denom.where(denom > 0))
    vaf = vaf.where((vaf > 0) & (vaf < 1))          # 无效/口径异常 -> NaN -> None

    for i, r in sv.iterrows():
        if pd.isna(r['Site1_Hugo_Symbol']):   # 50K 源存在 9 条无基因符号的 SV，跳过
            continue
        ev = r['Event_Info']
        if not isinstance(ev, str) or not ev:
            ev = f"{r['Site1_Hugo_Symbol']}-{r['Site2_Hugo_Symbol']} Fusion"
        v = None if pd.isna(vaf.iloc[i]) else round(float(vaf.iloc[i]), 4)
        out.setdefault(r['Sample_Id'], []).append(
            {'gene': r['Site1_Hugo_Symbol'], 'alteration': ev, 'vaf': v})
    return out


def build_cna(cna: pd.DataFrame) -> dict:
    out = {}
    for sid in cna.columns:
        col = cna[sid]
        vs = ([{'gene': g, 'alteration': 'Amplification', 'vaf': None}
               for g in col.index[col == 2].tolist()] +
              [{'gene': g, 'alteration': 'Deep Deletion', 'vaf': None}
               for g in col.index[col == -2].tolist()])
        if vs:
            out[sid] = vs
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('--src', default='data/msk_impact_50k_2026')
    ap.add_argument('--out', default='data/real/msk_impact_50k_2026')
    args = ap.parse_args()
    src, out = args.src, args.out
    os.makedirs(out, exist_ok=True)

    print(f"[1/6] 读取 cBioPortal 发布文件: {src}")
    pat, sam, mut, sv, cna = load_tables(src)
    print(f"      患者 {pat.PATIENT_ID.nunique():,} | 样本 {sam.SAMPLE_ID.nunique():,} | "
          f"突变行 {len(mut):,} | 融合 {len(sv):,} | CNA 基因 {cna.shape[0]}")

    print("[2/6] 组装三类变异 ...")
    snv = build_snv(mut)
    fus = build_fusion(sv)
    cnv = build_cna(cna)
    print(f"      SNV/InDel 样本 {len(snv):,} | 融合样本 {len(fus):,} | CNA 样本 {len(cnv):,}")

    print("[3/6] 合并临床表 ...")
    df = sam.merge(pat, on='PATIENT_ID', how='left', suffixes=('', '_pat'))

    print("[4/6] 生成流水线模式数据集 ...")
    df['somatic_variants'] = df['SAMPLE_ID'].map(
        lambda s: repr(snv.get(s, []) + fus.get(s, []) + cnv.get(s, [])))
    df['source_file'] = df['SAMPLE_ID']
    df['patient_info.age'] = pd.to_numeric(df['AGE_AT_DX'], errors='coerce')
    df['patient_info.gender'] = df['SEX'].map({'Female': '女', 'Male': '男'})
    df['patient_info.pathological_diagnosis'] = df['CANCER_TYPE_DETAILED']
    df['patient_info.sample_type'] = df['SAMPLE_TYPE']
    df['biomarkers.tmb'] = pd.to_numeric(df['TMB_SCORE'], errors='coerce').round(4)
    df['biomarkers.msi'] = df['MSI_TYPE'].map(MSI_MAP).fillna('')
    df['qc_metrics.overall_quality'] = df['FACETS_QC'].map(QC_MAP).fillna('UNKNOWN')
    df['qc_metrics'] = ''

    cols = ['somatic_variants', 'source_file', 'patient_info.age',
            'patient_info.gender', 'patient_info.pathological_diagnosis',
            'patient_info.sample_type', 'biomarkers.tmb', 'biomarkers.msi',
            'qc_metrics.overall_quality', 'qc_metrics']
    real = df[cols].copy()

    out_csv = os.path.join(out, 'real_dataset.csv')
    real.to_csv(out_csv, index=False, encoding='utf-8')
    print(f"      -> {out_csv}  ({len(real):,} 行 × {len(real.columns)} 列)")

    # ---- 研究级宽表 ----
    clin_cols = ['PATIENT_ID', 'SAMPLE_ID', 'SEX', 'AGE_AT_DX', 'OS_STATUS', 'OS_MONTHS',
                 'ANCESTRY_LABEL', 'SAMPLE_TYPE', 'DISEASE_STATUS', 'PRIMARY_SITE',
                 'METASTATIC_SITE', 'ONCOTREE_CODE', 'GENE_PANEL', 'SAMPLE_COVERAGE',
                 'TUMOR_PURITY', 'FACETS_PURITY', 'FACETS_PLOIDY', 'FACETS_WGD',
                 'MSI_SCORE', 'MSI_TYPE', 'TMB_SCORE', 'FACETS_QC',
                 'CANCER_TYPE', 'CANCER_TYPE_DETAILED']
    clin = df[[c for c in clin_cols if c in df.columns]].copy()
    for c in ['AGE_AT_DX', 'OS_MONTHS', 'MSI_SCORE', 'TMB_SCORE',
              'SAMPLE_COVERAGE', 'TUMOR_PURITY', 'FACETS_PURITY', 'FACETS_PLOIDY']:
        if c in clin.columns:
            clin[c] = pd.to_numeric(clin[c], errors='coerce')
    clin['n_snv'] = clin['SAMPLE_ID'].map(lambda s: len(snv.get(s, [])))
    clin['n_fusion'] = clin['SAMPLE_ID'].map(lambda s: len(fus.get(s, [])))
    clin['n_cna'] = clin['SAMPLE_ID'].map(lambda s: len(cnv.get(s, [])))
    out_clin = os.path.join(out, 'real_dataset_clinical_wide.csv')
    clin.to_csv(out_clin, index=False, encoding='utf-8')
    print(f"      -> {out_clin}  ({len(clin):,} 行)")

    # ---- 质量汇总 ----
    n_snv_tot = sum(len(v) for v in snv.values())
    n_fus_tot = sum(len(v) for v in fus.values())
    n_cnv_tot = sum(len(v) for v in cnv.values())
    vafs = [x['vaf'] for vs in snv.values() for x in vs if x['vaf'] is not None]
    vaf_series = pd.Series([len(snv.get(s, [])) + len(fus.get(s, [])) + len(cnv.get(s, []))
                            for s in sam['SAMPLE_ID']])
    age_cov = real['patient_info.age'].notna().mean() * 100
    msi_cov = (real['biomarkers.msi'] != '').mean() * 100
    report = {
        'source': ('cBioPortal DataHub msk_impact_50k_2026 '
                   '(CC BY-NC-ND 4.0; Bandlamudi et al., Cancer Cell 2026, PMID 41895280)'),
        'rows_samples': int(len(real)),
        'rows_patients': int(pat.PATIENT_ID.nunique()),
        'variants_total': int(n_snv_tot + n_fus_tot + n_cnv_tot),
        'variants_snv_indel': int(n_snv_tot),
        'variants_fusion': int(n_fus_tot),
        'variants_cna': int(n_cnv_tot),
        'variants_per_sample_median': float(vaf_series.median()),
        'vaf_computable_snv_pct': round(100 * len(vafs) / max(1, n_snv_tot), 2),
        'vaf_median': round(float(pd.Series(vafs).median()), 4),
        'age_coverage_pct': round(age_cov, 2),
        'msi_coverage_pct': round(msi_cov, 2),
        'limitations': [
            '许可由 2017 版 ODbL 变更为 CC BY-NC-ND 4.0：署名-非商业-禁止演绎，发布衍生/合成数据前需法务确认 ND 条款边界',
            'MSI_TYPE 为三态化映射: Stable->MSS, Instable->MSI-H, Indeterminate/Do not report->空 (连续 MSI_SCORE 保留在宽表)',
            'qc_metrics.overall_quality 采用 FACETS_QC (Pass/Fail)，不再是统一 PASS',
            'somatic_variants 仅含临床可报告类别 (TMB_SCORE 仍为官方全口径)',
            '融合 VAF 分母为 Split+Paired 支持读数 (此版无 Tumor_Read_Count)，仅 0<vaf<1 时给出',
        ],
    }
    with open(os.path.join(out, 'CONVERSION_REPORT.json'), 'w', encoding='utf-8') as f:
        json.dump(report, f, ensure_ascii=False, indent=2)
    print("[5/6] 转换汇总:")
    print(json.dumps({k: v for k, v in report.items() if k != 'limitations'},
                     ensure_ascii=False, indent=2))
    print("[6/6] 限制:")
    for lim in report['limitations']:
        print(f"  ⚠ {lim}")
    return 0


if __name__ == '__main__':
    sys.exit(main())

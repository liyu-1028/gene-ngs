#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Stage 5 v2: 路线参数化评估 (2017 发表路线 / 50K 扩展路线)
子命令: stats | tstr | privacy | all
指标口径与 src/stage5 v1 保持一致, 但:
  - 真实/合成路径与输出目录全部 CLI 参数化
  - 字段自适应: 全空列 (如 2017 的 age/msi) 自动剔除, 不进 KS/相关性/模型特征
输出: <eval_dir>/statistical_fidelity_metrics.json, manifold_overlap.png,
      utility_tstr_metrics.json, privacy_metrics.json
"""
import argparse, json, os, sys, ast
import numpy as np
import pandas as pd
from scipy.stats import ks_2samp, spearmanr

PROJECT_BASE = os.path.abspath(os.path.join(os.path.dirname(__file__), "../.."))
sys.path.insert(0, PROJECT_BASE)

from src.stage5.statistical_fidelity import load_data, compute_tvd, get_gene_frequencies
from src.stage5.privacy_assessment import compute_dcr_and_nps, membership_inference_attack

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.preprocessing import StandardScaler, OneHotEncoder
from sklearn.compose import ColumnTransformer
from sklearn.pipeline import Pipeline
from sklearn.manifold import TSNE
from sklearn.decomposition import PCA
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import roc_auc_score, f1_score, accuracy_score
from sklearn.model_selection import train_test_split

ALL_NUMERIC = ['patient_info.age', 'biomarkers.tmb']
ALL_CATEGORICAL = ['patient_info.gender', 'biomarkers.msi']


def detect_fields(real_df):
    """字段自适应: 仅保留真实数据中覆盖率>50%的字段"""
    numeric = [c for c in ALL_NUMERIC
               if c in real_df.columns and pd.to_numeric(real_df[c], errors='coerce').notna().mean() > 0.5]
    categorical = [c for c in ALL_CATEGORICAL
                   if c in real_df.columns and real_df[c].notna().mean() > 0.5]
    return numeric, categorical


def add_mut_columns(df, genes):
    for gene in genes:
        df[f'Mut_{gene}'] = df['somatic_variants'].apply(lambda s: int(
            any(isinstance(v, dict) and v.get('gene') == gene
                for v in safe_parse(s))))
    return df


def safe_parse(val):
    try:
        if pd.isna(val) or str(val) == 'nan' or not val:
            return []
        if isinstance(val, str):
            return ast.literal_eval(val)
        return val
    except Exception:
        return []


# ══════════════════════════ stats ══════════════════════════
def run_stats(real_df, syn_df, eval_dir, numeric, categorical):
    print("[stats] KS / TVD / gene-freq / corr-MAE / PCA-tSNE ...")
    metrics = {"Univariate_Marginal_Similarity": {},
               "Bivariate_Correlation": {}, "Somatic_Variant_Fidelity": {}}

    for col in numeric:
        r = pd.to_numeric(real_df[col], errors='coerce').dropna()
        s = pd.to_numeric(syn_df[col], errors='coerce').dropna()
        if len(r) and len(s):
            stat, p = ks_2samp(r, s)
            metrics["Univariate_Marginal_Similarity"][col] = {
                "KS_statistic": round(float(stat), 4), "p_value": round(float(p), 4)}

    for col in categorical:
        r = real_df[col].fillna("Missing").astype(str)
        s = syn_df[col].fillna("Missing").astype(str)
        metrics["Univariate_Marginal_Similarity"][col] = {"TVD": round(compute_tvd(r, s), 4)}

    real_freq = get_gene_frequencies(real_df)
    syn_freq = get_gene_frequencies(syn_df)
    top30 = sorted(real_freq, key=lambda g: real_freq[g], reverse=True)[:30]
    diffs = {g: abs(real_freq.get(g, 0.0) - syn_freq.get(g, 0.0)) for g in top30}
    metrics["Somatic_Variant_Fidelity"] = {
        "Average_Gene_Frequency_Difference": round(float(np.mean(list(diffs.values()))), 4) if diffs else 0.0,
        "Top_Gene_Differences": {g: {"Real_Frequency": round(real_freq.get(g, 0.0), 4),
                                     "Synthetic_Frequency": round(syn_freq.get(g, 0.0), 4),
                                     "Frequency_Difference_TVD": round(d, 4)} for g, d in diffs.items()}}

    if len(numeric) >= 2:
        rc = real_df[numeric].apply(pd.to_numeric, errors='coerce').corr().fillna(0).values
        sc = syn_df[numeric].apply(pd.to_numeric, errors='coerce').corr().fillna(0).values
        metrics["Bivariate_Correlation"]["Pearson_MAE"] = round(float(np.mean(np.abs(rc - sc))), 4)
    else:
        metrics["Bivariate_Correlation"]["Pearson_MAE"] = None

    # PCA / t-SNE
    print("[stats] PCA/t-SNE (syn 下采样至 2000) ...")
    top15 = sorted(real_freq, key=lambda g: real_freq[g], reverse=True)[:15]
    real_v, syn_v = real_df.copy(), syn_df.copy()
    real_v['DataType'], syn_v['DataType'] = 'Real', 'Synthetic'
    if len(syn_v) > 2000:
        syn_v = syn_v.sample(n=2000, random_state=42)
    combined = pd.concat([real_v, syn_v], ignore_index=True)
    combined = add_mut_columns(combined, top15)
    bin_cols = [f'Mut_{g}' for g in top15]
    feats = combined[numeric + categorical + bin_cols].copy()
    for c in numeric:
        feats[c] = pd.to_numeric(feats[c], errors='coerce').fillna(0)
    for c in categorical:
        feats[c] = feats[c].fillna('Missing').astype(str)

    transformers = []
    if numeric:
        transformers.append(('num', StandardScaler(), numeric))
    if categorical:
        transformers.append(('cat', OneHotEncoder(handle_unknown='ignore', sparse_output=False), categorical))
    transformers.append(('bin', 'passthrough', bin_cols))
    X = ColumnTransformer(transformers=transformers).fit_transform(feats)

    Xp = PCA(n_components=2, random_state=42).fit_transform(X)
    Xts = TSNE(n_components=2, perplexity=min(30, len(combined) - 1), random_state=42).fit_transform(X)
    plt.figure(figsize=(14, 6))
    for i, (emb, name) in enumerate([(Xp, 'PCA'), (Xts, 't-SNE')], 1):
        plt.subplot(1, 2, i)
        sns.scatterplot(x=emb[:, 0], y=emb[:, 1], hue=combined['DataType'], alpha=0.6,
                        palette={'Real': 'blue', 'Synthetic': 'orange'})
        plt.title(f'{name}: Real vs Synthetic')
    plt.tight_layout()
    plot_path = os.path.join(eval_dir, 'manifold_overlap.png')
    plt.savefig(plot_path, dpi=300); plt.close()

    with open(os.path.join(eval_dir, 'statistical_fidelity_metrics.json'), 'w') as f:
        json.dump(metrics, f, indent=4)
    print(f"[stats] ✅ -> {eval_dir}/statistical_fidelity_metrics.json + manifold_overlap.png")
    print(json.dumps({k: v for k, v in metrics.items() if k != "Somatic_Variant_Fidelity"}, indent=2))


# ══════════════════════════ tstr ══════════════════════════
def run_tstr(real_df, syn_df, eval_dir, numeric, categorical):
    if 'biomarkers.tmb' not in numeric:
        print("[tstr] ⚠ TMB 不可用, 跳过")
        return
    from shap import TreeExplainer
    print("[tstr] 目标: TMB>=10 | 特征: Mut_top10 + 分类/数值字段 ...")
    # 关键: biomarkers.tmb 是目标本身, 绝不能进特征 (防泄漏)
    numeric_safe = [c for c in numeric if c != 'biomarkers.tmb']
    for df in (real_df, syn_df):
        df['TMB_High'] = (pd.to_numeric(df['biomarkers.tmb'], errors='coerce').fillna(0) >= 10).astype(int)
    _freqs = get_gene_frequencies(real_df)  # 缓存: 避免排序时对每基因重复解析全部行
    top10 = sorted(_freqs, key=_freqs.get, reverse=True)[:10]
    real_df, syn_df = add_mut_columns(real_df, top10), add_mut_columns(syn_df, top10)
    bin_cols = [f'Mut_{g}' for g in top10]
    feat_cols = numeric_safe + categorical + bin_cols

    for df in (real_df, syn_df):
        for c in numeric_safe:
            df[c] = pd.to_numeric(df[c], errors='coerce').fillna(df[c].median() if df[c].notna().any() else 0)
        for c in categorical:
            df[c] = df[c].fillna('Missing').astype(str)

    Xr, yr = real_df[feat_cols], real_df['TMB_High']
    Xs, ys = syn_df[feat_cols], syn_df['TMB_High']

    transformers = []
    if numeric_safe:
        transformers.append(('num', StandardScaler(), numeric_safe))
    if categorical:
        transformers.append(('cat', OneHotEncoder(handle_unknown='ignore', sparse_output=False), categorical))
    transformers.append(('bin', 'passthrough', bin_cols))
    pre = ColumnTransformer(transformers=transformers)
    pre.fit(pd.concat([Xr, Xs], ignore_index=True))

    # 标准协议: 真实数据 70/30 切分, 在真实测试集上评估 TRTR/TSTR (避免训练集记忆饱和)
    Xr_tr, Xr_te, yr_tr, yr_te = train_test_split(Xr, yr, test_size=0.3,
                                                  random_state=42, stratify=yr)
    Xrp_tr, Xrp_te = pre.transform(Xr_tr), pre.transform(Xr_te)
    Xsp = pre.transform(Xs)

    clf_trtr = RandomForestClassifier(random_state=42, n_estimators=100).fit(Xrp_tr, yr_tr)
    clf_tstr = RandomForestClassifier(random_state=42, n_estimators=100).fit(Xsp, ys)
    try:
        auc_trtr = roc_auc_score(yr_te, clf_trtr.predict_proba(Xrp_te)[:, 1])
        auc_tstr = roc_auc_score(yr_te, clf_tstr.predict_proba(Xrp_te)[:, 1])
    except ValueError:
        auc_trtr = auc_tstr = 0.5
    f1_trtr = f1_score(yr_te, clf_trtr.predict(Xrp_te), zero_division=0)
    f1_tstr = f1_score(yr_te, clf_tstr.predict(Xrp_te), zero_division=0)

    # SHAP 特征重要性秩相关 (TRTR 在真实训练集, TSTR 在合成集)
    # SHAP_SAMPLE 环境变量可对大矩阵抽样 (默认 0=全量, 与历史行为一致)
    import os as _os
    _ss = int(_os.environ.get('SHAP_SAMPLE', '0'))
    if _ss and _ss < len(Xrp_tr):
        Xrp_sh = shap.sample(Xrp_tr, _ss, random_state=42)
    else:
        Xrp_sh = Xrp_tr
    if _ss and _ss < len(Xsp):
        Xsp_sh = shap.sample(Xsp, _ss, random_state=42)
    else:
        Xsp_sh = Xsp
    print(f"[tstr] SHAP importance ... (样本数: {len(Xrp_sh)}/{len(Xsp_sh)})")
    ex_r = TreeExplainer(clf_trtr).shap_values(Xrp_sh)
    ex_s = TreeExplainer(clf_tstr).shap_values(Xsp_sh)
    def cls1(sv):
        if isinstance(sv, list): return sv[1]
        if hasattr(sv, 'shape') and len(sv.shape) == 3: return sv[:, :, 1]
        return sv
    imp_r = np.abs(cls1(ex_r)).mean(axis=0)
    imp_s = np.abs(cls1(ex_s)).mean(axis=0)
    res = spearmanr(imp_r, imp_s)
    rho = float(res.statistic) if hasattr(res, 'statistic') else float(res[0])
    p = float(res.pvalue) if hasattr(res, 'pvalue') else float(res[1])

    results = {"Downstream_Task": "Predict TMB >= 10",
               "Features": feat_cols,
               "Performance": {"TRTR": {"AUROC": round(auc_trtr, 4), "F1": round(f1_trtr, 4)},
                               "TSTR": {"AUROC": round(auc_tstr, 4), "F1": round(f1_tstr, 4)},
                               "Delta_AUROC": round(abs(auc_trtr - auc_tstr), 4)},
               "Explainability": {"SHAP_Rank_Correlation": round(rho, 4), "p_value": round(p, 4)}}
    with open(os.path.join(eval_dir, 'utility_tstr_metrics.json'), 'w') as f:
        json.dump(results, f, indent=4)
    print(f"[tstr] ✅ -> {eval_dir}/utility_tstr_metrics.json")
    print(json.dumps(results, indent=2))


# ══════════════════════════ privacy ══════════════════════════
def run_privacy(real_df, syn_df, eval_dir, numeric, categorical):
    print("[privacy] exact-match / DCR-NPS / MIA ...")
    top10 = sorted(get_gene_frequencies(real_df), key=lambda g: get_gene_frequencies(real_df)[g],
                   reverse=True)[:10]
    real_df, syn_df = add_mut_columns(real_df, top10), add_mut_columns(syn_df, top10)
    bin_cols = [f'Mut_{g}' for g in top10]
    feat_cols = numeric + categorical + bin_cols

    for df in (real_df, syn_df):
        for c in numeric:
            df[c] = pd.to_numeric(df[c], errors='coerce').fillna(df[c].median() if df[c].notna().any() else 0)
        for c in categorical:
            df[c] = df[c].fillna('Missing').astype(str)

    # 1) 记录级完全匹配 (含原始10列) 与特征级匹配
    def exact(df_a, cols):
        ta = set(map(tuple, df_a[cols].to_numpy()))
        return ta
    raw_cols = [c for c in real_df.columns if c != 'somatic_variants' and c in syn_df.columns]
    m_raw = len(exact(real_df, raw_cols) & exact(syn_df, raw_cols))
    m_feat = len(exact(real_df, feat_cols) & exact(syn_df, feat_cols))

    # 2) DCR / NPS
    transformers = []
    if numeric:
        transformers.append(('num', StandardScaler(), numeric))
    if categorical:
        transformers.append(('cat', OneHotEncoder(handle_unknown='ignore', sparse_output=False), categorical))
    transformers.append(('bin', 'passthrough', bin_cols))
    pre = ColumnTransformer(transformers=transformers)
    pre.fit(pd.concat([real_df[feat_cols], syn_df[feat_cols]], ignore_index=True))
    dcr = compute_dcr_and_nps(pre.transform(real_df[feat_cols]), pre.transform(syn_df[feat_cols]))

    # 3) MIA
    mia_acc, mia_f1 = membership_inference_attack(real_df[feat_cols], syn_df[feat_cols], pre)

    results = {"Privacy_Assessment": {
        "Exact_Match": {"Record_Level_Matches": m_raw,
                        "Feature_Level_Matches": m_feat,
                        "Total_Synthetic_Records": len(syn_df),
                        "Status": "PASS" if m_raw == 0 else "FAIL"},
        "Distance_to_Closest_Record": dcr,
        "Membership_Inference_Attack": {
            "Attacker_Accuracy": round(float(mia_acc), 4),
            "Attacker_F1_Score": round(float(mia_f1), 4),
            "Baseline_Random_Guess_Accuracy": 0.5,
            "Status": "PASS" if mia_acc < 0.65 else "WARNING"}}}
    with open(os.path.join(eval_dir, 'privacy_metrics.json'), 'w') as f:
        json.dump(results, f, indent=4)
    print(f"[privacy] ✅ -> {eval_dir}/privacy_metrics.json")
    print(json.dumps(results, indent=2))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('task', choices=['stats', 'tstr', 'privacy', 'all'])
    ap.add_argument('--real', required=True)
    ap.add_argument('--syn', required=True)
    ap.add_argument('--eval-dir', required=True)
    args = ap.parse_args()
    os.makedirs(args.eval_dir, exist_ok=True)

    print(f"加载 real={args.real}\n      syn ={args.syn}")
    real_df, syn_df = load_data(args.real, args.syn)
    numeric, categorical = detect_fields(real_df)
    print(f"自适应字段: numeric={numeric}, categorical={categorical}")

    tasks = ['stats', 'tstr', 'privacy'] if args.task == 'all' else [args.task]
    for t in tasks:
        {'stats': run_stats, 'tstr': run_tstr, 'privacy': run_privacy}[t](
            real_df.copy(), syn_df.copy(), args.eval_dir, numeric, categorical)


if __name__ == "__main__":
    main()

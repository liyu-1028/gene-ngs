"""
Real-vs-synthetic discriminability control experiment.

Context: privacy_metrics.json reports "Membership_Inference_Attack" accuracy 0.679
(> 0.5 random-guess baseline). This experiment determines whether that signal comes
from (a) generator memorization/leakage, or (b) the metric's intrinsic background /
generator fidelity gap.

Controls (identical preprocessor & attacker as privacy_assessment.py):
  R0  real vs synthetic   -- reproduces the published 0.679 (anchor)
  C1  realA vs realB      -- disjoint halves of the REAL cohort (intrinsic background)
  C2  synA vs synB        -- attacker calibration (~0.5 expected)
  C3  C1 with 3 extra seeds -- background stability

Findings (2026-09-18): C1/C3 ~ 0.49-0.51 across seeds -> the 0.68 signal is not a
dataset artifact. Combined with Exact_Match=0 and DCR leakage=0.0 (both PASS), the
evidence supports the interpretation that the discriminability proxy quantifies the
generator fidelity gap, not membership leakage.

Usage:  python src/stage5/mia_control_experiment.py   (paths configurable below)
"""
import sys, os, json
import numpy as np
import pandas as pd
from sklearn.preprocessing import StandardScaler, OneHotEncoder
from sklearn.compose import ColumnTransformer

sys.path.insert(0, "/home/lyl/gene-ngs")
from src.stage5.privacy_assessment import prepare_data, membership_inference_attack

REAL = "/home/lyl/gene-ngs/data/real/msk_impact_2017/real_dataset.csv"
SYN  = "/home/lyl/gene-ngs/data/synthetic_2017/audited_synthetic_data.csv"

def make_preprocessor(df_a, df_b, top_genes):
    pre = ColumnTransformer(transformers=[
        ('num', StandardScaler(), ['patient_info.age', 'biomarkers.tmb']),
        ('cat', OneHotEncoder(handle_unknown='ignore', sparse_output=False), ['patient_info.gender', 'biomarkers.msi']),
        ('bin', 'passthrough', [f'Mut_{g}' for g in top_genes])])
    pre.fit(pd.concat([df_a, df_b], ignore_index=True))
    return pre

def half(df, seed=42):
    idx = np.random.RandomState(seed).permutation(len(df))
    h = len(idx) // 2
    return df.iloc[idx[:h]].reset_index(drop=True), df.iloc[idx[h:]].reset_index(drop=True)

real_df, syn_df, top_genes = prepare_data(REAL, SYN)
print(f"真实 {len(real_df)} 行, 合成 {len(syn_df)} 行, top_genes={len(top_genes)}")
results = {}

pre = make_preprocessor(real_df, syn_df, top_genes)
acc, f1 = membership_inference_attack(real_df, syn_df, pre)
results["R0_real_vs_syn"] = {"acc": round(acc,4), "f1": round(f1,4)}
print(f"R0 真实vs合成    : acc={acc:.4f} f1={f1:.4f}  (锚点, 应≈0.679)")

rA, rB = half(real_df, 42)
pre1 = make_preprocessor(rA, rB, top_genes)
acc, f1 = membership_inference_attack(rA, rB, pre1)
results["C1_realA_vs_realB"] = {"acc": round(acc,4), "f1": round(f1,4)}
print(f"C1 真实A vs真实B : acc={acc:.4f} f1={f1:.4f}  ← 数据集本底")

sA, sB = half(syn_df, 42)
pre2 = make_preprocessor(sA, sB, top_genes)
acc, f1 = membership_inference_attack(sA, sB, pre2)
results["C2_synA_vs_synB"] = {"acc": round(acc,4), "f1": round(f1,4)}
print(f"C2 合成A vs合成B : acc={acc:.4f} f1={f1:.4f}  (应≈0.5)")

accs = []
for sd in (7, 13, 99):
    rA, rB = half(real_df, sd)
    pre3 = make_preprocessor(rA, rB, top_genes)
    a, f = membership_inference_attack(rA, rB, pre3)
    accs.append(round(a,4))
results["C3_real_halves_multiseed"] = accs
print(f"C3 真实对半×3种子: {accs}")

out = os.environ.get("MIA_CONTROL_OUT", "evaluation/mia_control_results.json")
json.dump(results, open(out, "w"), indent=2)
print(f"已保存 {out}")

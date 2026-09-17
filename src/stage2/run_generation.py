# --- src/stage2/run_generation.py (v2, MSK-IMPACT 50K 量级) ---
"""
混合生成 Stage 2:
  1) CTGAN 学习真实临床特征的联合分布 (年龄/性别/诊断/样本类型/TMB/MSI/QC)
  2) Offspring Enhancer 通过双亲变异重组生成体细胞变异组合

v2 变更:
  - 输入默认切换为 MSK-IMPACT 50K 真实数据集 (54,331 样本)
  - 合成规模 N_SYNTHETIC 默认与真实数据 1:1 (原为硬编码 1000)
  - CTGAN 仅在关键字段完整的样本上训练 (避免 NaN 训练偏置), 生成仍为全量 N
"""
import pandas as pd
import ast
import logging
import os
import sys
import argparse
import numpy as np

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '../..')))
from src.stage2.ctgan_generator import ClinicalTabularGenerator

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

DEFAULT_INPUT = "data/real/msk_impact_50k_2026/real_dataset.csv"
OUTPUT_DIR = "data/synthetic"


def run_stage2_generation(input_csv: str = DEFAULT_INPUT,
                          n_synthetic: int = None,
                          seed: int = 42,
                          output_dir: str = OUTPUT_DIR):
    np.random.seed(seed)
    os.makedirs(output_dir, exist_ok=True)

    if not os.path.exists(input_csv):
        logging.error(f"Input file not found: {input_csv}")
        return

    df = pd.read_csv(input_csv)
    n_target = n_synthetic or len(df)          # 默认与真实数据量级 1:1
    logging.info(f"Loaded {len(df):,} real records from {input_csv}. "
                 f"Synthetic target = {n_target:,}")

    # 1. 分离临床特征与基因组特征
    clinical_cols = [c for c in df.columns if c not in ['somatic_variants', 'source_file', 'qc_metrics']]
    clinical_df = df[clinical_cols]
    variant_df = df[['somatic_variants']]

    # 2. CTGAN: 仅在实际有值且覆盖率>50%的关键字段上要求完整性 (自适应 2017/50K 源)
    key_fields = ['patient_info.age', 'biomarkers.tmb', 'patient_info.gender']
    usable = [c for c in key_fields
              if c in clinical_df.columns and clinical_df[c].notna().mean() > 0.5]
    train_df = clinical_df.dropna(subset=usable)
    # 全空列 (如 2017 源的 age) 不进 CTGAN，生成后补回保持列模式
    all_null_cols = [c for c in train_df.columns if train_df[c].isna().all()]
    train_df = train_df.drop(columns=all_null_cols)
    logging.info(f"CTGAN training subset (complete on {usable or 'no key fields'}): "
                 f"{len(train_df):,} ({100*len(train_df)/len(clinical_df):.1f}%); "
                 f"excluded all-null cols: {all_null_cols}")

    logging.info("Step 1: Training CTGAN for clinical features...")
    tabular_gen = ClinicalTabularGenerator()
    tabular_gen.fit(train_df)
    synthetic_clinical = tabular_gen.sample(n_target)
    for c in all_null_cols:
        synthetic_clinical[c] = pd.NA
    synthetic_clinical = synthetic_clinical[[c for c in clinical_df.columns]]

    # 3. Offspring Enhancer (双亲变异交叉重组)
    logging.info("Step 2: Running Offspring Enhancer for variants...")
    variant_df = variant_df.copy()
    variant_df['somatic_variants_obj'] = variant_df['somatic_variants'].apply(ast.literal_eval)

    real_variants = variant_df['somatic_variants_obj'].tolist()

    def hybridize_variants(v1, v2):
        all_v = v1 + v2
        if not all_v:
            return "[]"
        sample_size = np.random.randint(len(all_v) // 2 + 1, len(all_v) + 1)
        unique_v = {f"{v['gene']}_{v['alteration']}": v for v in all_v}.values()
        child_v = np.random.choice(list(unique_v), min(sample_size, len(unique_v)), replace=False)
        return str(list(child_v))

    synthetic_variants = []
    for i in range(n_target):
        idx1, idx2 = np.random.choice(len(real_variants), 2)
        synthetic_variants.append(hybridize_variants(real_variants[idx1], real_variants[idx2]))
        if (i + 1) % 5000 == 0:
            logging.info(f"  Offspring progress: {i+1:,}/{n_target:,}")

    synthetic_clinical['somatic_variants'] = synthetic_variants

    # 4. 保存
    output_path = os.path.join(output_dir, "raw_synthetic_data.csv")
    synthetic_clinical.to_csv(output_path, index=False, encoding="utf-8")
    logging.info(f"Hybrid Generation complete! {n_target:,} records saved to {output_path}")


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument('--input', default=DEFAULT_INPUT)
    ap.add_argument('--n', type=int, default=None,
                    help='合成样本数; 默认与真实数据集 1:1')
    ap.add_argument('--seed', type=int, default=42)
    ap.add_argument('--output-dir', default=OUTPUT_DIR)
    args = ap.parse_args()
    run_stage2_generation(args.input, args.n, args.seed, args.output_dir)

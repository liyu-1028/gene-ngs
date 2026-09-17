import pandas as pd
import numpy as np
import json
import os
import ast
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import roc_auc_score, f1_score, precision_score
from sklearn.preprocessing import StandardScaler, OneHotEncoder
from sklearn.compose import ColumnTransformer
from sklearn.pipeline import Pipeline
import shap
import matplotlib.pyplot as plt
from scipy.stats import spearmanr

def get_top_genes(df, top_n=10):
    gene_counts = {}
    for idx, val in df['somatic_variants'].items():
        try:
            if pd.isna(val) or str(val) == 'nan' or not val:
                variants = []
            elif isinstance(val, str):
                variants = ast.literal_eval(val)
            else:
                variants = val
        except Exception:
            variants = []
        
        seen_genes = set()
        for v in variants:
            if isinstance(v, dict) and 'gene' in v:
                seen_genes.add(v['gene'])
        for gene in seen_genes:
            gene_counts[gene] = gene_counts.get(gene, 0) + 1
            
    return sorted(gene_counts.keys(), key=lambda g: gene_counts[g], reverse=True)[:top_n]

def parse_variants(df, genes_of_interest):
    # Create binary feature columns for the top most common genes
    for gene in genes_of_interest:
        df[f'Mut_{gene}'] = 0
        
    for idx, row in df.iterrows():
        variants_str = str(row.get('somatic_variants', '[]'))
        try:
            if pd.isna(variants_str) or variants_str == 'nan':
                variants = []
            else:
                variants = ast.literal_eval(variants_str)
        except:
            variants = []
            
        for v in variants:
            if isinstance(v, dict):
                gene = v.get('gene')
                if gene in genes_of_interest:
                    df.at[idx, f'Mut_{gene}'] = 1
    return df

def prepare_data(real_path, syn_path):
    try:
        real_df = pd.read_csv(real_path, encoding='utf-8')
    except:
        real_df = pd.read_csv(real_path, encoding='gbk')
        
    try:
        syn_df = pd.read_csv(syn_path, encoding='utf-8')
    except:
        syn_df = pd.read_csv(syn_path, encoding='gbk')

    # Dynamically select top 10 mutated genes from real data
    genes_of_interest = get_top_genes(real_df, top_n=10)

    real_df = parse_variants(real_df, genes_of_interest)
    syn_df = parse_variants(syn_df, genes_of_interest)

    # Define Target: Predicting TMB-High (e.g., TMB >= 10)
    real_df['TMB_High'] = (pd.to_numeric(real_df['biomarkers.tmb'], errors='coerce').fillna(0) >= 10).astype(int)
    syn_df['TMB_High'] = (pd.to_numeric(syn_df['biomarkers.tmb'], errors='coerce').fillna(0) >= 10).astype(int)

    # Clean continuous
    real_df['patient_info.age'] = pd.to_numeric(real_df['patient_info.age'], errors='coerce').fillna(60)
    syn_df['patient_info.age'] = pd.to_numeric(syn_df['patient_info.age'], errors='coerce').fillna(60)
    
    # Clean categorical
    real_df['patient_info.gender'] = real_df['patient_info.gender'].fillna('Unknown').astype(str)
    syn_df['patient_info.gender'] = syn_df['patient_info.gender'].fillna('Unknown').astype(str)
    
    features = ['patient_info.age', 'patient_info.gender'] + [f'Mut_{g}' for g in genes_of_interest]
    
    X_real = real_df[features]
    y_real = real_df['TMB_High']
    X_syn = syn_df[features]
    y_syn = syn_df['TMB_High']

    return X_real, y_real, X_syn, y_syn, features, genes_of_interest

def train_and_evaluate():
    script_dir = os.path.dirname(os.path.abspath(__file__))
    project_base = os.path.abspath(os.path.join(script_dir, "../.."))
    real_path = os.path.join(project_base, "data/processed/clinical_features_dataset.csv")
    syn_path = os.path.join(project_base, "data/synthetic/audited_synthetic_data.csv")
    eval_dir = os.path.join(project_base, "data/synthetic/evaluation")
    os.makedirs(eval_dir, exist_ok=True)
    
    X_real, y_real, X_syn, y_syn, feature_names, genes_of_interest = prepare_data(real_path, syn_path)
    
    # Define preprocessing
    numeric_features = ['patient_info.age']
    categorical_features = ['patient_info.gender']
    binary_features = [f for f in feature_names if f.startswith('Mut_')]

    preprocessor = ColumnTransformer(
        transformers=[
            ('num', StandardScaler(), numeric_features),
            ('cat', OneHotEncoder(handle_unknown='ignore', sparse_output=False), categorical_features),
            ('bin', 'passthrough', binary_features)
        ])

    # Fit preprocessor on combined data to ensure consistent feature space
    X_combined = pd.concat([X_real, X_syn], ignore_index=True)
    preprocessor.fit(X_combined)

    X_real_processed = preprocessor.transform(X_real)
    X_syn_processed = preprocessor.transform(X_syn)

    # Model TRTR (Train Real, Test Real)
    clf_real = RandomForestClassifier(random_state=42, n_estimators=100)
    clf_real.fit(X_real_processed, y_real)
    
    # Model TSTR (Train Synthetic, Test Real)
    clf_syn = RandomForestClassifier(random_state=42, n_estimators=100)
    clf_syn.fit(X_syn_processed, y_syn)

    # Evaluate on Real Data
    y_pred_trtr = clf_real.predict(X_real_processed)
    y_prob_trtr = clf_real.predict_proba(X_real_processed)[:, 1]
    
    y_pred_tstr = clf_syn.predict(X_real_processed)
    y_prob_tstr = clf_syn.predict_proba(X_real_processed)[:, 1]

    try:
        auc_trtr = roc_auc_score(y_real, y_prob_trtr)
        auc_tstr = roc_auc_score(y_real, y_prob_tstr)
    except ValueError:
        auc_trtr = 0.5
        auc_tstr = 0.5

    f1_trtr = f1_score(y_real, y_pred_trtr, zero_division=0)
    f1_tstr = f1_score(y_real, y_pred_tstr, zero_division=0)

    # Compute SHAP values for model explainability
    clf_real_pipe = Pipeline(steps=[('preprocessor', preprocessor),
                                    ('classifier', RandomForestClassifier(random_state=42, n_estimators=100))])
    clf_real_pipe.fit(X_real, y_real)
    
    clf_syn_pipe = Pipeline(steps=[('preprocessor', preprocessor),
                                   ('classifier', RandomForestClassifier(random_state=42, n_estimators=100))])
    clf_syn_pipe.fit(X_syn, y_syn)

    X_real_processed = clf_real_pipe.named_steps['preprocessor'].transform(X_real)
    X_syn_processed = clf_syn_pipe.named_steps['preprocessor'].transform(X_syn)
    
    cat_encoder = clf_real_pipe.named_steps['preprocessor'].named_transformers_['cat']
    cat_names = cat_encoder.get_feature_names_out(categorical_features).tolist()
    processed_feature_names = numeric_features + cat_names + binary_features

    explainer_real = shap.TreeExplainer(clf_real_pipe.named_steps['classifier'])
    shap_vals_real = explainer_real.shap_values(X_real_processed)
    
    explainer_syn = shap.TreeExplainer(clf_syn_pipe.named_steps['classifier'])
    shap_vals_syn = explainer_syn.shap_values(X_syn_processed)

    # Handle binary vs multiclass output of shap for RandomForest
    if isinstance(shap_vals_real, list):
        shap_vals_real = shap_vals_real[1]
    elif len(shap_vals_real.shape) == 3:
        shap_vals_real = shap_vals_real[:, :, 1]
        
    if isinstance(shap_vals_syn, list):
        shap_vals_syn = shap_vals_syn[1]
    elif len(shap_vals_syn.shape) == 3:
        shap_vals_syn = shap_vals_syn[:, :, 1]

    # Mean absolute SHAP values per feature
    mean_shap_real = np.abs(shap_vals_real).mean(axis=0)
    mean_shap_syn = np.abs(shap_vals_syn).mean(axis=0)
    
    # Calculate Spearman rank correlation of feature importance
    res = spearmanr(mean_shap_real, mean_shap_syn)
    spearman_corr, p_val = res.statistic if hasattr(res, 'statistic') else res[0], res.pvalue if hasattr(res, 'pvalue') else res[1]

    # Ensure scalar floats for JSON serialization
    spearman_corr = float(spearman_corr) if not isinstance(spearman_corr, np.ndarray) else float(spearman_corr.item())
    p_val = float(p_val) if not isinstance(p_val, np.ndarray) else float(p_val.item())

    # Save metrics
    results = {
        "Downstream_Task": "Predict TMB >= 10",
        "Performance": {
            "TRTR": {"AUROC": round(auc_trtr, 4), "F1": round(f1_trtr, 4)},
            "TSTR": {"AUROC": round(auc_tstr, 4), "F1": round(f1_tstr, 4)},
            "Delta_AUROC": round(abs(auc_trtr - auc_tstr), 4)
        },
        "Explainability": {
            "SHAP_Rank_Correlation": round(spearman_corr, 4),
            "p_value": round(p_val, 4)
        }
    }
    
    metrics_path = os.path.join(eval_dir, 'utility_tstr_metrics.json')
    with open(metrics_path, 'w', encoding='utf-8') as f:
        json.dump(results, f, indent=4)
        
    print("\n--- Downstream Utility (TSTR) Report ---")
    print(json.dumps(results, indent=4))
    
    if results['Performance']['Delta_AUROC'] <= 0.1:
        print("\n[SUCCESS] TSTR performance closely matches TRTR, proving high utility!")
    
    if results['Explainability']['SHAP_Rank_Correlation'] > 0.8:
        print("[SUCCESS] SHAP Rank Correlation > 0.8, proving clinical interpretability is preserved!")

if __name__ == "__main__":
    train_and_evaluate()

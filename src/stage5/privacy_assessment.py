import pandas as pd
import numpy as np
import json
import os
import ast
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score, f1_score
from sklearn.preprocessing import StandardScaler, OneHotEncoder
from sklearn.compose import ColumnTransformer
from sklearn.neighbors import NearestNeighbors

def get_gene_frequencies(df):
    gene_counts = {}
    total_samples = len(df)
    if total_samples == 0:
        return gene_counts
    
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
            
    return {gene: count / total_samples for gene, count in gene_counts.items()}

def prepare_data(real_path, syn_path):
    try:
        real_df = pd.read_csv(real_path, encoding='utf-8')
    except:
        real_df = pd.read_csv(real_path, encoding='gbk')
        
    try:
        syn_df = pd.read_csv(syn_path, encoding='utf-8')
    except:
        syn_df = pd.read_csv(syn_path, encoding='gbk')

    # Extract top 30 genes from real data to construct binary mutation features
    real_gene_freqs = get_gene_frequencies(real_df)
    top_genes = sorted(real_gene_freqs.keys(), key=lambda g: real_gene_freqs[g], reverse=True)[:30]

    # Pre-clean continuous & categorical columns
    for df in [real_df, syn_df]:
        df['patient_info.age'] = pd.to_numeric(df['patient_info.age'], errors='coerce').fillna(60)
        df['biomarkers.tmb'] = pd.to_numeric(df['biomarkers.tmb'], errors='coerce').fillna(2.0)
        df['patient_info.gender'] = df['patient_info.gender'].fillna('Unknown').astype(str)
        df['biomarkers.msi'] = df['biomarkers.msi'].fillna('MSS').astype(str)
        
        # Binarize mutations for the top genes
        for gene in top_genes:
            col_name = f'Mut_{gene}'
            df[col_name] = 0
            for idx, row in df.iterrows():
                val = row.get('somatic_variants', '[]')
                try:
                    if pd.isna(val) or str(val) == 'nan' or not val:
                        variants = []
                    elif isinstance(val, str):
                        variants = ast.literal_eval(val)
                    else:
                        variants = val
                except:
                    variants = []
                if any(isinstance(v, dict) and v.get('gene') == gene for v in variants):
                    df.at[idx, col_name] = 1
                    
    feature_cols = ['patient_info.age', 'patient_info.gender', 'biomarkers.tmb', 'biomarkers.msi'] + [f'Mut_{g}' for g in top_genes]
    
    return real_df[feature_cols].copy(), syn_df[feature_cols].copy(), top_genes

def compute_exact_match(real_df, syn_df):
    """
    Check if any synthetic record is an exact duplicate of a real record.
    Returns the count of exact matches.
    """
    real_tuples = set([tuple(x) for x in real_df.to_numpy()])
    syn_tuples = set([tuple(x) for x in syn_df.to_numpy()])
    
    matches = real_tuples.intersection(syn_tuples)
    return len(matches), len(syn_df)

def compute_dcr_and_nps(X_real, X_syn):
    """
    Computes Distance to Closest Record (DCR) and Nearest Neighbor Privacy Score (NPS).
    """
    # 1. Cross-set DCR (Synthetic to Real)
    nn_cross = NearestNeighbors(n_neighbors=1, metric='euclidean')
    nn_cross.fit(X_real)
    distances_cross, _ = nn_cross.kneighbors(X_syn)
    dcr_cross = distances_cross.flatten()
    
    # 2. Intra-set DCR (Real to other Real)
    nn_intra = NearestNeighbors(n_neighbors=min(2, len(X_real)), metric='euclidean')
    nn_intra.fit(X_real)
    distances_intra, _ = nn_intra.kneighbors(X_real)
    
    if distances_intra.shape[1] > 1:
        dcr_intra = distances_intra[:, 1]
    else:
        dcr_intra = distances_intra[:, 0]
        
    # Calculate NPS & Leakage Rate
    # Privacy Leakage Rate is the proportion of synthetic records whose DCR to real data 
    # is smaller than the 5th percentile of real-to-real nearest neighbor distances.
    threshold_5th = np.percentile(dcr_intra, 5) if len(dcr_intra) > 0 else 0.0
    leakage_rate = np.mean(dcr_cross < threshold_5th) if len(dcr_cross) > 0 else 0.0
    
    dcr_results = {
        "Median_Cross_Set_DCR": round(float(np.median(dcr_cross)), 4),
        "Median_Intra_Set_DCR": round(float(np.median(dcr_intra)), 4),
        "Min_Cross_Set_DCR": round(float(np.min(dcr_cross)), 4),
        "Privacy_Threshold_5th_Percentile": round(float(threshold_5th), 4),
        "Privacy_Leakage_Rate": round(float(leakage_rate), 4),
        "Status": "PASS" if leakage_rate < 0.05 else "WARNING"
    }
    
    return dcr_results

def membership_inference_attack(real_df, syn_df, preprocessor):
    """
    Simulate a Membership Inference Attack (MIA).
    """
    # Transform data
    X_real = preprocessor.transform(real_df)
    X_syn = preprocessor.transform(syn_df)
    
    # Setup the Attack Scenario
    # Label: 1 = Real, 0 = Synthetic (proxy for holdout)
    real_target = np.ones(len(X_real))
    syn_target = np.zeros(len(X_syn))
    
    # Handle class imbalance if synthetic data is large
    # Attacker trains on a balanced subset of members and non-members
    sample_size = min(len(X_real), len(X_syn))
    
    real_indices = np.random.choice(len(X_real), sample_size, replace=False)
    syn_indices = np.random.choice(len(X_syn), sample_size, replace=False)
    
    X_attack = np.vstack([X_real[real_indices], X_syn[syn_indices]])
    y_attack = np.hstack([real_target[real_indices], syn_target[syn_indices]])
    
    # Train/Test split for the attacker
    X_train, X_test, y_train, y_test = train_test_split(X_attack, y_attack, test_size=0.5, random_state=42, stratify=y_attack)
    
    # Train Attacker Model
    attacker = RandomForestClassifier(n_estimators=100, max_depth=5, random_state=42)
    attacker.fit(X_train, y_train)
    
    # Evaluate Attack Success
    y_pred = attacker.predict(X_test)
    acc = accuracy_score(y_test, y_pred)
    f1 = f1_score(y_test, y_pred)
    
    return acc, f1

def evaluate_privacy():
    script_dir = os.path.dirname(os.path.abspath(__file__))
    project_base = os.path.abspath(os.path.join(script_dir, "../.."))
    real_path = os.path.join(project_base, "data/processed/clinical_features_dataset.csv")
    syn_path = os.path.join(project_base, "data/synthetic/audited_synthetic_data.csv")
    eval_dir = os.path.join(project_base, "data/synthetic/evaluation")
    os.makedirs(eval_dir, exist_ok=True)
    
    print("Loading data for privacy assessment...")
    real_df, syn_df, top_genes = prepare_data(real_path, syn_path)
    
    print("Computing Exact Match Score...")
    match_count, syn_total = compute_exact_match(real_df, syn_df)
    exact_match_rate = match_count / syn_total if syn_total > 0 else 0
    
    # Fit preprocessor on combined data
    numeric_features = ['patient_info.age', 'biomarkers.tmb']
    categorical_features = ['patient_info.gender', 'biomarkers.msi']
    binary_features = [f'Mut_{g}' for g in top_genes]
    
    preprocessor = ColumnTransformer(
        transformers=[
            ('num', StandardScaler(), numeric_features),
            ('cat', OneHotEncoder(handle_unknown='ignore', sparse_output=False), categorical_features),
            ('bin', 'passthrough', binary_features)
        ])
    
    combined_df = pd.concat([real_df, syn_df], ignore_index=True)
    preprocessor.fit(combined_df)
    
    print("Computing Distance to Closest Record (DCR)...")
    X_real = preprocessor.transform(real_df)
    X_syn = preprocessor.transform(syn_df)
    dcr_metrics = compute_dcr_and_nps(X_real, X_syn)
    
    print("Simulating Membership Inference Attack (MIA)...")
    mia_acc, mia_f1 = membership_inference_attack(real_df, syn_df, preprocessor)
    
    results = {
        "Privacy_Assessment": {
            "Exact_Match": {
                "Match_Count": match_count,
                "Total_Synthetic_Records": syn_total,
                "Exact_Match_Rate": round(exact_match_rate, 4),
                "Status": "PASS" if match_count == 0 else "FAIL"
            },
            "Distance_to_Closest_Record": dcr_metrics,
            "Membership_Inference_Attack": {
                "Attacker_Accuracy": round(mia_acc, 4),
                "Attacker_F1_Score": round(mia_f1, 4),
                "Baseline_Random_Guess_Accuracy": 0.5000,
                "Status": "PASS" if mia_acc < 0.65 else "WARNING"
            }
        }
    }
    
    metrics_path = os.path.join(eval_dir, 'privacy_metrics.json')
    with open(metrics_path, 'w', encoding='utf-8') as f:
        json.dump(results, f, indent=4)
        
    print("\n--- Privacy Preservation Report ---")
    print(json.dumps(results, indent=4))
    
    if results["Privacy_Assessment"]["Exact_Match"]["Status"] == "PASS":
        print("\n[SUCCESS] No exact matches found. Patient records are not duplicated.")
    else:
        print("\n[WARNING] Exact matches found! Synthesizer is memorizing data.")
        
    if results["Privacy_Assessment"]["Distance_to_Closest_Record"]["Status"] == "PASS":
        print("[SUCCESS] DCR leakage rate is within safe limits (leakage < 5%).")
    else:
        print("[WARNING] DCR leakage rate is elevated. Check for overfitting.")
        
    if results["Privacy_Assessment"]["Membership_Inference_Attack"]["Status"] == "PASS":
        print("[SUCCESS] MIA accuracy is near random guessing (~0.5). High privacy threshold achieved!")
    else:
        print("[WARNING] MIA accuracy is high. Members are distinguishable from non-members.")

if __name__ == "__main__":
    evaluate_privacy()

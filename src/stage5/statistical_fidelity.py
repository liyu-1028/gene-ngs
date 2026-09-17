import pandas as pd
import numpy as np
from scipy.stats import ks_2samp
from sklearn.preprocessing import StandardScaler, OneHotEncoder
from sklearn.compose import ColumnTransformer
from sklearn.pipeline import Pipeline
from sklearn.manifold import TSNE
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import seaborn as sns
import json
import os
import ast

def load_data(real_path, syn_path):
    try:
        real_df = pd.read_csv(real_path, encoding='utf-8')
    except UnicodeDecodeError:
        real_df = pd.read_csv(real_path, encoding='gbk')
        
    try:
        syn_df = pd.read_csv(syn_path, encoding='utf-8')
    except UnicodeDecodeError:
        syn_df = pd.read_csv(syn_path, encoding='gbk')
        
    return real_df, syn_df

def compute_tvd(real_series, syn_series):
    # Calculate Total Variation Distance for categorical variables
    real_counts = real_series.value_counts(normalize=True).to_dict()
    syn_counts = syn_series.value_counts(normalize=True).to_dict()
    
    all_keys = set(real_counts.keys()).union(set(syn_counts.keys()))
    tvd = 0.5 * sum(abs(real_counts.get(k, 0) - syn_counts.get(k, 0)) for k in all_keys)
    return tvd

def compute_correlation_mae(real_df, syn_df, continuous_cols):
    real_corr = real_df[continuous_cols].corr(method='pearson').fillna(0).values
    syn_corr = syn_df[continuous_cols].corr(method='pearson').fillna(0).values
    mae = np.mean(np.abs(real_corr - syn_corr))
    return mae

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
        
        # Deduplicate genes per patient so we get patient-level mutation frequency
        seen_genes = set()
        for v in variants:
            if isinstance(v, dict) and 'gene' in v:
                seen_genes.add(v['gene'])
        for gene in seen_genes:
            gene_counts[gene] = gene_counts.get(gene, 0) + 1
            
    # Convert to frequencies
    return {gene: count / total_samples for gene, count in gene_counts.items()}

def evaluate_statistical_fidelity():
    script_dir = os.path.dirname(os.path.abspath(__file__))
    project_base = os.path.abspath(os.path.join(script_dir, "../.."))
    real_path = os.path.join(project_base, "data/processed/clinical_features_dataset.csv")
    syn_path = os.path.join(project_base, "data/synthetic/audited_synthetic_data.csv")
    eval_dir = os.path.join(project_base, "data/synthetic/evaluation")
    os.makedirs(eval_dir, exist_ok=True)
    
    print("Loading data for statistical fidelity assessment...")
    real_df, syn_df = load_data(real_path, syn_path)
    
    continuous_cols = ['patient_info.age', 'biomarkers.tmb']
    categorical_cols = ['patient_info.gender', 'biomarkers.msi']
    
    metrics = {
        "Univariate_Marginal_Similarity": {}, 
        "Bivariate_Correlation": {},
        "Somatic_Variant_Fidelity": {}
    }
    
    # 1. KS-Test for continuous variables
    for col in continuous_cols:
        real_data = real_df[col].dropna()
        syn_data = syn_df[col].dropna()
        if len(real_data) > 0 and len(syn_data) > 0:
            stat, p_value = ks_2samp(real_data, syn_data)
            metrics["Univariate_Marginal_Similarity"][col] = {"KS_statistic": round(stat, 4), "p_value": round(p_value, 4)}
            
    # 2. TVD for categorical variables
    for col in categorical_cols:
        real_data = real_df[col].fillna("Missing").astype(str)
        syn_data = syn_df[col].fillna("Missing").astype(str)
        tvd = compute_tvd(real_data, syn_data)
        metrics["Univariate_Marginal_Similarity"][col] = {"TVD": round(tvd, 4)}
        
    # 3. Somatic Variant Mutation Frequency TVD
    real_gene_freqs = get_gene_frequencies(real_df)
    syn_gene_freqs = get_gene_frequencies(syn_df)
    
    # Sort genes by frequency in real data
    top_genes = sorted(real_gene_freqs.keys(), key=lambda g: real_gene_freqs[g], reverse=True)[:30]
    
    gene_tvd_metrics = {}
    total_gene_diff = 0.0
    for gene in top_genes:
        f_real = real_gene_freqs.get(gene, 0.0)
        f_syn = syn_gene_freqs.get(gene, 0.0)
        diff = abs(f_real - f_syn)
        gene_tvd_metrics[gene] = {
            "Real_Frequency": round(f_real, 4),
            "Synthetic_Frequency": round(f_syn, 4),
            "Frequency_Difference_TVD": round(diff, 4)
        }
        total_gene_diff += diff
        
    avg_gene_diff = total_gene_diff / len(top_genes) if top_genes else 0.0
    metrics["Somatic_Variant_Fidelity"] = {
        "Average_Gene_Frequency_Difference": round(avg_gene_diff, 4),
        "Top_Gene_Differences": gene_tvd_metrics
    }
        
    # 4. Correlation MAE (PCD)
    mae = compute_correlation_mae(real_df, syn_df, continuous_cols)
    metrics["Bivariate_Correlation"]["Pearson_MAE"] = round(float(mae), 4)
    
    # 5. Dimensionality Reduction (t-SNE/PCA)
    print("Computing PCA and t-SNE embeddings for visualization...")
    real_df_copy = real_df.copy()
    syn_df_copy = syn_df.copy()
    real_df_copy['DataType'] = 'Real'
    syn_df_copy['DataType'] = 'Synthetic'
    
    # Downsample synthetic data for visualization to prevent performance issues
    MAX_VIZ_POINTS = 2000
    if len(syn_df_copy) > MAX_VIZ_POINTS:
        syn_df_viz = syn_df_copy.sample(n=MAX_VIZ_POINTS, random_state=42)
        print(f"Downsampled synthetic data from {len(syn_df_copy)} to {MAX_VIZ_POINTS} for visualization.")
    else:
        syn_df_viz = syn_df_copy
        
    combined = pd.concat([real_df_copy, syn_df_viz], ignore_index=True)
    
    # Extract top 15 mutated genes to include in visualization feature space
    top_viz_genes = sorted(real_gene_freqs.keys(), key=lambda g: real_gene_freqs[g], reverse=True)[:15]
    for gene in top_viz_genes:
        col_name = f'Mut_{gene}'
        combined[col_name] = 0
        for idx, row in combined.iterrows():
            val = row.get('somatic_variants', '[]')
            try:
                if pd.isna(val) or str(val) == 'nan' or not val:
                    variants = []
                elif isinstance(val, str):
                    variants = ast.literal_eval(val)
                else:
                    variants = val
            except Exception:
                variants = []
            if any(isinstance(v, dict) and v.get('gene') == gene for v in variants):
                combined.at[idx, col_name] = 1
                
    # Select features for dimensionality reduction
    binary_cols = [f'Mut_{g}' for g in top_viz_genes]
    features = combined[continuous_cols + categorical_cols + binary_cols].copy()
    for col in continuous_cols:
        features[col] = pd.to_numeric(features[col], errors='coerce').fillna(0)
    for col in categorical_cols:
        features[col] = features[col].astype(str)
        
    preprocessor = ColumnTransformer(
        transformers=[
            ('num', StandardScaler(), continuous_cols),
            ('cat', OneHotEncoder(handle_unknown='ignore', sparse_output=False), categorical_cols),
            ('bin', 'passthrough', binary_cols)
        ])
    
    X_processed = preprocessor.fit_transform(features)
    
    # PCA
    pca = PCA(n_components=2)
    X_pca = pca.fit_transform(X_processed)
    
    # t-SNE (limit perplexity to be smaller than dataset size)
    tsne = TSNE(n_components=2, perplexity=min(30, len(combined)-1), random_state=42)
    X_tsne = tsne.fit_transform(X_processed)
    
    # Plotting
    plt.figure(figsize=(14, 6))
    
    plt.subplot(1, 2, 1)
    sns.scatterplot(x=X_pca[:, 0], y=X_pca[:, 1], hue=combined['DataType'], alpha=0.6, palette={'Real':'blue', 'Synthetic':'orange'})
    plt.title('PCA: Real vs Synthetic')
    
    plt.subplot(1, 2, 2)
    sns.scatterplot(x=X_tsne[:, 0], y=X_tsne[:, 1], hue=combined['DataType'], alpha=0.6, palette={'Real':'blue', 'Synthetic':'orange'})
    plt.title('t-SNE: Real vs Synthetic')
    
    plt.tight_layout()
    plot_path = os.path.join(eval_dir, 'manifold_overlap.png')
    plt.savefig(plot_path, dpi=300)
    plt.close()
    print(f"Saved visualization to {plot_path}")
    
    # Save metrics
    metrics_path = os.path.join(eval_dir, 'statistical_fidelity_metrics.json')
    with open(metrics_path, 'w', encoding='utf-8') as f:
        json.dump(metrics, f, indent=4)
    print(f"Saved metrics report to {metrics_path}")
    
    # Print summary
    print("\n--- Statistical Fidelity Report Summary ---")
    print(f"Univariate marginal TMB: {metrics['Univariate_Marginal_Similarity'].get('biomarkers.tmb', 'N/A')}")
    print(f"Somatic Mutation TVD (Mean across top 30 genes): {metrics['Somatic_Variant_Fidelity']['Average_Gene_Frequency_Difference']}")
    print(f"Bivariate Correlation (Pearson MAE): {metrics['Bivariate_Correlation']['Pearson_MAE']}")
    
if __name__ == "__main__":
    evaluate_statistical_fidelity()

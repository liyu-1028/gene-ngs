import pandas as pd
import json
import os
import ast
import re

class SimulationOrchestrator:
    def __init__(self, data_path, output_dir, config_path="configs/mutation_coordinates.json"):
        script_dir = os.path.dirname(os.path.abspath(__file__))
        project_base = os.path.abspath(os.path.join(script_dir, "../.."))
        
        self.data_path = data_path if os.path.isabs(data_path) else os.path.join(project_base, data_path)
        self.output_dir = output_dir if os.path.isabs(output_dir) else os.path.join(project_base, output_dir)
        self.config_path = config_path if os.path.isabs(config_path) else os.path.join(project_base, config_path)
        
        os.makedirs(self.output_dir, exist_ok=True)
        self.mutation_db = self._load_config()
        
    def _load_config(self):
        with open(self.config_path, 'r') as f:
            return json.load(f)

    def load_data(self):
        print(f"Loading audited data from {self.data_path}")
        return pd.read_csv(self.data_path)

    def process_variants(self, df):
        patient_data = {}
        snv_db = self.mutation_db.get("SNV", {})
        cna_db = self.mutation_db.get("CNA", {})
        
        for idx, row in df.iterrows():
            patient_id = f"virtual_pt_{idx:04d}"
            variants_str = row['somatic_variants']
            pm, indel, cna = [], [], []
            
            try:
                variants = ast.literal_eval(variants_str)
            except Exception as e:
                print(f"Error parsing variants for {patient_id}: {e}")
                continue
                
            for v in variants:
                gene, alt = v.get('gene'), v.get('alteration')
                vaf = v.get('vaf', 0.1) or 0.1
                key = f"{gene}:{alt}"
                
                if key in snv_db:
                    chrom, pos, alt_base = snv_db[key]
                    if alt_base in ["A", "C", "G", "T"]:
                        pm.append([chrom, pos, alt_base, vaf, 1, 1])
                elif any(kw in alt for kw in ["Amplification", "扩增", "拷贝数"]):
                    if gene in cna_db:
                        chrom, start, end = cna_db[gene]
                        cna.append([chrom, start, end, 5, 1, vaf])
                elif any(kw in alt for kw in ["del", "dup", "ins"]):
                    if gene in cna_db:
                        chrom, start, _ = cna_db[gene]
                        indel.append([chrom, start + 500, 1 if "ins" in alt or "dup" in alt else 0, "A", vaf, 1, 1])
                else:
                    if gene in cna_db:
                        chrom, start, _ = cna_db[gene]
                        pm.append([chrom, start + 100, "T", vaf, 1, 1])

            patient_data[patient_id] = {'pm': pm, 'indel': indel, 'cna': cna}
        return patient_data

    def write_synggen_files(self, patient_data):
        var_base_dir = os.path.join(self.output_dir, "variant_specs")
        os.makedirs(var_base_dir, exist_ok=True)
        patient_list = []
        for pt_id, data in patient_data.items():
            pt_dir = os.path.join(var_base_dir, pt_id)
            os.makedirs(pt_dir, exist_ok=True)
            pd.DataFrame(data['pm']).to_csv(os.path.join(pt_dir, "variants.pm"), sep='\t', header=False, index=False)
            pd.DataFrame(data['indel']).to_csv(os.path.join(pt_dir, "variants.indel"), sep='\t', header=False, index=False)
            pd.DataFrame(data['cna']).to_csv(os.path.join(pt_dir, "variants.cna"), sep='\t', header=False, index=False)
            patient_list.append(pt_id)
        print(f"Generated variant specs for {len(patient_list)} patients.")
        return var_base_dir, patient_list

    def generate_bash_script(self, var_base_dir, patient_list):
        # Generate the bash script using dynamic PROJECT_BASE resolution at runtime
        # This keeps the script fully dynamic and safe for public open-source release, 
        # while resolving to absolute paths at runtime to satisfy synggen's requirements.
        script_content = f"""#!/bin/bash
# Batch simulation script with unique seeds (Dynamic Absolute Paths)
# Resolves the absolute path to the project root directory at runtime
PROJECT_BASE=$(cd "$(dirname "${{BASH_SOURCE[0]}}")/../.." && pwd)

PATIENTS=({ " ".join(patient_list) }) # Run all generated patients

for ((i=0; i<${{#PATIENTS[@]}}; i++)); do
    PT_ID="${{PATIENTS[i]}}"
    SEED=$((42 + i))
    echo "Simulating $PT_ID with seed $SEED..."
    OUT_DIR="$PROJECT_BASE/data/synthetic/reads/$PT_ID"
    mkdir -p "$OUT_DIR"
    VAR_DIR="$PROJECT_BASE/data/synthetic/variant_specs/$PT_ID"
    
    "$PROJECT_BASE/bin/synggen" mode=1 \\
        bed="$PROJECT_BASE/data/raw/reference/target_regions.bed" \\
        fasta="$PROJECT_BASE/data/raw/reference/GRCh38_full_analysis_set_plus_decoy_hla.fa" \\
        rdm="$PROJECT_BASE/data/raw/reference/models/wes_model/target_regions.rdm.gz" \\
        pbe="$PROJECT_BASE/data/raw/reference/models/wes_model/target_regions.pbe.gz" \\
        qm="$PROJECT_BASE/data/raw/reference/models/wes_model/target_regions.qm.gz" \\
        pm="$VAR_DIR/variants.pm" \\
        indel="$VAR_DIR/variants.indel" \\
        cna="$VAR_DIR/variants.cna" \\
        threads=4 \\
        nreads=200000 \\
        tc=0.3 \\
        out="$OUT_DIR" \\
        seed=$SEED
done
echo "Batch simulation completed."
"""
        script_path = os.path.join(self.output_dir, "run_simulation.sh")
        with open(script_path, 'w', newline='\n') as f:
            f.write(script_content)
        print(f"Generated script: {script_path}")

if __name__ == "__main__":
    orchestrator = SimulationOrchestrator(data_path="data/synthetic/audited_synthetic_data.csv", output_dir="data/synthetic")
    df = orchestrator.load_data()
    patient_data = orchestrator.process_variants(df)
    var_base_dir, pt_list = orchestrator.write_synggen_files(patient_data)
    orchestrator.generate_bash_script(var_base_dir, pt_list)

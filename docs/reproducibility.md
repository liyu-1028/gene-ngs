# Reproducibility Guide

This document describes **exactly** how the archived synthetic cohorts were produced
and how to reproduce the pipeline end-to-end. Every command below reflects the
code as committed in this repository.

---

## 1. System requirements

| Item | Requirement |
|---|---|
| OS | Ubuntu 22.04 LTS (any modern Linux works) |
| Python | 3.10.x |
| RAM | ≥ 16 GB (CTGAN training on 10k rows is the peak) |
| Disk | ~10 GB (code + models) + ~40 GB if you regenerate the full 2,000-patient FASTQ corpus |
| Compiler | gcc / make (to build synggen) |
| External tool | `samtools` (only for anchor-coverage validation utilities) |

## 2. Python environment (verified versions)

The pipeline was run with the exact package versions below
(`pip install -r requirements.txt` installs these pins):

```
pandas==2.3.3          numpy==1.26.4        scipy==1.15.3
scikit-learn==1.7.2    matplotlib==3.10.9   seaborn==0.13.2
shap==0.49.1           sdv==1.11.0          torch==2.5.1
requests==2.34.2       python-dotenv==1.2.2 openai==2.38.0
pillow==12.2.0
```

```bash
python3.10 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

> `torch` is required by SDV/CTGAN; CPU-only wheels are sufficient
> (`pip install torch --index-url https://download.pytorch.org/whl/cpu`).

## 3. Build synggen (read simulator)

```bash
git clone https://bitbucket.org/CibioBCG/synggen.git
cd synggen/src && make            # produces the `synggen` binary
mkdir -p /path/to/project/bin && cp synggen /path/to/project/bin/
```

Notes learned during production:

* **Always invoke the binary via an absolute path** — relative invocation can
  segfault (the binary resolves resources relative to its own location).
* The pipeline drives synggen in `mode=1` (matched tumor-normal style, BED-guided).

## 4. Reference data

| File | Source |
|---|---|
| `GRCh38_full_analysis_set_plus_decoy_hla.fa` | UCSC Genome Browser (GCA_000001405.15 + decoy + HLA) |
| `target_regions.bed` | WES target manifest (e.g., IDT/Illumina exome BED, sorted, merged) |
| `target_regions.{rdm,pbe,qm}.gz` | synggen model files — generate with synggen's bundled
  model tooling from cohort-level mutation statistics (see synggen README), or use the
  WES example model shipped in the synggen repository |

Index the FASTA once: `samtools faidx GRCh38_....fa`

## 5. API key (Stage 3 LLM audit only)

Stage 3 combines deterministic hard rules with an optional LLM plausibility audit.
The LLM path needs a DeepSeek-compatible key in a `.env` file at the project root:

```
DEEPSEEK_API_KEY=sk-xxxxxxxx
```

Without a key, run the rules-only audit path (the `audit_stage` column will record
`stage3_hardrule` instead of `stage3_llm`).

Optional environment knobs:

| Variable | Default | Purpose |
|---|---|---|
| `AUDIT_SEED` | — | Fixed seed for LLM-audit sampling order |
| `LLM_AUDIT_N` | all | Cap the number of LLM-audited records |
| `SHAP_SAMPLE` | 0 (full) | Subsample size for SHAP in Stage 5 TSTR utility |

## 6. Stage-by-stage commands

Inputs/outputs are the files used for the archived cohorts; adapt paths to your layout.

### Stage 1 — data governance
Ingest the real cohort (cBioPortal DataHub `msk_impact_2017.tar.gz`, ODbL) into the
pipeline schema. See `src/stage1/` and the converter `src/stage0/convert_msk_impact_to_pipeline.py`.

### Stage 2 — tabular generation (CTGAN)
```bash
python src/stage2/run_generation.py    # trains CTGANSynthesizer on the real-cohort
                                       # dataframe, emits synthetic tables
```
Model class: `sdv.single_table.CTGANSynthesizer` (see `src/stage2/ctgan_generator.py`).

### Stage 3 — clinical audit
```bash
python src/stage3/clinical_auditor.py  # hard rules (+ LLM if key present)
```
Produces the 9-column audited CSV (incl. `audit_stage`).

### Stage 4 — read simulation
```bash
python src/stage4/simulation_orchestrator.py
```
Per patient *i* the orchestrator writes variant specs (`variants.pm/.indel/.cna`)
and calls:
```
bin/synggen mode=1 bed=... fasta=... rdm=... pbe=... qm=... \
  pm=variants.pm indel=variants.indel cna=variants.cna \
  threads=4 nreads=200000 tc=0.3 out=out/virtual_pt_XXXX seed=$((42+i))
```
Determinism: `seed = 42 + patient_index` for every patient.
If the run is interrupted, delete the patient's incomplete output directory and
re-run — the orchestrator skips completed patients.

After simulation, merge per-thread FASTQs, gzip, and verify:

```bash
# per patient: expect exactly 800,000 lines (= 200,000 reads)
zcat out/virtual_pt_0000.fastq.gz | wc -l
```

### Stage 5 — technical validation
```bash
python src/stage5/statistical_fidelity.py REAL.csv SYN.csv   # KS / TVD
python src/stage5/utility_tstr.py REAL.csv SYN.csv           # TRTR vs TSTR AUC
python src/stage5/privacy_assessment.py REAL.csv SYN.csv     # exact match / DCR / MIA
```
Anchor-level VAF recovery (hi-coverage physical validation) is orchestrated by
`src/stage0/run_hicov_minisim_validation.py` (uses `samtools faidx`).

## 7. Expected results

With the archived inputs the validation suite reproduces (within Monte-Carlo
tolerance for the generative stage):

* TMB KS = 0.115, gender TVD = 0.113 (2017 cohort)
* TRTR AUC 0.791 / TSTR AUC 0.425 (2017), TRTR 0.751 / TSTR 0.565 (50K cohort)
* Privacy: 0 exact record matches, DCR leakage 0.0, MIA AUC ≈ 0.68
* Hi-coverage anchor study: 102/124 anchors detected (ALT≠REF), Pearson r = 0.851
  between specified and observed VAF

### Privacy metric control experiment

The "Membership_Inference_Attack" entry in `privacy_metrics.json` (attacker accuracy
0.679 > 0.5 baseline) is a real-vs-synthetic discriminability proxy, not a
memorization test. `src/stage5/mia_control_experiment.py` establishes its background:

| Experiment | Accuracy | Meaning |
|---|---|---|
| R0 real vs synthetic | 0.684 | reproduces published 0.679 |
| C1 realA vs realB (disjoint REAL halves) | 0.496 | intrinsic background ≈ chance |
| C2 synA vs synB | 0.503 | attacker calibrated |
| C3 C1 × 3 seeds | 0.490–0.511 | background stable |

Since perfect memorization would *reduce* discriminability toward 0.5, a value of
0.68 quantifies the generator fidelity gap (consistent with the disclosed TRTR/TSTR
utility gap), while membership privacy is covered by Exact_Match = 0 and
DCR leakage = 0.0 (both PASS). Results: `evaluation/mia_control_results.json`.

## 8. Known quirks

* synggen emits one FASTQ per thread (`target_regions_N.R1.fastq`) — concatenate
  thread files, then gzip, then validate the 200,000-read count.
* synggen model files (`rdm`) must be regenerated if you change the target BED;
  row order/accumulation follows synggen's internal expectations (see its docs).
* CPU-only torch is sufficient; GPU does not change outputs but speeds up CTGAN.

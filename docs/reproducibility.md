# V2 reproducibility and verification guide

## Version and access status

The dataset record is [ScienceDB DOI 10.57760/sciencedb.0131a](https://doi.org/10.57760/sciencedb.0131a). Version 1 is public. Version 2 (V2) has been submitted for ScienceDB review and is the package described by the revised manuscript. V2 should be verified after it is approved and publicly downloadable; a local or reviewer-visible package is not a substitute for that final external check.

The V2 package retains the V1 CSV and FASTQ tar. It adds `sample_index.tsv`, `fastq_patient_mapping_v2.csv`, classified spike-in truth and specification files, quality flags, corrected analyses, V1 archive-audit outputs, independent public-read recount outputs, figure data and `MANIFEST_V2.md5`.

## Source bundle

`gene-ngs-v2.0.0-source.tar.gz` contains the source scripts used to prepare the V2 data package and manuscript-derived outputs. Its contents are listed in [V2_CODE_RELEASE.md](../V2_CODE_RELEASE.md). It does not include source clinical records, the 36.46-GB FASTQ archive, API credentials or any unavailable historical generator asset.

Install only the packages needed for the bounded V2 checks:

```bash
python3 -m pip install -r requirements-v2-quick.txt
```

## Verification workflow

1. Download the approved V2 package and retain its `MANIFEST_V2.md5`.
2. Verify all downloaded small objects against that manifest; verify the retained CSV and FASTQ tar against the identities recorded in `V2_RETAINED_CORE_ASSETS.tsv`.
3. Run `scripts/audit_fastq_archive.py` against the full retained FASTQ tar and V2 mapping. This checks archive identity, all 2,000 member paths/MD5s, gzip integrity, FASTQ structure, read counts and read lengths.
4. Run `scripts/verify_fastq_anchor_signatures.py` using the V2 independent hg38 contexts and mapping, then use `scripts/prepare_public_read_validation_figure_data.py` to produce the two public-read figure panels.
5. Run `scripts/check_v2_release_candidate.py` using the archive and recount reports. It checks the index, mapping, manifest and release structure.

The 21-mer recount evaluates a preselected diagnostic subset. It demonstrates that the retained file-level mapping and saved counts agree for 183 distinguishable REF/ALT anchors; it does not establish that every listed clinical-table variant was realized as an observed read-level event.

## Interpretation limits

The `fastq_spikein_truth.tsv` table records inputs and transformations into simulator specifications, not an observed call set. Only `exact_snv` has an alteration-specific coordinate. `gene_fallback_pm`, `heuristic_indel`, conflict, simulator-error and independent-hg38 REF=ALT statuses must be excluded from unambiguous variant-positive scoring.

The exact historical reference FASTA, target BED, WES error models, simulator binary and high-depth mini-simulation reads are unavailable. V2 therefore supports controlled mapping, integrity and bounded read-signature checks. It does not make a byte-for-byte historical regeneration claim or provide an all-variant benchmark.

The corrected V2 analyses identify one normalized full-record real--synthetic match. DCR/NNDR are descriptive diagnostics, and the historical real-versus-synthetic discriminator is not a membership-inference attack. The V2 documentation records the corrected values and their scope.

# V2 code release

## Scope

`gene-ngs-v2.0.0-source.tar.gz` is the source bundle synchronized with the submitted ScienceDB V2 package. It contains scripts for V2 metadata/index creation, corrected analyses, manuscript tables and figures, FASTQ archive auditing, selected public-read recounts and release checks. It is source code only; data are available from the ScienceDB record.

- Version: `v2.0.0-submitted`
- ScienceDB record: <https://doi.org/10.57760/sciencedb.0131a>
- V2 public status: submitted for ScienceDB review; not yet approved
- Archive SHA-256: `c44c69fdfb7d6f0913002159180b93fadc77ea6162a607bb99c571290e41b550`

## Archive layout

```text
gene-ngs-v2.0.0-source/
├── README.md
├── requirements-v2-quick.txt
├── scripts/
│   ├── prepare_sciencedb_v2_metadata.py
│   ├── prepare_sciencedb_v2_release_indexes.py
│   ├── reanalyse_2017_quick.py
│   ├── prepare_stage5_figure_data.py
│   ├── generate_manuscript_results_tables.py
│   ├── render_stage5_quantitative_figures.py
│   ├── render_manuscript_overview.py
│   ├── audit_fastq_archive.py
│   ├── verify_fastq_anchor_signatures.py
│   ├── correct_delivery_anchor_validation.py
│   ├── prepare_public_read_validation_figure_data.py
│   ├── check_v2_release_candidate.py
│   ├── check_numbers.py
│   ├── check_manuscript.py
│   ├── build_sciencedb_v2_upload_package.py
│   └── package_v2_verification_code.py
└── v2_verification_code/
    ├── README.md
    ├── requirements.txt
    └── bounded audit and recount scripts
```

## Release boundary

A user must download the approved V2 package and compare it with `MANIFEST_V2.md5` before describing this code release as a verified public-data release. The V2 source bundle does not contain raw source clinical records, the original simulation environment, credentials or licence-restricted material.

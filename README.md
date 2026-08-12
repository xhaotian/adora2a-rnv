# ADORA2A transcript evidence in retinal neovascularization

This repository is the version-locked reproducibility release for the PLOS ONE submission. It regenerates the deposited-sample/context-level synthesis of 16 mouse cohorts (14/16 positive), tissue-compartment sensitivity, exploratory Pustejovsky–Rodgers funnel-asymmetry diagnostic, donor-level human sensitivity models, non-pooled human context tables, Source Data, PRISMA checklist, four main figures, and five supplementary figures.

The scientific boundary is narrow: the mouse synthesis estimates the average standardized ADORA2A direction across eligible P17 OIR retinal transcriptomic contexts. The historical-list-based six-gene human score is sensitivity-only. Human contexts remain separate and are not pooled. Outputs do not establish receptor activation, causality, or therapeutic efficacy.

## Reproduce

Use Python 3.13 and R 4.5 with the versions in `requirements.lock.txt` and `R-packages.lock.txt`, then run:

```bash
python -m pip install -r requirements.lock.txt
Rscript -e 'install.packages("metafor", repos="https://cloud.r-project.org")'
bash run_all.sh
```

Generated files are written under `review_remediation/` and `final_submission/`. The run ends with numerical, boundary, figure, and checksum validation.

## Public-data retrieval

No new sequencing data were generated. GEO datasets can be obtained from `https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=<ACCESSION>`. Non-GEO mouse cohorts are NCBI BioProjects PRJNA1227125 and PRJNA1347474; SRA run identifiers and source URLs are listed in `review_remediation/03_systematic_search/non_geo_fastq_manifest.tsv`. The complete item-level search audit is `review_remediation/03_systematic_search/Search_All_Sources.tsv`.

The release-distributed tables are audited public-data extracts sufficient for deterministic reproduction of all submitted numerical results and displays. Full raw-data re-quantification instructions and accession checksums are in `docs/DATA_RETRIEVAL.md`.

## License

Code is released under the MIT License. Public source datasets remain governed by their originating repositories and publications.

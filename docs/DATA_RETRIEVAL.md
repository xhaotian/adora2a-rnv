# Public-data retrieval

## Mouse GEO cohorts

The 14 included GEO accessions are listed in `MOUSE_COMPARTMENT_REGISTRY.tsv`. Download processed matrices and SOFT metadata from each accession record. `Search_All_Sources.tsv` preserves the source URL and decision for every retrieved object.

## Mouse non-GEO cohorts

- PRJNA1227125: use the SRA run identifiers in `non_geo_fastq_manifest.tsv`.
- PRJNA1347474: use the listed SRA runs; PRJNA483866 is its duplicate-cohort counterpart and is not counted separately.

The original workflow used `fasterq-dump`, GENCODE mouse vM36 transcripts, Salmon 1.10.3, transcript-to-gene summation for *Adora2a*, and edgeR log2 CPM with prior count 0.5. Retrieval and verification fields are retained in the manifests.

## Human donor reconstruction

Official count matrices and metadata are available from GSE165784 and GSE245561. Counts are summed within donor for quality-controlled endothelial cells, converted to log2(CPM+0.5), and used by the included donor-model script. Context-specific accessions and biological-unit boundaries are listed in `CONTEXTUAL_HUMAN_REGISTRY.tsv`.

## Voigt 2026 MNV input

The analysis used bioRxiv version 1, supplementary file `media-2.xlsx`, downloaded 22 July 2026, SHA-256 `3549a93926b4ea689fad4bfec54e8981f59c129d4f608ce18eda86fdc487d529`. Version-specific source: https://www.biorxiv.org/content/10.64898/2026.03.30.714946v1.

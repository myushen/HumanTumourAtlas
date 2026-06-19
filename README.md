# HumanTumourAtlas

An R-based data engineering pipeline that ingests raw single-cell RNA-seq count
data from the [Human Tumour Atlas Network (HTAN)](https://humantumoratlas.org)
and produces harmonised, analysis-ready artefacts consumable by
[cellNexus](https://github.com/MangiolaLaboratory/cellNexus) and
[HPCell](https://github.com/MangiolaLaboratory/HPCell).

---

## Background

HTAN is a multi-centre NCI-funded initiative that profiled tumour ecosystems
across cancer types and disease stages using single-cell and spatial genomics.
Data are deposited on [Synapse](https://www.synapse.org) per centre and are
publicly accessible via the
[HTAN Data Portal](https://data.humantumoratlas.org).

### The core challenge

Each contributing centre chose its own file organisation, naming conventions,
and processing strategy independently.  There is **no universal parse
strategy** across centres:

| Centre | Atlas ID | Format | L4 strategy |
|--------|----------|--------|-------------|
| HTAN HTAPP | HTA1 | L3 MTX triplet | L4 TSV batch file with global cell index |
| HTAN BU | HTA3 | L3 dense CSV | L4 CSV cell filter |
| HTAN OHSU | HTA9 | L3 10X HDF5 | No L4 — all barcodes used |
| HTAN WUSTL | HTA12 | L4 per-sample Seurat RDS | RNA assay extracted; SCT is active assay |
| HTAN Duke | HTA6 | L4 merged Seurat RDS (22 samples) | L3 MTX fallback for 8 L4-absent samples |
| HTAN DFCI | HTA5 | L3 pooled HDF5 | Pooled lanes — demultiplexing pending |
| HTAN MSK | HTA8 | L3 MTX / H5 | Per-sample L4 filtered objects |

The pipeline therefore applies a **centre-by-centre, organ-by-organ**
classification strategy: each centre has a dedicated parser registered in
`dev/0_center_parser_registry.R`, and a unified cell index
(`*_cell_index_df`) is built and stored as a persistent targets object for
every centre to ensure reproducible, auditable barcode → integer cell ID
mapping.

An additional challenge is scale: some samples exceed **1 million cells**,
read from raw 10X barcodes / features / matrix triplets.  The pipeline
addresses this via lazy chunked loading, SLURM-backed
[crew](https://wlandau.github.io/crew/) elastic worker pools, and
[targets](https://docs.ropensci.org/targets/) dependency tracking so only
invalidated samples are ever re-processed.

---

## Pipeline steps

```
1  Download source data       dev/1_download_synapse_counts_target_script.R
2  Parse counts → h5ad        dev/2_HTAN_parse_all_centers_targets.R
3  Harmonise metadata         dev/4_HTA_generate_cell_and_sample_donor_metadata_v0.2.0.R
4  Process HPCell             dev/4_execute_hpcell_for_hta_samples.R
5  Append QC annotations      dev/4_HTA_generate_cell_and_sample_donor_metadata_v0.2.0.R
6  Compute count caches       dev/5_prepare_pseudobulk_local_counts_cache.R
```

### Step 1 — Retrieve metadata from HTAN

Metadata (files, samples, donors, biospecimens) is exported from the
[HTAN Data Portal](https://data.humantumoratlas.org/explore) and stored under
`inst/extdata/`.  The portal supports filtering by organ type (Lung, Breast,
…), assay (scRNA-seq), level, and download source (Synapse) before export.

### Step 2 — Download source data from Synapse

Raw count files are downloaded programmatically via
[synapser](https://r-docs.synapse.org/articles/synapser.html).  A Synapse
personal access token is required:

```r
# Add to ~/.Renviron
SYNAPSE_TOKEN="your_token_here"
```

Downloads are orchestrated with a `targets` pipeline backed by a SLURM crew
controller so multiple files are fetched in parallel and already-downloaded
files are skipped automatically.

### Step 3 — Parse counts to h5ad (`dev/2_HTAN_parse_all_centers_targets.R`)

The central parsing step.  For each centre:

1. A `*_cell_index_df` target is built (persistent, queryable with
   `tar_read()`) that maps every barcode to a deterministic integer
   `cell_index`.  For L4-backed centres indices are derived from sorted
   barcodes in the L4 file; for L3-only samples indices are assigned
   sequentially over the barcodes file.
2. Each sample is parsed by the centre-specific function in
   `dev/0_center_parser_registry.R` (e.g. `parse_HTAPP`, `parse_OHSU`,
   `parse_Duke_seurat`).  Only the `counts` assay and `colData` are retained;
   reduced dims and auxiliary slots are stripped to keep objects lean.
3. Gene IDs are normalised to Ensembl identifiers via `convert_gene_to_ensembl`.
4. Each sample is written as a gzip-compressed `.h5ad` file.

Workers are dispatched across an elastic SLURM pool (5 GB → 160 GB with
automatic fallback) so memory-hungry samples are handled without manual
intervention.

### Step 4 — Harmonise metadata (`dev/4_HTA_generate_cell_and_sample_donor_metadata_v0.2.0.R`)

Reads all `.h5ad` `colData` sections (memory-efficient; counts not loaded) and
joins them with the HTAN files / samples / donors metadata to produce
standardised cell-level and sample-level metadata tables aligned with the
[cellNexus](https://github.com/MangiolaLaboratory/cellNexus) schema:

- `cell_id`, `sample_id`, `file_id_cellNexus_single_cell`,
  `file_id_cellNexus_pseudobulk`
- Donor fields: `donor_id`, `sex`, `diagnosis_age`, `ethnicity`
- Disease fields: `tissue`, `disease`, `primary_diagnosis`
- Preserved cell type columns (`cell_type_*`) where present in source data

Output: `hta_cell_metadata_v*.parquet` + `hta_sample_metadata_v*.parquet`
(ZSTD-compressed).

### Step 5 — Process HPCell (`dev/4_execute_hpcell_for_hta_samples.R`)

Runs [HPCell](https://github.com/MangiolaLaboratory/HPCell) on each sample via the
harmonised metadata.  HPCell performs empty droplet, doublet, dead cells detection,
and cell-type harmonised annotation (`cell_type_unified_ensemble`) in a SLURM-parallelised targets pipeline.

### Step 6 — Append QC annotations

HPCell QC flags and cell-type calls are joined back into the harmonised
metadata, producing the *ready* version of `hta_cell_metadata` with
`is_cell`, `cell_type_harmonised`, and per-sample QC summaries.

### Step 7 — Compute count caches (`dev/5_prepare_pseudobulk_local_counts_cache.R`)

For each sample, a targets pipeline computes and caches:

| Representation | Description |
|---------------|-------------|
| Raw counts | Integer UMI counts per gene × cell |
| CPM | Counts per million (library-size normalised) |
| Rank | Per-cell rank transformation |
| SCT | Regularised negative binomial normalisation (Seurat SCTransform) |
| Pseudobulk | Summed raw counts per sample × cell_type (`cell_type_unified_ensemble`) for bulk-style DE |

---

## Repository layout

```
dev/
  0_center_parser_registry.R          Parser functions for every centre
  1_download_synapse_counts_*.R        Metadata retrieval + Synapse download
  2_HTAN_parse_all_centers_targets.R   Main targets pipeline (all centres)
  2_HTAN_parse_*_targets.R             Centre-specific legacy pipelines
  3_HTA_generate_*_metadata*.R         Metadata generation (v0.1)
  4_HTA_generate_*_metadata_v0.2.0.R   Metadata harmonisation (v0.2)
  4_execute_hpcell_for_hta_samples.R   HPCell execution
  5_prepare_pseudobulk_*.R             Count cache generation
inst/extdata/
  files_metadata_*.tsv                 HTAN file-level metadata exports
  samples_metadata_*.tsv               Sample metadata
  donors_metadata_*.tsv                Donor/patient metadata
  {Centre}/                            Centre-specific supplementary files
```

---

## Infrastructure

| Tool | Role |
|------|------|
| [targets](https://docs.ropensci.org/targets/) | Reproducible pipeline DAG, smart re-execution |
| [tarchetypes](https://docs.ropensci.org/tarchetypes/) | `tar_group` / branching utilities |
| [crew](https://wlandau.github.io/crew/) + [crew.cluster](https://wlandau.github.io/crew.cluster/) | Elastic SLURM worker pool (5–160 GB) |
| [synapser](https://r-docs.synapse.org) | Synapse authentication and file download |
| [Seurat](https://satijalab.org/seurat/) / [SingleCellExperiment](https://bioconductor.org/packages/SingleCellExperiment/) | Count matrix handling |
| [zellkonverter](https://theislab.github.io/zellkonverter/) | SCE ↔ AnnData (h5ad) conversion |
| [DropletUtils](https://bioconductor.org/packages/DropletUtils/) | HDF5 10X file reading |
| [arrow](https://arrow.apache.org/docs/r/) | Parquet I/O for metadata tables |
| [HPCell](https://github.com/MangiolaLaboratory/HPCell) | QC, doublet detection, cell-type annotation |
| [cellNexus](https://github.com/MangiolaLaboratory/cellNexus) | Downstream single-cell data access layer |

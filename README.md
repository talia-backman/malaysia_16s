# Malaysia Marine 16S

Analysis of PacBio HiFi 16S rRNA sequencing data from Malaysian coastal marine habitats (coral, mangrove, and seagrass).

## Dataset

Final dataset:

* 107 biological samples
* 1,806 ASVs
* Environmental data from Bio-ORACLE
* Sample metadata in `MetaData_eDNA_Malaysia_withQC.tsv`

## Repository structure

* `data/` – external datasets and Bio-ORACLE data
* `scripts/` – bash scripts for computationally intensive processing steps
* `scripts/R/` – R processing and analysis scripts
* `output/` – intermediate and statistical outputs
* `figures/` – manuscript figures

## Analysis workflow

R scripts are located in `scripts/R/` and are numbered in approximate workflow order:

```text
01_trim_filter_reads.R
02_dada2_asv_inference.R
03_build_phyloseq_objects.R
04_clean_phyloseq.R
05_build_phylogeny.R
06_alpha_diversity.R
07_beta_diversity.R
08_spatial_connectivity.R
09_important_taxa.R
10_gdm_connectivity.R
11_environmental_analysis.R
12_habitat_vs_geography.R
13_nestedness_turnover.R
14_core_microbiome_comparison.R
15_habitat_core_overlap.R
```

Shared R functions are stored in `scripts/R/functions.R`.

The first two computationally intensive processing steps can be run on a SLURM-based computing cluster using the corresponding bash scripts in `scripts/`.

## Figure outputs

* Figure 1: Beta diversity
* Figure 2: Spatial connectivity
* Figure 3: Habitat- and coast-associated taxa
* Figure 4: GDM and connectivity analyses
* Figure 5: Environmental analyses
* Figure S1: Sequencing characteristics and taxonomic composition
* Figure S2: Alpha diversity
* Figure S3: Environmental PCA
* Figure S4: Habitat versus geographic effects
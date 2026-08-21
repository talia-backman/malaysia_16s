# Malaysia Marine 16S

Analysis of PacBio HiFi full-length 16S rRNA sequencing data from coral reef, mangrove, and seagrass-associated seawater communities across the east and west coasts of Peninsular Malaysia.

## Dataset

Final analytical dataset:

- 107 environmental samples
- 1,806 ASVs
- Coral reef, mangrove, and seagrass habitats
- East and west coasts of Peninsular Malaysia
- Environmental data from Bio-ORACLE v2.2
- Sample metadata in `MetaData_eDNA_Malaysia_withQC.tsv`

Raw sequence data are available through NCBI BioProject `PRJNA1516003`.

## Repository structure

- `data/` – external datasets and Bio-ORACLE data
- `scripts/` – bash scripts for computationally intensive processing steps
- `scripts/R/` – R processing and analysis scripts
- `output/` – intermediate data and statistical outputs
- `figures/` – main and supplemental manuscript figures
- `manuscript/` – manuscript-related files

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

## Manuscript figures

| Figure | Description | Script |
|---|---|---|
| Figure 1 | Beta diversity and habitat/coast structure | `07_beta_diversity.R` |
| Figure 2 | Spatial connectivity and distance-decay | `08_spatial_connectivity.R` |
| Figure 3 | Habitat- and coast-associated bacterial taxa | `09_important_taxa.R` |
| Figure 4 | GDM and ocean-connectivity analyses | `10_gdm_connectivity.R` |
| Figure 5 | Environmental gradients and microbial community structure | `11_environmental_analysis.R` |
| Figure S1 | Sequencing characteristics and taxonomic composition | `04_clean_phyloseq.R` |
| Figure S2 | Alpha diversity | `06_alpha_diversity.R` |
| Figure S3 | Bio-ORACLE environmental PCA | `11_environmental_analysis.R` |
| Figure S4 | Habitat versus coast effects on community dissimilarity | `12_habitat_vs_geography.R` |

## Manuscript tables

| Table | Description | Source |
|---|---|---|
| Table 1 | Sampling distribution across coast, habitat, and protected-area status | Sample metadata |
| Table S1 | Sample metadata, sequencing depth, and QC status | `04_clean_phyloseq.R` |
| Table S2 | Alpha diversity mixed-model results | `06_alpha_diversity.R` |
| Table S3 | PERMANOVA, geographic envfit, and dispersion statistics | `07_beta_diversity.R` |
| Table S4A | Sample-level habitat versus coast comparisons | `12_habitat_vs_geography.R` |
| Table S4B | Site-level habitat versus coast robustness analysis | `12_habitat_vs_geography.R` |
| Table S4C | Habitat centroid distances and seagrass intermediacy | `07_beta_diversity.R` |
| Table S4D | Sørensen dissimilarity, turnover, and nestedness among habitat pairs | `13_nestedness_turnover.R` |
| Table S4E | Habitat core ASV overlap | `15_habitat_core_overlap.R` |
| Table S5 | Mantel distance-decay tests | `08_spatial_connectivity.R` |
| Table S6 | GDM and MRM model results | `10_gdm_connectivity.R` |
| Table S7 | Bio-ORACLE environmental envfit statistics | `11_environmental_analysis.R` |
| Table S8 | Habitat- and coast-associated bacterial taxa | `09_important_taxa.R` |
| Table S9 | Overlap with the published mangrove core microbiome | `14_core_microbiome_comparison.R` |

## Reproducibility

Analyses were conducted in R 4.4.3. Major packages include DADA2, phyloseq, vegan, corncob, gdm, ecodist, betapart, and sdmpredictors.

The repository version associated with the manuscript is permanently archived on Zenodo:

**DOI:** https://doi.org/10.5281/zenodo.22032741
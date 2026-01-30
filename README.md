# targeted-sci-Plex
# Analysis and Plotting: TRTL & Targeted sci-Plex

This directory contains the necessary code for running analysis and plotting figures for the TRTL & targeted-sci-Plex paper.

## Directory Structure

* **processFromRaw**: Contains the initial processing scripts used to handle raw sequencing data and generate preliminary count matrices.
* **annotationScripts**: Contains scripts to convert unfiltered data to filtered, annotated data as well as some minimal analysis. The scripts are intended to be run in the order of the number shown.
* **plotScripts**: Contains scripts that do the bulk of the plotting and analysis of the processed data.
* **helperScripts**: Contains a tools script that has helper and plotting functions.
* **graphs**: Contains sub-figures used to construct figures in the main text and supplement.


---

## Data Availability

The following datasets are available via the Gene Expression Omnibus (GEO).

| Accession | Title |
| :--- | :--- |
| **GSE317074** | A modular transcript enrichment strategy for scalable, atlas-aligned, and clonotype-resolved single-cell transcriptomics |
| **GSM9464724** | sci-Plex-mouse-brain |
| **GSM9464725** | sci-RNA-TRTL-mouse-1 |
| **GSM9464726** | sci-RNA-TRTL-mouse-2 |
| **GSM9464727** | targeted-sci-Plex-Jurkat |

---

## Analysis Pipeline

The workflow follows a sequential progression through the directories:

1. **Initial Processing**: The pipeline begins with `processFromRaw` to go from fastq files to cds objects (containing cell, gene, and counts information), TCR counts, and hash counts.
2. **Data Refinement**: The scripts in `annotationScripts` perform filtering and cell-type annotation, resulting in the final processed objects.
3. **Visualization and Analysis**: The `plotScripts` utilize the annotated objects to generate the final analytical results and the figures stored in `graphs`.
4. **Shared Utilities**: The `helperScripts` provide the necessary functions and used in multiple scripts

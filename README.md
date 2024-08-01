# Impact of sampling strategies on taxonomic classification using codon metrics : applications to Ascomycetes

## Overview

In this study we compiled a dataset of coding sequences from 197 species within the Ascomycota phylum. We computed 26 characteristics and physico-chemical properties of these sequences. To address the imbalance in species representation, we evaluated the impact of different sampling techniques on multi-level classification using classical machine learning methods. We aimed to assess how these techniques influence the performance of classification models when faced with imbalanced data.

## CDS extraction and metrics computation

This pipeline is designed to process genomic annotations, extract coding sequences (CDS), filter CDS based on specific criteria, compute various codon usage and sequence-based measures, and generate final datasets for taxonomic classification studies. All the scripts are in **CDS_extract_metrics** repository.

#### Requirements

- **Operating System:** MacOS X Yosemite or later
- **Python:** Version 3.8.2 or later
- **Python Modules:**
  - argparse
  - BioPython
  - math
  - numpy
  - optparse
  - os
  - re
  - subprocess
  - sys
- **Other Tools:**
  - gffread v0.11.8
  - Perl v5.18.2
  - codonW version 1.3

#### Installation

Ensure all dependencies are installed:

```bash
pip install -r requirements.txt
```
#### Data Preparation

The pipeline supports genome annotations from the following sources:

- RIKEN
- JGI
- iGenolevures
- NCBI

#### Steps to Prepare Data:

##### Genome Annotation and CDS Extraction:
- Convert GFF3, GFF2, or EMBL files to CDS FASTA files using provided scripts.
- Ensure CDS names are unique and verify the integrity of extracted sequences.

##### CDS Filtration:
Apply filtration criteria to retain only valid CDS, ensuring they meet length, start codon, stop codon, and nucleotide composition requirements.

##### Codon Counting:
Count codons for both normal and wobble codon positions.

##### Sequence-Based Measures:
Calculate GC content, nucleotide counts, and other sequence-based metrics at the genome and CDS levels.

##### Codon Usage and Optimal Codons:
Determine optimal codons using codonW and calculate related measures such as CAI, CBI, and Fop.

##### Codon Context Measures:
Analyze codon context and compute related metrics like Fpc, Fav, Boc, Ipc, Iav, and Bic.

##### Protein Composition:
Calculate physico-chemical properties of proteins translated from CDS, including Gravy, Aromo, pI, and II.

#### Pipeline Steps

The pipeline consists of the following steps:

##### Annotation Conversion:
Use gffread or custom scripts to convert annotation files to FASTA format.

##### CDS Extraction and Filtration:
Run checkCDS.py to filter CDS based on defined criteria.

##### Codon Counting and Usage Analysis:
Run Codon_count_V3.py and codonW to compute codon counts and usage metrics.

##### Context Census and Sequence Measures:
Execute contextCensus.py to analyze codon contexts and compute related metrics.

##### Final Measures Aggregation:
Aggregate all computed measures into final tables for further analysis.

#### Output Files

The pipeline generates the following output files:

- Filtered CDS Sequences: One FASTA file per species.
- Codon Counts: Two files per species (complete and wobble).
- CodonW Analyses: One summary file per species.
- Context Counts: Two files per species (complete and wobble).
- Measures: Three files per species (complete, wobble, aminoacid).
- Taxonomy File: One file (genome-species.csv).

For more detailed information on using the pipeline (command line examples), please refer to the procedure_CDS.md

## Codon Metrics Classification and Analysis Pipeline

The **sampling_classification** repository contains a set of R scripts designed to perform codon metrics classification, evaluate model performance, and generate figures for a research paper. The pipeline includes data preparation, classification, and visualization steps, as well as utilities for working with graphics and essential functions for analysis.

### Overview

The pipeline is structured into several key components:

- Data Preparation: Splitting datasets into training and test sets using various sampling strategies.
- Classification: Applying multiple classification methods to the prepared datasets and evaluating their performance.
- Visualization: Generating plots to visualize the results of the classification tasks.
- Utilities: Helper functions for data handling, classification, and visualization.

### Requirements

- **R: Version 4.1 or later**
- **R Packages:**
    - ggplot2
    - caret
    - glmnet
    - data.table
    - e1071
    - class
    - randomForest
    - themis
    - groupdata2
    - sail
    - pheatmap

### Installation

Install the required R packages:

```r
install.packages(c("ggplot2", "caret", "glmnet", "data.table", "e1071", "class", "randomForest", "themis", "groupdata2", "sail", "pheatmap"))
```

### Usage
#### Data Preparation

To prepare the datasets:

Split the dataset into 4 folds and generate training and test sets based on different sampling strategies.

```r
Rscript build_datasets.R
```

This will create the necessary datasets for classification.

#### Running Classifications

To run the classification tasks, use the run_classification.R script. Example usage:

```bash
Rscript run_classification.R 1 V1 class rapide TRUE
```

This command runs the classification for fold 1 using the "V1" sampling strategy at the "class" level with the "rapide" method (logistic regression with lasso constraint) and computes Shapley values.

#### Generating Figures

To generate the figures used in the paper:

```r
Rscript plots_for_figures.R
```

Ensure that the classification results are available as this script will process and visualize the data.

#### Scripts Description

- myUtils.R: Contains utility functions for exporting graphics.
- functions_for_classif.R: Includes all functions necessary for dataset preparation and running classification tasks.
- build_datasets.R: Splits the dataset into training and test sets using various sampling strategies (V1, V3, SMOTE).
- run_classification.R: Executes classification methods (Rapid, KNN, Random Forest, SVM) on the prepared datasets, evaluates performance, and optionally computes Shapley values for explainability.
- functions_for_figures.R: Functions to assist in generating the figures from the classification results.
- plots_for_figures.R: Generates plots used in the publication.
- cmd_launch_analysis.sh: Example bash script to launch build_datasets.R and run_classification.R on a cluster using SLURM.








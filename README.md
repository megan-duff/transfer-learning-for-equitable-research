# Transfer Learning for Equitable Research

This repository contains code, tutorials, and data simulations to support transfer learning approaches in genetics research. The tools provided here help process genetic data, simulate phenotype data, and analyze simulated datasets. The repository also includes a toy example and UK Biobank data analysis tools.

## Table of Contents
- [Overview](#overview)
- [Installation](#installation)
- [Directory Structure](#directory-structure)
- [Usage](#usage)
- [Tutorials](#tutorials)

## Overview

This repository includes:
- **Genetic and phenotype data simulation** tools
- **Data processing** scripts for genetic and phenotypic data
- **Simulated data analysis** comparing transfer learning against benchmarking models
- **UK Biobank analysis** comparing transfer learning against benchmarking models
- **Tutorials** on using the provided tools

## Installation

### Requirements

The repository requires the following software and packages:

- `R` (for phenotype simulation)
- `PLINK` (for genetic data analysis)
- `Python 3.8+` (for machine learning models and data processing)
- `TensorFlow` (for deep learning models)
- `pandas`, `numpy`, `scikit-learn` (for data handling)

To clone the repository:

```bash
git clone https://github.com/megan-duff/transfer-learning-for-equitable-research.git
```

Ensure all dependencies are installed before running the code:

```bash
pip install -r requirements.txt
```

## Directory Structure

```plaintext
.
├── README.md
├── Softwares
│   └── genomic-cnn-transfer-learning
│   └── genomic-cnn
│   └── genomic-mlp-transfer-learning
│   └── genomic-mlp
│   └── phenotype_simulators
├── data-processing
│   └── archive
│   └── compute-allele-frequencies
│   └── create-genome-wide-files
├── genetic-data-simulation
├── phenotype-data-simulation
│   └── additive-and-interaction-effects-phenotypes
│   └── additive-effects-only-phenotypes
│   └── archive
│   └── train-val-test-file-generation
├── simulated-data-analysis
│   └── pre-process-data
├── toy_example
├── tutorials
├── ukb-data-analysis
│   └── Scripts and analysis for working with the UK Biobank dataset

```

### Description of Key Directories

- **Softwares**: Contains various softwares developed from this work.
  - **genomic-cnn-transfer-learning**: Software to train the transfer learning model with a CNN architecture. 
  - **genomic-cnn**: Software to train the source model with a CNN architecture. 
  - **genomic-mlp-transfer-learning**: Software to train the transfer learning model with a MLP architecture. 
  - **genomic-mlp**: Software to train the source model with a MLP architecture. 
  - **phenotype_simulators**: Contains R scripts to simulate additive and interaction effects phenotypes for both the source and target set.

- **data-processing**: Contains scripts for cleaning and computing necessary metrics, primarily from PLINK binary files.
  - **archive**: Rough-draft/extra scripts used to develop this section of the work.  
  - **compute-allele-frequencies**: Scripts used to compute allele frequencies for a variety of datasets. 
  - **create-genome-wide-files**: Scripts used to convert chromosome specific files to genome wide files. 

- **genetic-data-simulation**: Scripts to simulate genetic datasets, used to train and validate models.

- **phenotype-data-simulation**: R scripts to simulate two types of phenotype based on genetic input, phenotypes with additive effects only and phenotypes with additive and interaction effects.
  - **additive-and-interaction-effects-phenotypes**: Scripts used to simulate phenotypes with additive and interaction effects. 
  - **additive-effects-only-phenotypes**: Scripts used to simulate phenotypes with additive effects. 
  - **archive**: Rough-draft/extra scripts used to develop this section of the work. 
  - **train-val-test-file-generation**: Scripts used to break up files into training, validation, and test set specific files. 

- **simulated-data-analysis**: Scripts for trainnig the source and transfer learning models, along with scripts to produce necessary metrics and plots. 
  - **pre-process-data**: Scripts used to run GWAS analysis. 

- **toy_example**: A simple, ready-to-run example. One of the tutoritals uses this dataset to demonstrate the entire pipeline from phenotype simulation to training source and transfer learning models.

- **tutorials**: Markdown tutorials providing detailed instructions for data simulation and analysis.

- **ukb-data-analysis**: Scripts to run the pipeline on the UK Biobank data.

## Usage

To simulate data and run the analysis pipeline:

1. **Simulate Genetic and Phenotype Data**
   - Use scripts in `genetic-data-simulation` and `phenotype-data-simulation` to create a simulated dataset.

2. **Process the Data**
   - Use the `data-processing` scripts to prepare your data for analysis.

3. **Run the Toy Example**
   - Navigate to the `tutorial` folder and follow the tutorial "pheno-sim-transfer-learning-models.ipynb" to run a complete simulation-to-analysis workflow.

4. **Run source and transfer learning method**
   - Use the general scripts in `simulated-data-analysis` or the softwares provided in `Softwares`.

## Tutorials

Detailed tutorials for data processing and genetic data simulation are available in the `tutorials` directory. These include:

- **Data Processing and Genetic Data Simulation Tutorial**: Explains how to clean and prepare genetic datasets, then using HAPGEN2 to simulate genetic data. 

- **Phenotypic Simulation and Source/Transfer Learning Tutorial**: Walkthrough on generating phenotypic data. Using simulated genetic and phenotypic data, source and transfer learning methods are trained. 

# Oral Microbiome Dysbiosis, Salivary Oestradiol Deficiency, and Genetic Risk Factors in Early Alzheimer’s Disease in Postmenopausal Women

This repository contains the datasets, bioinformatics pipelines, processing scripts, and statistical models associated with the study of the "Triple Hit" multifactorial aetiology of early Alzheimer's disease (AD) in postmenopausal women. The study integrates clinical periodontology, localised endocrinology, multi-omic profiling, and genetic susceptibility screening.

## Repository Overview

- **`data/raw/`**: Raw sequence tables, ASV count matrices, and taxonomy assignments generated via QIIME2 / DADA2 processing pipelines.
- **`data/processed/`**: Curated and aggregated clinical, demographic, endocrinological, and periodontal baseline parameters.
- **`scripts/`**: R scripts utilised for baseline clinical analyses, microbiome compositional clustering (PCoA, PERMANOVA, DESeq2), and multivariable logistic regression models evaluating intersectional risk factors.

## Methodological Summary

### 1. Microbiome Analysis (`02_microbiome_analysis.R`)
16S rRNA gene sequencing (V3-V4 hypervariable regions) was conducted on subgingival plaque samples. Bioinformatics steps, including quality filtering, denoising, and merging, were completed using QIIME2 to produce high-resolution Amplicon Sequence Variants (ASVs). Downstream analyses involved Principal Coordinate Analysis (PCoA) using Bray-Curtis dissimilarity matrices and PERMANOVA for cohort differentiation. Differential abundance was determined via `DESeq2`.

### 2. Endocrine Quantification
Free, biologically active 17β-oestradiol concentrations were extracted from unstimulated whole saliva. Oestradiol was quantified utilising Isotope-Dilution Liquid Chromatography-Tandem Mass Spectrometry (LC-MS/MS) with Dansyl Chloride derivation to improve ionization and sensitivity (LLOQ 0.5 pg/mL).

### 3. Clinical and Genetic Assessment (`01_baseline_analysis.R` & `03_multivariable_model.R`)
Genomic DNA isolated from buccal epithelial cells was used for qPCR-based genotyping of the APOE ε4 allele (rs429358 and rs7412) and the TREM2 R47H variant (rs75932628). Multivariable logistic regression modelling evaluated the interactive (Triple Hit) effects of:
1. Localised oestradiol withdrawal
2. High relative abundance of *Porphyromonas gingivalis*
3. The presence of the APOE ε4 allele

## Data Availability
The provided clinical metadata, microbiome taxonomy, and count matrices serve to ensure methodological transparency and computational reproducibility.

## Contact
Yusuf Efe (Pasha) Kivilcim 
contact@kivilcim.me - yek22@cam.ac.uk

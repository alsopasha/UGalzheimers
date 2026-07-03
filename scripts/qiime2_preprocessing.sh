#!/bin/bash
# ==============================================================================
# Early Alzheimer's Data Analysis - 16S Microbiome Pipeline
# Author: Yusuf Efe Kivilcim
# Institution: University of Cambridge
#
# Description: QIIME 2 (v2023.5) pipeline for processing raw V3-V4 Illumina 
# amplicon data into an ASV feature table. 
# ==============================================================================

# 1. Import raw FASTQ data into QIIME 2 artifact
qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path ../data/raw/raw_fastq_manifest.csv \
  --output-path demux_sequences.qza \
  --input-format PairedEndFastqManifestPhred33V2

# 2. Denoising and ASV generation using DADA2
# (Parameters adjusted for V3-V4 amplicons 2x300bp)
qiime dada2 denoise-paired \
  --i-demultiplexed-seqs demux_sequences.qza \
  --p-trunc-len-f 280 \
  --p-trunc-len-r 260 \
  --p-trim-left-f 0 \
  --p-trim-left-r 0 \
  --o-table asv_table.qza \
  --o-representative-sequences rep_seqs.qza \
  --o-denoising-stats dada2_stats.qza \
  --p-n-threads 8

# 3. Taxonomy Classification 
# Using a Naive Bayes classifier trained on SILVA v138 (V3-V4 region)
qiime feature-classifier classify-sklearn \
  --i-classifier silva-138-99-v3-v4-classifier.qza \
  --i-reads rep_seqs.qza \
  --o-classification taxonomy.qza

# 4. Exporting data for R analysis
qiime tools export \
  --input-path asv_table.qza \
  --output-path exported_feature_table

# The resulting biom file is converted to CSV (microbiome_asv_table.csv) for R integration.

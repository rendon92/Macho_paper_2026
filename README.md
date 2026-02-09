# Computational Analysis Pipeline: Zebrafish Eye Development Multiomics Study

This repository contains the custom computational pipelines used for the analysis of single-cell and single-nucleus multiome data in the study:

**"Single-cell multiomics reconstructs gene regulatory networks during early zebrafish eye development"**  
*bioRxiv 2026.01.31.703039*

## Repository Overview

This repository is organized into two main analytical pipelines:

### 1. **Single-cell RNA-seq Analysis** (`/single_cell_pipeline/`)
- **Input**: Processed gene expression matrices from 10X Genomics Cell Ranger
- **Methods**: Standard Seurat workflow for quality control, normalization, integration, clustering, and visualization
- **Key Feature**: Anatomical annotation system using zebrafish anatomical ontology (ZFIN) with hypergeometric testing
- **Output**: Annotated cell clusters, UMAP visualizations, marker genes, and developmental trajectories

### 2. **Single-nucleus Multiome Analysis** (`/multiome_pipeline/`)
- **Input**: Paired RNA-seq and ATAC-seq data from 10X Genomics Multiome platform
- **Methods**: Integration of transcriptomic and chromatin accessibility data using Seurat and Signac
- **Output**: Joint embeddings, differentially accessible regions, motif enrichment, and coordinated gene-regulatory analyses

### **Note on Gene Regulatory Network Inference**
The gene regulatory network (GRN) inference presented in the manuscript was performed using **MatchaiRen**, a separate tool developed for this study. The MatchaiRen pipeline is available in its own repository:  
🔗 **[MatchaiRen_v0.1](https://github.com/rendon92/MatchaiRen_v0.1)**

## Quick Start

### Prerequisites
- R (≥ 4.1.0)
- Key R packages: Seurat, Signac, ggplot2, dplyr
- Python (for some preprocessing scripts)
- 10X Genomics Cell Ranger outputs

### Basic Usage
1. **Single-cell RNA-seq analysis**:
   ```bash
   cd single_cell_pipeline/downstream_analysis_R
   Rscript main.R --input path/to/cellranger/output --output results/

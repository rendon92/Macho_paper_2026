
---

```markdown
# Single-cell RNA-seq Analysis Pipeline

This directory contains the downstream analysis pipeline for single-cell RNA-seq data from zebrafish embryonic eye and forebrain development.

## Analysis Workflow

### Input Requirements
- **Processed Cell Ranger outputs**: `filtered_feature_bc_matrix/` directories for each sample
- **Sample metadata**: CSV file with sample information (stage, genotype, etc.)

### Pipeline Steps

1. **Quality Control & Filtering** (`main.R`)
   - Filter cells based on: nFeature_RNA, nCount_RNA, percent.mt
   - Remove low-quality cells

2. **Normalization & Integration** (`main.R`)
   - SCTransform normalization for each sample
   - Integration across samples/timepoints using Seurat anchors
   - Batch effect correction

3. **Dimensionality Reduction & Clustering** (`main.R`)
   - PCA on variable features
   - UMAP for visualization
   - Louvain clustering at multiple resolutions

4. **Cell Type Annotation** (`annotation.R`)
   - **Anatomical ontology-based annotation** using ZFIN database
   - Hypergeometric testing for enriched anatomical terms
   - Automatic assignment of the most enriched anatomical structure

5. **Differential Expression & Marker Identification**
   - Find markers for each cluster
   - Identify stage-specific genes
   - Export results for downstream analysis

## Anatomical Annotation System

The key feature in this pipeline is our anatomical annotation approach implemented in `annotation.R`:

### How It Works
1. **Input**: List of significantly overexpressed genes for each cluster
2. **Reference Database**: Zebrafish anatomical ontology from ZFIN, linking genes to anatomical structures where they are expressed
3. **Statistical Test**: Hypergeometric test evaluates enrichment of anatomical terms in each cluster's gene set
4. **Scoring**: Each anatomical term receives a p-value (adjusted for multiple testing)
5. **Assignment**: The cluster is annotated with the anatomical term showing the most significant enrichment

### Output
- Primary annotation: The most enriched anatomical structure for each cluster
- Full ranking: All anatomical terms with their enrichment scores (useful for ambiguous clusters)
- Visualization: UMAP plots colored by anatomical annotation


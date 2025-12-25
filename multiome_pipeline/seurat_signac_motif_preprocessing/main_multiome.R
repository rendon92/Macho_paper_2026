#loading packages

library(Seurat)
library(Signac)
library(GenomicRanges)
library(MACSr)
library(dplyr)
library(cowplot)
library(lattice)
library(SeuratData)
library(patchwork)
library(ggalluvial)
library(future)
library(stringr)
library(limma)
library(Matrix)
library(qlcMatrix)
library(JASPAR2020)
library(TFBSTools)
options(future.globals.maxSize = 8000 * 1024^2)
set.seed(2807)


# -----------------------------
# Paths
# -----------------------------

counts <- Read10X("/data/scMultiome/filtered_feature_bc_matrix/")
fragpath <- "/data/scMultiome/atac_fragments.tsv.gz"
fa <-  Rsamtools::FaFile("/data/scMultiome/genomes/danio_rerio/GRCz11_NCBI/fasta/genome.fa")

annotation <- rtracklayer::import('/data/scMultiome/filtered_Danio_rerio.GRCz11.107.gtf')

#ZFIN ANNOTATIONS
zf_annot <- read.delim("/data/scRNAseq/seurat/zf_annot/wt_expression_fish.txt", header=FALSE)
zf_comb <- subset(zf_annot,
                  (V3 %in% c("WT","AB","AB/TU","TU")) &
                    (V8 %in% c("Segmentation:14-19 somites","Segmentation:20-25 somites")) &
                    (V9 %in% c("Segmentation:14-19 somites","Segmentation:20-25 somites")))
fc_cut <- 0.5; p_cut <- 1e-5
zf_comb <- distinct(zf_comb, V2, V5, V8, V9, .keep_all = TRUE)
zf_comb <- zf_comb[, c("V5", "V2")]
colnames(zf_comb) <- c("Anatomical", "Genename")

## Running the pipeline with the whole dataset

dr18 <- CreateSeuratObject(
  counts = counts$`Gene Expression`,
  assay = "RNA"
)

dr18[["ATAC"]] <- CreateChromatinAssay(
  counts = counts$Peaks,
  sep = c(":", "-"),
  fragments = fragpath,
  annotation = annotation
)

#quality control ATAC

DefaultAssay(dr18) <- "ATAC"
dr18 <- NucleosomeSignal(dr18)
dr18 <- TSSEnrichment(dr18)

# filter out low quality cells
dr18 <- subset(
  x = dr18,
  subset = nCount_ATAC < 100000 &
    nCount_RNA < 25000 &
    nCount_ATAC > 1000 &
    nCount_RNA > 1000 &
    nucleosome_signal < 2 &
    TSS.enrichment > 1
)

gc()

DefaultAssay(dr18) <- "RNA"
dr18 <- NormalizeData(object = dr18, normalization.method = "LogNormalize", scale.factor = 10000)
dr18 <- FindVariableFeatures(object = dr18, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5, nfeatures = 2000)


all.genes_combined <- rownames(dr18)
dr18 <- ScaleData(dr18, verbose = FALSE, features = all.genes_combined)
dr18 <- RunPCA(dr18, features = VariableFeatures(object = dr18))
dr18 <- FindNeighbors(dr18, reduction = "pca", dims = 1:30)
dr18 <- FindClusters(dr18, resolution = 5)
dr18 <- RunUMAP(dr18, reduction = "pca", dims = 1:30)

# Doublet identification
library(SingleCellExperiment)
library(scDblFinder)

dr18_sce <- as.SingleCellExperiment(dr18)
colLabels(dr18_sce) <- dr18_sce$seurat_clusters

set.seed(2807)
dr18_sce <- scDblFinder(dr18_sce)


dr18$doublet_class <- dr18_sce$scDblFinder.class


# eliminar doublets
dr18 <- subset(dr18, subset = doublet_class == "singlet")

# Invoke anatomical enrichment function


cluster_anatomical_enrichment <- function(
    seurat_object,
    annotation_data,
    log2FC = 0.5,
    exclude_terms = NULL
) {
  
  library(dplyr)
  
  ## -----------------------------
  ## Prepare data
  ## -----------------------------
  genes_presentes <- rownames(seurat_object)
  
  annotation_data <- annotation_data %>% 
    filter(Genename %in% genes_presentes) %>% 
    distinct(Genename, Anatomical)
  
  background_genes <- unique(annotation_data$Genename)
  tfreq <- table(annotation_data$Anatomical)
  total_background <- length(background_genes)
  
  results_list <- list()
  
  ## -----------------------------
  ## Loop per cluster
  ## -----------------------------
  for (cluster_id in levels(seurat_object)) {
    
    message("Procesando cluster: ", cluster_id)
    
    markers <- FindMarkers(
      seurat_object,
      ident.1 = cluster_id,
      logfc.threshold = log2FC,
      min.pct = 0.05
    )
    
    markers$pct_dif <- markers$pct.1 - markers$pct.2
    markers <- markers[markers$pct_dif >= 0.15, ]
    if (nrow(markers) < 8) next
    
    markers$genename <- rownames(markers)
    
    markers <- markers %>%
      filter(avg_log2FC > 0.5) %>%
      mutate(corrected = avg_log2FC * pct.1) %>%
      filter(corrected > 0.5)
    
    gene_set <- unique(markers$genename)
    
    enriched_terms <- annotation_data %>%
      filter(Genename %in% gene_set)
    
    if (nrow(enriched_terms) == 0) next
    
    efreq <- table(enriched_terms$Anatomical)
    term_names <- names(efreq)
    
    df <- data.frame(
      cluster = cluster_id,
      Anatomical = term_names,
      count = as.integer(efreq),
      geneRatio = as.numeric(efreq) / length(gene_set),
      BgRatio = as.numeric(tfreq[term_names]) / total_background,
      stringsAsFactors = FALSE
    )
    
    df$pvalue <- mapply(
      function(k, K) {
        phyper(k - 1, K, total_background - K, length(gene_set), lower.tail = FALSE)
      },
      k = df$count,
      K = tfreq[df$Anatomical]
    )
    
    if (!is.null(exclude_terms)) {
      df <- df[!df$Anatomical %in% exclude_terms, ]
    }
    
    df$padj <- p.adjust(df$pvalue, method = "BH")
    
    genes_by_term <- enriched_terms %>%
      group_by(Anatomical) %>%
      summarise(genes = paste(unique(Genename), collapse = "/"),
                .groups = "drop")
    
    df <- left_join(df, genes_by_term, by = "Anatomical")
    
    results_list[[cluster_id]] <- df[order(df$padj), ]
  }
  
  ## -----------------------------
  ## Final Table
  ## -----------------------------
  final_result <- bind_rows(results_list)
  
  ## Best term per cluster
  final_min_padj <- final_result %>%
    arrange(cluster, padj) %>%
    distinct(cluster, .keep_all = TRUE)
  
  ## -----------------------------
  ## Annotate Seurat
  ## -----------------------------
  annotation_df <- final_min_padj %>%
    select(cluster, Anatomical) %>%
    mutate(cluster = as.character(cluster))
  
  seurat_object@meta.data$annotation <-
    annotation_df$Anatomical[
      match(
        as.character(seurat_object@meta.data$RNA_snn_res.5),
        annotation_df$cluster
      )
    ]
  
  ## -----------------------------
  ## Return everything
  ## -----------------------------
  return(list(
    seurat_object = seurat_object,
    enrichment_table = final_result,
    best_annotation_per_cluster = final_min_padj
  ))
}

# Annotate cells

res <- cluster_anatomical_enrichment(
  seurat_object = dr18,
  annotation_data = zf_comb
)
## -----------------------------
# Recover objects 
## -----------------------------
dr18 <- res$seurat_object
final_result <- res$enrichment_table
final_min_padj <- res$best_annotation_per_cluster


## -----------------------------
### Filter out cells annotated as eye and neural origin. Obtain cell barcodes
## -----------------------------
idents <- c("hindbrain", "immature eye", "midbrain hindbrain boundary", "forebrain",
            "diencephalon", "telencephalon", "optic vesicle", "optic stalk", "neural tube", "lens placode")
dr18_subset <- subset(dr18, annotation %in% idents)

ids_subset <- colnames(dr18_subset)


#####Rerun downstream analysis with multimodal mode



# Rebuild initial object
dr18 <- CreateSeuratObject(
  counts = counts$`Gene Expression`,
  assay = "RNA"
)

dr18[["ATAC"]] <- CreateChromatinAssay(
  counts = counts$Peaks,
  sep = c(":", "-"),
  fragments = fragpath,
  annotation = annotation
)
#subset the data with eye and neural cells
dr18_s <- subset(dr18, cells = ids_subset ) 


DefaultAssay(dr18_s) <- "RNA"
dr18_s <- NormalizeData(object = dr18_s, normalization.method = "LogNormalize", scale.factor = 10000)
dr18_s <- FindVariableFeatures(object = dr18_s, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5, nfeatures = 2000)
all.genes_combined <- rownames(dr18_s)
dr18_s <- ScaleData(dr18_s, verbose = FALSE, features = all.genes_combined)
dr18_s <- RunPCA(dr18_s, features = VariableFeatures(object = dr18_s))



DefaultAssay(dr18_s) <- "ATAC"
dr18_s <- NucleosomeSignal(dr18_s)
dr18_s <- TSSEnrichment(dr18_s)
dr18_s <- FindTopFeatures(dr18_s, min.cutoff = 5)
dr18_s <- RunTFIDF(dr18_s)
dr18_s <- RunSVD(dr18_s)

dr18_s <- FindMultiModalNeighbors(
    object = dr18_s,
    reduction.list = list("pca", "lsi"), 
    dims.list = list(1:50, 2:40),
    modality.weight.name = "RNA.weight",
    verbose = TRUE
)

DefaultAssay(dr18_s) <- "RNA"

dr18_s <- FindNeighbors(dr18_s, reduction = "pca", dims = 1:30)
dr18_s <- FindClusters(dr18_s, resolution = 5)


dr18_s <- RunUMAP(
    object = dr18_s,
    nn.name = "weighted.nn",
    assay = "RNA",
    verbose = TRUE
)


DimPlot(dr18_s, label = TRUE, repel = TRUE, reduction = "umap", group.by = "seurat_clusters", pt.size = 1) + NoLegend()



## Differenrial expression analysis

res_subset <- cluster_anatomical_enrichment(
  seurat_object = dr18_s,
  annotation_data = zf_comb
)

# Recover objects
dr18_s <- res_subset$seurat_object
final_result_subset <- res_subset$enrichment_table
final_min_padj_subset <- res_subset$best_annotation_per_cluster
DimPlot(dr18_s, label = TRUE)

#Find differentially expressed genes per cluster
all_DEGs <- FindAllMarkers(dr18_s, logfc.threshold = 0.25, min.pct = 0.05)

#Refine annotations
dr18_s@meta.data$annotation[dr18_s@meta.data$RNA_snn_res.5 %in% as.character(c(1,5,6,8,20,53))] <- "neural retina"
dr18_s@meta.data$annotation[dr18_s@meta.data$RNA_snn_res.5 %in% as.character(c(3,35,45))] <- "RPE"


# Find differentially expressed genes per domain
Idents(dr18_s) <- dr18_s@meta.data$annotation

all_DEGs_domains <- FindAllMarkers(dr18_s, logfc.threshold = 0.25, min.pct = 0.05)

DefaultAssay(dr18_s) <- "ATAC"
all_DOCRs_domains <- FindAllMarkers(dr18_s, logfc.threshold = 0.25, min.pct = 0.05, test.use = "LR", latent.vars = "nCount_ATAC")




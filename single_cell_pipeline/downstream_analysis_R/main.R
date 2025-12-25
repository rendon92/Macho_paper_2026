# main_scRNA_pipeline.R

library(Seurat)
library(dplyr)
library(SCpubr)
library(SeuratData)
library(SeuratDisk)
library(Matrix.utils)

source("annotation.R")

set.seed(2807)

# -----------------------------
# Paths
# -----------------------------
zf_annot <- read.delim("wt_expression_fish.txt", header=FALSE)

dataset_dirs <- list(
  "15hpf" = "run_count_JR_SCR_15/outs/filtered_feature_bc_matrix/",
  "18hpf" = "run_count_JR_SCR_18/outs/filtered_feature_bc_matrix/",
  "23hpf" = "run_count_JR_SCR_23/outs/filtered_feature_bc_matrix/"
)

# -----------------------------
# Hard-coded cluster annotation per time
# -----------------------------
cluster_annotations <- list(
  "15hpf" = list(
    "RPE" = c("15","17","26"),
    "Neural retina" = c("0","1","4","7","9","10","19","20","22","23","25","30"),
    "Telencephalon" = c("8","18","28","34"),
    "Optic Stalk" = c("16","29"),
    "Diencephalon" = c("2","5","11","24"),
    "Forebrain" = c("33"),
    "Neural tube" = c("32"),
    "Trigeminal placode" = c("37"),
    "Lens placode" = c("3","36"),
    "Periderm" = c("12","14","31","38","39"),
    "Proliferative region" = c("6","13","21","27","35","40")
  ),
  "18hpf" = list(
  "Neural retina" = c("0","2","3","9","13","18","20","23","33"),
  "RPE" = c("8","19"),
  "Telencephalon" = c("11"),
  "Optic Stalk" = c("10","15","17","30","38"),
  "Diencephalon" = c("4","7","12","14"),
  "Lens placode" = c("24","28","34"),
  "Trigeminal placode" = c("31","32"),
  "Midbrain" = c("39"),
  "Periderm" = c("16","25","27","35","36","37"),
  "Proliferative region" = c("5","6","21", "22"),
  "Hindbrain" = c("1","29"),
  "Primary neuron" = c("26")
),
"23hpf" = list(
    "RPE" = c("6","9","19","20"),
    "Neural retina" = c("0","3","8","10","15","16","17"),
    "Telencephalon" = c("14"),
    "Optic Stalk" = c("11"),
    "Diencephalon" = c("1","4","12"),
    "Lens placode" = c("18"),
    "Proliferative region" = c("2","5","13","21","22"),
    "Neural crest" = c("7")
  )
)

# -----------------------------
# Function to run pipeline
# -----------------------------
process_timepoint <- function(timepoint, dataset_dir, zf_annot, cluster_annotations) {
  message("Processing ", timepoint)
  
  # -----------------------------
  # Subset annotation
  # -----------------------------
  if (timepoint == "15hpf") {
    zf_comb <- subset(zf_annot,
                      (V3 %in% c("WT","AB","AB/TU","TU")) &
                        (V8 %in% c("Segmentation:10-13 somites","Segmentation:14-19 somites")) &
                        (V9 %in% c("Segmentation:10-13 somites","Segmentation:14-19 somites")))
    fc_cut <- 0.5; p_cut <- 1
  } else if (timepoint == "18hpf") {
    zf_comb <- subset(zf_annot,
                      (V3 %in% c("WT","AB","AB/TU","TU")) &
                        (V8 %in% c("Segmentation:14-19 somites","Segmentation:20-25 somites")) &
                        (V9 %in% c("Segmentation:14-19 somites","Segmentation:20-25 somites")))
    fc_cut <- 0.5; p_cut <- 1e-5
  } else if (timepoint == "23hpf") {
    zf_comb <- subset(zf_annot,
                      (V3 %in% c("WT","AB","AB/TU","TU")) &
                        (V8 %in% c("Segmentation:26+ somites","Pharyngula:Prim-5")) &
                        (V9 %in% c("Segmentation:26+ segmentation","Pharyngula:Prim-5")))
    fc_cut <- 1; p_cut <- 1e-10
  }
  
  zf_comb <- distinct(zf_comb, V2, V5, V8, V9, .keep_all = TRUE)
  
  # -----------------------------
  # Load 10X data
  # -----------------------------
  data <- Read10X(data.dir = dataset_dir)
  seurat_obj <- CreateSeuratObject(counts = data, min.cells = 3, min.features = 200,
                                   project = paste0(timepoint, "_10x"), assay = "RNA")
  seurat_obj@meta.data$barcode <- rownames(seurat_obj@meta.data)
  seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^mt-")
  
  # -----------------------------
  # Preprocessing
  # -----------------------------
  seurat_obj <- SCTransform(seurat_obj, vars.to.regress = "percent.mt", verbose = FALSE)
  seurat_obj <- RunPCA(seurat_obj)
  seurat_obj <- RunUMAP(seurat_obj, dims = 1:30, verbose = FALSE)
  seurat_obj <- FindNeighbors(seurat_obj, dims = 1:30, verbose = FALSE)
  seurat_obj <- FindClusters(seurat_obj, verbose = FALSE, resolution = 2)
  
  DimPlot(seurat_obj, label = TRUE)
  
  # -----------------------------
  # Select clusters with marker genes
  # -----------------------------
  all_markers <- FindAllMarkers(seurat_obj, logfc.threshold = 0.5, min.pct = 0.1)
  genes_interest <- c("vsx2","bhlhe40","vax1","emx3","nkx2.4a","foxe3")
  
  clusters <- unique(as.numeric(as.character(
    filter(all_markers, avg_log2FC > fc_cut & p_val_adj <= p_cut & gene %in% genes_interest)$cluster
  )))
  
  seurat_obj <- subset(seurat_obj, cells = colnames(subset(seurat_obj, seurat_clusters %in% clusters)))
  
  # -----------------------------
  # Re-run clustering on subset
  # -----------------------------
  seurat_obj <- SCTransform(seurat_obj, vars.to.regress = "percent.mt", verbose = FALSE)
  seurat_obj <- RunPCA(seurat_obj)
  seurat_obj <- RunUMAP(seurat_obj, dims = 1:30, verbose = FALSE)
  seurat_obj <- FindNeighbors(seurat_obj, dims = 1:30, verbose = FALSE)
  seurat_obj <- FindClusters(seurat_obj, verbose = FALSE, resolution = 5)
  
  DimPlot(seurat_obj, label = TRUE)
  
  # -----------------------------
  # Annotate cells
  # -----------------------------
  seurat_obj_a <- cell_annotation(seurat_obj, log2FC = 0.5, annotation_data = zf_comb)[[1]]
  
  # Hard-coded cluster names
  for (cl_name in names(cluster_annotations[[timepoint]])) {
    seurat_obj_a@meta.data$cluster[seurat_obj_a@meta.data$SCT_snn_res.5 %in% cluster_annotations[[timepoint]][[cl_name]]] <- cl_name
  }
  
  Idents(seurat_obj_a) <- seurat_obj_a@meta.data$cluster
  
  DimPlot(seurat_obj_a, label = TRUE, repel = TRUE, reduction = "umap", pt.size = 1) + NoLegend()
  FeaturePlot(seurat_obj_a, features = genes_interest)
  
  markers <- FindAllMarkers(seurat_obj_a, logfc.threshold = 0.5, min.pct = 0.1)
  markers$orig.ident <- timepoint
  
  return(list(seurat = seurat_obj_a, markers = markers))
}

# -----------------------------
# Run pipeline for all timepoints
# -----------------------------
results <- lapply(names(dataset_dirs), function(tp) {
  process_timepoint(tp, dataset_dirs[[tp]], zf_annot, cluster_annotations)
})
names(results) <- names(dataset_dirs)

# -----------------------------
# Prepare combined annotations
# -----------------------------
zf_annot <- read.delim("wt_expression_fish.txt", header=FALSE)

# Subset combining all somites/stages of interest
zf_comb <- subset(zf_annot,
                  (V3 %in% c("WT","AB","AB/TU","TU")) &
                    (V8 %in% c("Segmentation:10-13 somites",
                               "Segmentation:14-19 somites",
                               "Segmentation:20-25 somites",
                               "Segmentation:26+ somites",
                               "Pharyngula:Prim-5")) &
                    (V9 %in% c("Segmentation:10-13 somites",
                               "Segmentation:14-19 somites",
                               "Segmentation:20-25 somites",
                               "Segmentation:26+ somites",
                               "Pharyngula:Prim-5")))

zf_comb <- zf_comb[, c("V5", "V2")]
colnames(zf_comb) <- c("Anatomical", "Genename")

# -----------------------------
# Function to run enrichment by timepoint
# -----------------------------
run_anatomical_for_timepoints <- function(results_list, annotation_data) {
  res_ao <- list()
  
  for (tp in names(results_list)) {
    message("Processing anatomical enrichment for", tp)
    seurat_obj <- results_list[[tp]][["seurat"]]
    Idents(seurat_obj) <- seurat_obj@meta.data$seurat_clusters
    
    df <- cluster_anatomical_enrichment(seurat_obj, annotation_data)
    res_ao[[tp]] <- df
  }
  
  return(res_ao)
}

# -----------------------------
# Execute for all timepoints
# -----------------------------
anatomical_results <- run_anatomical_for_timepoints(results, zf_comb)

# -----------------------------
# Access to results
# -----------------------------

head(anatomical_results[["15hpf"]])
head(anatomical_results[["18hpf"]])
head(anatomical_results[["23hpf"]])
# -----------------------------
# Save markers
# -----------------------------
for (tp in names(results)) {
  write.csv(results[[tp]]$markers, paste0("markers_", tp, ".csv"), row.names = FALSE)
}



dr15_a2 <- results[["15hpf"]][["seurat"]]
dr18_a2 <- results[["18hpf"]][["seurat"]]
dr23_a2 <- results[["23hpf"]][["seurat"]]


df15 <- FindAllMarkers(dr15_a2, logfc.threshold = 1, min.pct = 0.1)
df18 <- FindAllMarkers(dr18_a2, logfc.threshold = 1, min.pct = 0.1)
df23 <- FindAllMarkers(dr23_a2, logfc.threshold = 1, min.pct = 0.1)
df15$orig.ident <- "15hpf"
df18$orig.ident <- "18hpf"
df23$orig.ident <- "23hpf"

rownames(df15) <- NULL
rownames(df18) <- NULL
rownames(df23) <- NULL
df <- rbind(df15, df18, df23)

colors_clusters <- c(
  "Diencephalon" = "#ffd7a8a1",
  "Forebrain" = "#97a900ff",
  "Lens placode" = "#2fb600ff",
  "Neural retina" = "#569cc9ff",
  "Neural tube" = "#00a6ffff",
  "Optic Stalk" = "#00c0b7ff",
  "Periderm" = "#c8ffc3ff",
  "Proliferative region" = "#db8e00ff",
  "RPE" = "#ad16c699",
  "Telencephalon" = "#fba6a7ff",
  "Trigeminal placode" = "#ff63b6ff",
  "Neural crest" = "#6fa690ff",
  "Primary neuron" = "#b19641ff",
  "Midbrain" = "#b19641ff",
  "Hindbrain" = "#b19641ff"
)


meta.data <- rbind(dr15_a2@meta.data, dr18_a2@meta.data, dr23_a2@meta.data)

# Seurat's object list by stage
stages_list <- list(
  "15hpf" = dr15_a2,
  "18hpf" = dr18_a2,
  "23hpf" = dr23_a2
)

write.table(df, "SupplementaryDataset_1.tsv", sep = "\t", row.names = FALSE, quote = FALSE)
retinal_markers <- as.data.frame(unique(subset(df, df$cluster == "Neural retina" & df$avg_log2FC >= 0 & df$p_val_adj <= 0.05))$gene)
write.table(retinal_markers, "/data/scMultiome/paper/figure_all_data/figure1/retinal_markers.tsv", sep = "\t", row.names = FALSE, quote = FALSE, col.names = FALSE)

# UMAP by cell type
SCpubr::do_DimPlot(dr15_a2, colors.use = colors_clusters, label = TRUE) + NoLegend()
SCpubr::do_DimPlot(dr18_a2, colors.use = colors_clusters, label = TRUE) + NoLegend()
SCpubr::do_DimPlot(dr23_a2, colors.use = colors_clusters, label = TRUE) + NoLegend()

# UMAP of marker genes per cell type
FeaturePlot(dr15_a2, features =  c("vsx2", "bhlhe40", "vax1", "emx3", "nkx2.4a", "foxe3"))
FeaturePlot(dr18_a2, features =  c("vsx2", "bhlhe40", "vax1", "emx3", "nkx2.4a", "foxe3"))
FeaturePlot(dr23_a2, features =  c("vsx2", "bhlhe40", "vax1", "emx3", "nkx2.4a", "foxe3"))


# transform to seurat v4

convert_to_v4 <- function(obj) {
  DefaultAssay(obj) <- "RNA"
  
  # Extraer counts y metadata
  counts <- GetAssayData(obj, layer = names(obj[["RNA"]]@layers)[1])  # toma la capa original
  meta <- obj@meta.data
  
  # Crear Seurat v4 clásico
  obj_v4 <- CreateSeuratObject(counts = counts, meta.data = meta)
  return(obj_v4)
}

dr15_v4 <- convert_to_v4(dr15_a2)
dr18_v4 <- convert_to_v4(dr18_a2)
dr23_v4 <- convert_to_v4(dr23_a2)

#Merge

dr15_v4$stage <- "15hpf"
dr18_v4$stage <- "18hpf"
dr23_v4$stage <- "23hpf"

seurat_merged_v4 <- merge(dr15_v4, y = list(dr18_v4, dr23_v4), 
                          add.cell.ids = c("15hpf", "18hpf", "23hpf"), 
                          project = "zebrafish_dev")

# Standard preprocessing
seurat_merged_v4 <- SCTransform(seurat_merged_v4, verbose = FALSE)
seurat_merged_v4 <- RunPCA(seurat_merged_v4)
seurat_merged_v4 <- RunUMAP(seurat_merged_v4, dims = 1:30)
seurat_merged_v4 <- FindNeighbors(seurat_merged_v4, dims = 1:30)
seurat_merged_v4 <- FindClusters(seurat_merged_v4, resolution = 0.5)

# Plot merged clusters
SCpubr::do_DimPlot(seurat_merged_v4, colors.use = colors_clusters) + NoLegend()

# Integrate
dr15_v4 <- SCTransform(dr15_v4, verbose = FALSE)
dr18_v4 <- SCTransform(dr18_v4, verbose = FALSE)
dr23_v4 <- SCTransform(dr23_v4, verbose = FALSE)

seurat_list <- list(dr15_v4, dr18_v4, dr23_v4)
features <- SelectIntegrationFeatures(object.list = seurat_list, nfeatures = 3000)
seurat_list <- PrepSCTIntegration(object.list = seurat_list, anchor.features = features)
anchors <- FindIntegrationAnchors(object.list = seurat_list, normalization.method = "SCT", anchor.features = features)
seurat_integrated_v4 <- IntegrateData(anchorset = anchors, normalization.method = "SCT")

# PCA, UMAP, clustering
seurat_integrated_v4 <- RunPCA(seurat_integrated_v4)
seurat_integrated_v4 <- RunUMAP(seurat_integrated_v4, dims = 1:30)
seurat_integrated_v4 <- FindNeighbors(seurat_integrated_v4, dims = 1:30)
seurat_integrated_v4 <- FindClusters(seurat_integrated_v4, resolution = 0.5)

seurat_integrated_v4@meta.data$cluster[seurat_integrated_v4@meta.data$cluster == "proliferative region"] <- "Proliferative region"
Idents(seurat_integrated_v4) <- seurat_integrated_v4@meta.data$cluster

SCpubr::do_DimPlot(sample = seurat_integrated_v4, colors.use = c("Diencephalon" = "#ffd7a8a1", "Forebrain" = "#97a900ff", "Lens placode" = "#2fb600ff", "Neural retina" = "#569cc9ff", "Neural tube" = "#00a6ffff", "Optic Stalk" = "#00c0b7ff", "Periderm" = "#c8ffc3ff", "Proliferative region" = "#db8e00ff", "RPE" = "#ad16c699", "Telencephalon" = "#fba6a7ff", "Trigeminal placode" = "#ff63b6ff", "Neural crest" = "#6fa690ff", "primary neuron" = "#b19641ff", "Midbrain" = "#b19641ff", "Hindbrain" = "#b19641ff")) + NoLegend()





#### MERGING AND SUBSETTING NEURAL RETINA AND PROLIFERATIVE REGION

set.seed(2807)

dr15.data <- Read10X(data.dir = "run_count_JR_SCR_15/outs/filtered_feature_bc_matrix/")
dr15 <- CreateSeuratObject(counts = dr15.data, min.cells = 3, min.features  = 200, project = "dr15_10x", assay = "RNA")

dr18.data <- Read10X(data.dir = "run_count_JR_SCR_18/outs/filtered_feature_bc_matrix/")
dr18 <- CreateSeuratObject(counts = dr18.data, min.cells = 3, min.features  = 200, project = "dr18_10x", assay = "RNA")

dr23.data <- Read10X(data.dir = "run_count_JR_SCR_23/outs/filtered_feature_bc_matrix/")
dr23 <- CreateSeuratObject(counts = dr23.data, min.cells = 3, min.features  = 200, project = "dr23_10x", assay = "RNA")


#Merged the datasets

dr_merged <- merge(dr15, y = c(dr18, dr23), project = "dr_merged")

meta.data_merged <- dr_merged@meta.data

## Load metadata with information of each cell and its annotated type
meta.data <- readRDS("cell_annotation_dr15_dr18_dr23.rds")
meta.data$barcode <- rownames(meta.data)
rownames(meta.data) <- paste(meta.data$barcode, meta.data$orig.ident, sep = "_")
meta.data_merged$barcode <- rownames(meta.data_merged)
meta.data_merged$barcode <- gsub("_.$", "", meta.data_merged$barcode)
rownames(meta.data_merged) <- paste(meta.data_merged$barcode, meta.data_merged$orig.ident, sep = "_")
colnames(dr_merged) <- rownames(meta.data_merged)
dr_merged$barcode <- colnames(dr_merged)


cells <- intersect(rownames(meta.data), rownames(dr_merged@meta.data))
dr_merged_subset <- subset(dr_merged, barcode %in% cells)
meta.data <- subset(meta.data, rownames(meta.data) %in% cells)
dr_merged_subset@meta.data$cluster <- meta.data$cluster


DefaultAssay(dr_merged_subset) <- "RNA"
dr_merged_subset <- NormalizeData(dr_merged_subset, normalization.method = "LogNormalize", scale.factor = 10000)
dr_merged_subset <- FindVariableFeatures(dr_merged_subset, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(dr_merged_subset)
dr_merged_subset <- ScaleData(dr_merged_subset)
dr_merged_subset <- RunPCA(dr_merged_subset, features = VariableFeatures(object = dr_merged_subset))
dr_merged_subset <- RunUMAP(dr_merged_subset, dims = 1:30, verbose = FALSE)
dr_merged_subset <- FindNeighbors(dr_merged_subset, dims = 1:30, verbose = FALSE)
dr_merged_subset <- FindClusters(dr_merged_subset, verbose = FALSE, resolution = 2)
DimPlot(dr_merged_subset, label = TRUE, group.by = "orig.ident")
dr_merged_subset <- JoinLayers(dr_merged_subset)
#DimPlot(dr_merged_subset, label = TRUE, group.by = "cluster")
Idents(dr_merged_subset) <- dr_merged_subset@meta.data$cluster
dr_markers <- FindAllMarkers(dr_merged_subset, logfc.threshold = 0.5, min.pct = 0.05)

dr_merged_subset[["RNA3"]] <- as(object = dr_merged_subset[["RNA"]], Class = "Assay")
DefaultAssay(dr_merged_subset) <- "RNA3"

dr_merged_subset[["RNA"]] <- NULL
dr_merged_subset <- RenameAssays(object = dr_merged_subset, RNA3 = 'RNA')

## In case we wanted to create subsets of just the same annotated cell type through different time points:

sub_dr <- subset(dr_merged_subset, idents = c("Neural retina", "Proliferative region"))

sub_dr <- NormalizeData(object = sub_dr, normalization.method = "LogNormalize", scale.factor = 10000)

sub_dr <- FindVariableFeatures(object = sub_dr, mean.function = ExpMean, dispersion.function = LogVMR, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5, nfeatures = 2000)


all.genes_combined <- rownames(sub_dr)
sub_dr <- ScaleData(sub_dr, verbose = FALSE)


sub_dr <- RunPCA(sub_dr, features = VariableFeatures(object = sub_dr))
gc()
sub_dr <- FindNeighbors(sub_dr, reduction = "pca", dims = 1:30)
gc()
sub_dr <- FindClusters(sub_dr, resolution = 2)
gc()
sub_dr <- RunUMAP(sub_dr, reduction = "pca", dims = 1:30)

DimPlot(sub_dr, label = TRUE, pt.size = 1, repel = TRUE, group.by = "orig.ident") + NoLegend()
Idents(sub_dr) <- sub_dr@meta.data$orig.ident
all_markers <- FindAllMarkers(sub_dr, logfc.threshold = 0.5, min.pct = 0.05)
driver_DEGs_ordered <- read.table("driver_DEGs_ordered.csv", quote="\"", comment.char="")

FeaturePlot(sub_dr, features = c("otx1", "vox", "vsx2", "pax6b", "dlc", "fabp7a", "stm", "msx1a", "rx1", "hmx4", "ascl1a", "wfdc1"), ncol = 6)



### EXPORT SEURAT OBJECTS TO PYTHON

library(Matrix)
####### DR15 TO PYTHON
### Extract seurat data to export to python for scvelo, BYPASS seurat version incompat with h5seurat
matrix <- GetAssayData(dr15_a2, assay = "SCT", slot = "data")
# Save as .mtx
writeMM(matrix, "dr15.mtx")
# Save gene names
write.csv(rownames(matrix), "dr15_genes.csv", row.names = FALSE)
# Extract clean metadata (keep only atomic columns)
meta_clean <- dr15_a2@meta.data[, sapply(dr15_a2@meta.data, function(x) is.numeric(x) || is.character(x) || is.factor(x))]
# Convert factors to characters
meta_clean <- data.frame(lapply(meta_clean, function(x) if (is.factor(x)) as.character(x) else x),
                         row.names = rownames(meta_clean))
write.csv(meta_clean, "dr15_meta.csv", row.names = TRUE)
#save PCA and UMAP embeddings
write.csv(Embeddings(dr15_a2, "pca"), "dr15_pca.csv", row.names = TRUE)
write.csv(Embeddings(dr15_a2, "umap"), "dr15_umap.csv", row.names = TRUE)

####### DR18 TO PYTHON
### Extract seurat data to export to python for scvelo, BYPASS seurat version incompat with h5seurat
matrix <- GetAssayData(dr18_a2, assay = "SCT", slot = "data")
# Save as .mtx
writeMM(matrix, "dr18.mtx")
# Save gene names
write.csv(rownames(matrix), "dr18_genes.csv", row.names = FALSE)
# Extract clean metadata (keep only atomic columns)
meta_clean <- dr18_a2@meta.data[, sapply(dr18_a2@meta.data, function(x) is.numeric(x) || is.character(x) || is.factor(x))]
# Convert factors to characters
meta_clean <- data.frame(lapply(meta_clean, function(x) if (is.factor(x)) as.character(x) else x),
                         row.names = rownames(meta_clean))
write.csv(meta_clean, "dr18_meta.csv", row.names = TRUE)
#save PCA and UMAP embeddings
write.csv(Embeddings(dr18_a2, "pca"), "dr18_pca.csv", row.names = TRUE)
write.csv(Embeddings(dr18_a2, "umap"), "dr18_umap.csv", row.names = TRUE)


####### DR23 TO PYTHON
### Extract seurat data to export to python for scvelo, BYPASS seurat version incompat with h5seurat
matrix <- GetAssayData(dr23_a2, assay = "SCT", slot = "data")
# Save as .mtx
writeMM(matrix, "dr23.mtx")
# Save gene names
write.csv(rownames(matrix), "dr23_genes.csv", row.names = FALSE)
# Extract clean metadata (keep only atomic columns)
meta_clean <- dr23_a2@meta.data[, sapply(dr23_a2@meta.data, function(x) is.numeric(x) || is.character(x) || is.factor(x))]
# Convert factors to characters
meta_clean <- data.frame(lapply(meta_clean, function(x) if (is.factor(x)) as.character(x) else x),
                         row.names = rownames(meta_clean))
write.csv(meta_clean, "dr23_meta.csv", row.names = TRUE)
#save PCA and UMAP embeddings
write.csv(Embeddings(dr23_a2, "pca"), "dr23_pca.csv", row.names = TRUE)
write.csv(Embeddings(dr23_a2, "umap"), "dr23_umap.csv", row.names = TRUE)



####### INTEGRATED TO PYTHON
# Extract the integrated assay and metadata
library(Matrix)

# Extract integrated assay
integrated_matrix <- GetAssayData(seurat_integrated_v4, assay = "integrated", slot = "data")

# Save as .mtx
writeMM(integrated_matrix, "integrated.mtx")

# Save gene names
write.csv(rownames(integrated_matrix), "integrated_genes.csv", row.names = FALSE)

# Extract clean metadata (keep only atomic columns)
meta_clean <- seurat_integrated_v4@meta.data[, sapply(seurat_integrated_v4@meta.data, function(x) is.numeric(x) || is.character(x) || is.factor(x))]

# Convert factors to characters
meta_clean <- data.frame(lapply(meta_clean, function(x) if (is.factor(x)) as.character(x) else x),
                         row.names = rownames(meta_clean))

write.csv(meta_clean, "meta.csv", row.names = TRUE)

# Optional: save PCA and UMAP embeddings
write.csv(Embeddings(seurat_integrated_v4, "pca"), "pca.csv", row.names = TRUE)
write.csv(Embeddings(seurat_integrated_v4, "umap"), "umap.csv", row.names = TRUE)


### Merge retinal cells

#Retina
# Subset por cluster en cada dataset
ret15 <- subset(dr15_a2, idents = c("Neural retina", "Proliferative region"))
ret18 <- subset(dr18_a2, idents = c("Neural retina", "Proliferative region"))
ret23 <- subset(dr23_a2, idents = c("Neural retina", "Proliferative region"))
ret_all <- merge(ret15, y = c(ret18, ret23), add.cell.ids = c("15hpf","18hpf","23hpf"))
ret_all <- SCTransform(ret_all, vars.to.regress = "percent.mt", verbose = FALSE)
ret_all <- RunPCA(ret_all)
ret_all <- RunUMAP(ret_all, dims = 1:30, verbose = FALSE)
ret_all <- FindNeighbors(ret_all, dims = 1:30, verbose = FALSE)
ret_all <- FindClusters(ret_all, verbose = FALSE, resolution = 2)
DimPlot(ret_all, label = TRUE, group.by = "orig.ident")
Idents(ret_all) <- ret_all@meta.data$orig.ident
ret_all <- PrepSCTFindMarkers(ret_all)

# Calculate gene markers
retina_time_markers <- FindAllMarkers(
  ret_all,
  logfc.threshold = 0.5,
  min.pct = 0.1,
  only.pos = TRUE
)
retina_time_markers <- FindAllMarkers(ret_all, logfc.threshold = 0.5, min.pct = 0.1)
retina_time_markers$pct.diff <- retina_time_markers$pct.1 - retina_time_markers$pct.2


# Integrate retinal cells
Idents(dr15_v4) <- dr15_v4@meta.data$cluster
Idents(dr18_v4) <- dr18_v4@meta.data$cluster
Idents(dr23_v4) <- dr23_v4@meta.data$cluster
retina15 <- subset(dr15_v4, idents = c("Neural retina", "Proliferative region"))
retina18 <- subset(dr18_v4, idents = c("Neural retina", "Proliferative region"))
retina23 <- subset(dr23_v4, idents = c("Neural retina", "Proliferative region"))

retina_list <- list(retina15, retina18, retina23)
features <- SelectIntegrationFeatures(object.list = retina_list, nfeatures = 3000)
retina_list <- PrepSCTIntegration(object.list = retina_list, anchor.features = features)
anchors <- FindIntegrationAnchors(object.list = retina_list, normalization.method = "SCT", anchor.features = features)
retina_integrated <- IntegrateData(anchorset = anchors, normalization.method = "SCT")
retina_integrated <- RunPCA(retina_integrated)
retina_integrated <- RunUMAP(retina_integrated, dims = 1:30)
retina_integrated <- FindNeighbors(retina_integrated, dims = 1:30)
retina_integrated <- FindClusters(retina_integrated, resolution = 0.5)
Idents(retina_integrated) <- retina_integrated@meta.data$orig.ident

DefaultAssay(retina_integrated) <- "SCT"

# Prepara el objeto para DE con SCT
retina_integrated <- PrepSCTFindMarkers(retina_integrated)

Idents(retina_integrated) <- retina_integrated@meta.data$orig.ident
# Ahora sí puedes buscar DEGs por timepoint
retina_integrated_time <- FindAllMarkers(
  retina_integrated,
  logfc.threshold = 0.5,
  min.pct = 0.1
)

Idents(retina_integrated) <- retina_integrated@meta.data$cluster
retina_integrated_clusters <- FindAllMarkers(retina_integrated, logfc.threshold = 0.5, min.pct = 0.1)

DimPlot(retina_integrated, label = TRUE, group.by = "cluster", cols = colors_clusters)

retina_integrated@meta.data$cluster_time <- paste(retina_integrated@meta.data$cluster, retina_integrated@meta.data$stage)
Idents(retina_integrated) <- retina_integrated@meta.data$cluster_time
retina_integrated_cluster_time <- FindAllMarkers(retina_integrated, logfc.threshold = 0.5, min.pct = 0.1)


####### retina integrated TO PYTHON
### Extract seurat data to export to python for scvelo, BYPASS seurat version incompat with h5seurat
matrix <- GetAssayData(retina_integrated, assay = "SCT", slot = "data")
# Save as .mtx
writeMM(matrix, "retina_integrated.mtx")
# Save gene names
write.csv(rownames(matrix), "retina_integrated_genes.csv", row.names = FALSE)
# Extract clean metadata (keep only atomic columns)
meta_clean <- retina_integrated@meta.data[, sapply(retina_integrated@meta.data, function(x) is.numeric(x) || is.character(x) || is.factor(x))]
# Convert factors to characters
meta_clean <- data.frame(lapply(meta_clean, function(x) if (is.factor(x)) as.character(x) else x),
                         row.names = rownames(meta_clean))
write.csv(meta_clean, "retina_integrated_meta.csv", row.names = TRUE)
#save PCA and UMAP embeddings
write.csv(Embeddings(retina_integrated, "pca"), "retina_integrated_pca.csv", row.names = TRUE)
write.csv(Embeddings(retina_integrated, "umap"), "retina_integrated_umap.csv", row.names = TRUE)


# functions_annotations.R

library(Seurat)
library(dplyr)

# --------------------------------------------------------------------------------
# Function: cell_annotation
# --------------------------------------------------------------------------------
# This function takes a Seurat object and an annotation dataset (zf_comb_annot)
# and assigns tags to cells according to log2FC and matches with annotations.
# --------------------------------------------------------------------------------
cell_annotation <- function(seurat_obj, log2FC = 0.5, annotation_data) {
  
  # Copy original object
  seurat_annot <- seurat_obj
  
  # Calculate markers for all clusters
  all_markers <- FindAllMarkers(seurat_annot, logfc.threshold = log2FC, min.pct = 0.1)
  
  # For each cluster, assign annotations based on genes from the dataset
  cluster_labels <- sapply(unique(all_markers$cluster), function(cl) {
    cluster_genes <- filter(all_markers, cluster == cl & avg_log2FC >= log2FC)$gene
    matched <- annotation_data[annotation_data$V5 %in% cluster_genes, "V2"]
    if(length(matched) == 0) return(NA)
    # Take the most frequent
    label <- names(sort(table(matched), decreasing = TRUE))[1]
    return(label)
  })
  
  # Save tags in meta.data
  seurat_annot@meta.data$SCT_cluster_annotation <- NA
  for(cl in names(cluster_labels)) {
    cells <- colnames(subset(seurat_annot, seurat_clusters == cl))
    seurat_annot@meta.data[cells, "SCT_cluster_annotation"] <- cluster_labels[cl]
  }
  
  return(list(seurat_annot, all_markers))
}

# --------------------------------------------------------------------------------
# Function: cluster_anatomical_enrichment
# --------------------------------------------------------------------------------
# This function perform a hypergeometric test for each cluster,
cluster_anatomical_enrichment <- function(seurat_object, annotation_data, log2FC = 0.5, exclude_terms = NULL) {
  
  # Filtrar anotaciones solo para genes detectados en Seurat
  genes_presentes <- rownames(seurat_object)
  annotation_data <- annotation_data %>% 
    filter(Genename %in% genes_presentes) %>% 
    distinct(Genename, Anatomical)
  
  # Definir background REAL
  background_genes <- unique(annotation_data$Genename)
  tfreq <- table(annotation_data$Anatomical)
  total_background <- length(background_genes)
  
  # Lista de resultados
  results_list <- list()
  
  for (cluster_id in levels(seurat_object)) {
    message("Procesando cluster: ", cluster_id)
    
    markers <- FindMarkers(seurat_object, ident.1 = cluster_id, logfc.threshold = log2FC, min.pct = 0.05)
    markers$pct_dif <- markers$pct.1 - markers$pct.2
    markers <- markers[markers$pct_dif >= 0.15, ]
    if (nrow(markers) < 8) next
    
    markers$genename <- rownames(markers)
    markers <- markers %>% 
      filter(avg_log2FC > 0.5) %>%
      mutate(corrected = avg_log2FC * pct.1) %>%
      filter(corrected > 0.5)
    
    gene_set <- unique(markers$genename)
    enriched_terms <- annotation_data %>% filter(Genename %in% gene_set)
    
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
    
    # Hipergeométrica con universo real
    df$pvalue <- mapply(function(k, K) {
      phyper(k - 1, K, total_background - K, length(gene_set), lower.tail = FALSE)
    }, k = df$count, K = tfreq[df$Anatomical])
    
    if (!is.null(exclude_terms)) {
      df <- df[!df$Anatomical %in% exclude_terms, ]
    }
    
    df$padj <- p.adjust(df$pvalue, method = "BH")
    
    genes_by_term <- enriched_terms %>% 
      group_by(Anatomical) %>%
      summarise(genes = paste(unique(Genename), collapse = "/"))

    
    df <- left_join(df, genes_by_term, by = "Anatomical")
    
    results_list[[cluster_id]] <- df[order(df$pvalue), ]
  }
  
  final_result <- do.call(rbind, results_list)
  rownames(final_result) <- NULL
  
  final_min_padj <- final_result %>%
    ungroup() %>%                 # por si acaso
    arrange(cluster, padj) %>%
    distinct(cluster, .keep_all = TRUE)
  
  annotation_df <- final_min_padj[, c("cluster", "Anatomical")]
  colnames(annotation_df) <- c("RNA_snn_res.5", "annotation")
  
  seurat_object@meta.data$annotation <-
    annotation_df$annotation[
      match(
        seurat_object@meta.data$RNA_snn_res.5,
        annotation_df$RNA_snn_res.5
      )
    ]
  
  
  return(final_result)
}


idents <- c("hindbrain", "optic vesicle", "Unknown", "diencephalon", "neural tube",
            "midbrain hindbrain boundary", "telencephalon", "muscle pioneer", "immature eye",
            "optic stalk", "neural crest telencephalon", "rhombomere 7", "rhombomere 3",
            "rhombomere 7", "central nervous system", "cranial ganglion",
            "alar plate midbrain region", "midbrain", "lens placode")


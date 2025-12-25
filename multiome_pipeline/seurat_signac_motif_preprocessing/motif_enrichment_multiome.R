# Loading libraries

library(Signac)
library(TFBSTools)
library(JASPAR2020)
library(RSQLite)
library(dplyr)
library(stringr)
library(BSgenome.Hsapiens.UCSC.hg38)

# Precomputed Seurat scMultiome dataset
#load("/data/scMultiome/paper/files/scMultiome.RData")

zebra_human_orthologs_zfin_gene_name <- read.delim("/data/GRN_software/dir_20_07_2023/orth_database/zebra_human_orthologs_zfin_gene_name.tsv", header=FALSE)

head(zebra_human_orthologs_zfin_gene_name)
#       V1      V2
#1    a1cf    A1CF
#2    aaas    AAAS
#3    aacs    AACS
#4   aadac   AADAC
#5 aadacl4 AADACL4
#6   aadat   AADAT


colnames(zebra_human_orthologs_zfin_gene_name) <- c("zebra_genename", "human_genename")


# Pipeline taken from: https://stuartlab.org/signac/articles/motif_vignette

#motif enrichment analysis

pfm <- getMatrixSet(
  x = JASPAR2020,
  opts = list(collection = "CORE", tax_group = 'vertebrates', all_versions = FALSE)
)

dr18_s <- AddMotifs(
  object = dr18_s,
  genome = fa,
  pfm = pfm
)


open.peaks <- AccessiblePeaks(dr18_s, idents = c("neural retina", "RPE", "lens placode", "optic stalk", "diencephalon", "telencephalon"))

get_domain_motifs <- function(
  domain_name,
  object,
  da_peaks,
  meta_feat,
  open_peaks,
  ortholog_table,
  deg_table,
  pval_cutoff = 0.05,
  log2fc_cutoff = 1,
  pct_ratio_cutoff = 1.5,
  top_n_motifs = 16
) {
  
  # 1. Filter fifferentially open peaks
  top.da.peaks <- da_peaks %>%
    filter(cluster == domain_name,
           p_val <= pval_cutoff,
           avg_log2FC >= log2fc_cutoff,
           pct.ratio >= pct_ratio_cutoff) %>%
    pull(gene)
  
  if(length(top.da.peaks) == 0){
    message(paste0("No differential open peaks for ", domain_name))
    return(NULL)
  }
  
  # 2. Calculate background peaks
  peaks.matched <- MatchRegionStats(
    meta.feature = meta_feat[open_peaks, ],
    query.feature = meta_feat[top.da.peaks, ],
    n = length(top.da.peaks) * 10
  )
  
  # 3. Search for enriched motifs
  enriched.motifs <- FindMotifs(object = object, features = top.da.peaks, background = peaks.matched)
  
  # 4. Clean names
  enriched.motifs <- enriched.motifs %>%
    mutate(
      motif_clean = str_replace(motif.name, "\\(var\\.\\d+\\)", ""),
      motif_clean = str_split(motif_clean, "::", simplify = TRUE)[, 1]
    ) %>%
    filter(!str_detect(motif.name, "::"))
  
  # 5. Overexpressed genes
  de_genes <- deg_table %>%
    filter(cluster == domain_name, avg_log2FC >= 1, p_val <= 0.05) %>%
    pull(gene)
  
  # 6. Intersect with orthologs
  motifs_with_orthologs <- enriched.motifs %>%
    left_join(ortholog_table, by = c("motif_clean" = "human_genename")) %>%
    filter(!is.na(zebra_genename), zebra_genename %in% de_genes)
  
  # 7. Select top motifs
  motifs_to_plot <- head(unique((motifs_with_orthologs %>% filter(pvalue <= 0.05))$motif.name), n = 20)
  
  if(length(motifs_to_plot) == 0){
    message(paste0("Significant motifs not found for ", domain_name))
    return(list(
      motifs_with_orthologs = motifs_with_orthologs,
      motifs = character(0),
      plot = NULL
    ))
  }
  
  # 8. Generate plot
  message(paste0("Generating MotifPlot for ", domain_name, " (", length(motifs_to_plot), " available motifs)..."))
  p <- MotifPlot(object = object, motifs = head(motifs_to_plot, n = top_n_motifs))
  
  return(list(
    motifs_with_orthologs = motifs_with_orthologs,
    motifs = motifs_to_plot,
    plot = p
  ))
}

domains <- c("neural retina", "RPE", "optic stalk", "lens placode", "diencephalon", "telencephalon")

results <- lapply(domains, function(dom) {
  get_domain_motifs(
    domain_name = dom,
    object = dr18_s,
    da_peaks = all_DOCRs_domains,
    meta_feat = meta.feature,
    open_peaks = open.peaks,
    ortholog_table = zebra_human_orthologs_zfin_gene_name,
    deg_table = all_DEGs_domains,
    pval_cutoff = 0.05,
    log2fc_cutoff = 1,
    pct_ratio_cutoff = 1.5,
    top_n_motifs = 16
  )
})

names(results) <- domains


## Create venn diagrams for more than 3 datasets


gene_list <- list( Retina = results[[1]][[2]], RPE = results[[2]][[2]], 
		Optic = results[[3]][[2]], Lens = results[[4]][[2]], 
		Diencephalon = results[[5]][[2]], Telencephalon = results[[6]][[2]])
		
		
all_genes <- unique(unlist(gene_list))
df <- tibble(Gene = all_genes)
for (tissue in names(gene_list)) { df[[tissue]] <- df$Gene %in% gene_list[[tissue]] }
upset( df, intersect = names(gene_list))


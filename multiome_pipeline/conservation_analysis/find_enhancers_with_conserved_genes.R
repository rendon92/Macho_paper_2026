

library(dplyr)
library(tidyr)
library(GenomicRanges)
library(stringr)

retina_hits_DOCRs_vs_human_len50 <- read.delim("/data/scMultiome/paper/files/enhancer_sequence_identity/retina_hits_DOCRs_vs_human_len50.tsv", header=FALSE)
zebra_genes <- read.delim("/data/scMultiome/paper/files/enhancer_sequence_identity/zebra_genes.bed", header=FALSE)
human_genes <- read.delim("/data/scMultiome/paper/files/enhancer_sequence_identity/human_genes.bed", header=FALSE)

 head(retina_hits_DOCRs_vs_human_len50)
#                    V1 V2     V3 V4 V5 V6  V7  V8        V9       V10      V11  V12  V13       V14 V15
#1 17:31664934-31666128 14 91.228 57  5  0 642 698  74214314  74214258 1.24e-12 81.3 1194 107043718   -
#2 22:32432165-32433083 10 88.136 59  7  0 851 909  27002431  27002373 1.93e-10 73.6  918 133797422   -
#3 22:32432165-32433083 12 82.432 74 13  0 838 911  50522527  50522600 1.04e-08 67.8  918 133275309   +
#4 22:32432165-32433083  1 85.246 61  9  0 851 911 244519530 244519470 3.94e-08 65.9  918 248956422   -
#5 22:32432165-32433083  8 82.432 74 11  2 838 910  99899645  99899573 2.13e-06 60.2  918 145138636   -
#6  9:39549841-39550769 14 86.538 52  7  0 877 928  74868377  74868326 2.15e-06 60.2  928 107043718   -

head(zebra_genes)
#  V1       V2       V3         V4
#1  4  1722898  1730920   fgfr1op2
#2  4  8213062  8515824       ERC1
#3  4 30454495 30460765 zgc:173716
#4  4  1732824  1757460     TM7SF3
#5  4 40880432 40883162 zgc:174650
#6  4 19700307 19704628       pax4


head(human_genes)
#  V1       V2       V3     V4
#1  1  3069167  3438621 PRDM16
#2  1  2403963  2413797  PEX10
#3  1 10472287 10630758  PEX14
#4  1  2425979  2505532  PLCH2
#5  1  9292893  9369532  SPSB1
#6  1  9035105  9088478 SLC2A5




annotate_hits_with_genes <- function(hits, zebra_genes, human_genes,
                                     flank_zebra = 100000, flank_human = 200000) {
    # Rename gene columns
    colnames(zebra_genes) <- c("zebra_chr", "zebra_start", "zebra_end", "gene_name")
    colnames(human_genes) <- c("human_chr", "human_start", "human_end", "gene_name")
    
    # ZEBRAFISH HITS
    zebra_hits <- hits %>%
        separate(V1, into = c("zebra_chr", "coords"), sep = ":") %>%
        separate(coords, into = c("zebra_start", "zebra_end"), sep = "-") %>%
        mutate(zebra_start = as.numeric(zebra_start), zebra_end = as.numeric(zebra_end)) %>%
        dplyr::select(zebra_chr, zebra_start, zebra_end)
    
    zebra_hits_expanded <- zebra_hits %>%
        mutate(zebra_start = pmax(0, zebra_start - flank_zebra),
               zebra_end   = zebra_end + flank_zebra)
    
    # HUMAN HITS
    human_hits <- hits %>%
        transmute(human_chr   = V2,
                  human_start = V9,
                  human_end   = V10,
                  human_strand = V15)
    
    human_hits_expanded <- human_hits %>%
        mutate(human_start = pmax(0, human_start - flank_human),
               human_end   = human_end + flank_human)
    
    # GenomicRanges
    zebra_gr <- GRanges(seqnames = zebra_genes$zebra_chr,
                        ranges   = IRanges(zebra_genes$zebra_start, zebra_genes$zebra_end),
                        gene     = zebra_genes$gene_name)
    
    human_gr <- GRanges(seqnames = human_genes$human_chr,
                        ranges   = IRanges(human_genes$human_start, human_genes$human_end),
                        gene     = human_genes$gene_name)
    
    hits_zebra_gr <- GRanges(seqnames = zebra_hits_expanded$zebra_chr,
                             ranges   = IRanges(zebra_hits_expanded$zebra_start,
                                                zebra_hits_expanded$zebra_end))
    
    hits_human_gr <- GRanges(seqnames = human_hits_expanded$human_chr,
                             ranges   = IRanges(human_hits_expanded$human_start,
                                                human_hits_expanded$human_end))
    
    # Overlaps
    zebra_ov <- findOverlaps(hits_zebra_gr, zebra_gr)
    human_ov <- findOverlaps(hits_human_gr, human_gr)
    
    zebra_hits_with_genes <- data.frame(
        hit_id = queryHits(zebra_ov),
        gene   = zebra_gr$gene[subjectHits(zebra_ov)]
    )
    
    human_hits_with_genes <- data.frame(
        hit_id = queryHits(human_ov),
        gene   = human_gr$gene[subjectHits(human_ov)]
    )
    
    # Colapsar genes por hit
    zebra_hits_genes <- zebra_hits_with_genes %>%
        group_by(hit_id) %>%
        summarise(genes = paste(unique(gene), collapse = ","), .groups = "drop")
    
    human_hits_genes <- human_hits_with_genes %>%
        group_by(hit_id) %>%
        summarise(genes = paste(unique(gene), collapse = ","), .groups = "drop")
    
    # ---- Añadir a tabla original ----
    hits$zebra_genes <- zebra_hits_genes$genes[match(1:nrow(hits), zebra_hits_genes$hit_id)]
    hits$human_genes <- human_hits_genes$genes[match(1:nrow(hits), human_hits_genes$hit_id)]
    
    return(hits)
}





annotated_hits <- annotate_hits_with_genes(retina_hits_DOCRs_vs_human_len50,
                                           zebra_genes,
                                           human_genes)





# Retain hits with at least one ortholog gene nearby in zebrafish and human (100kb of distance in zebrafish, 200kb in human) 

orthologs <- read.delim("/data/GRN_software/dir_20_07_2023/orth_database/zebra_human_orthologs_zfin_gene_name.tsv", header=FALSE, )


colnames(orthologs) <- c("zebra", "human")


orthologs$zebra <- toupper(orthologs$zebra)
orthologs$human <- toupper(orthologs$human)

# function to match orthologs
annotated_hits$ortholog_match <- apply(annotated_hits, 1, function(row) {
  zebra <- unlist(strsplit(row[["zebra_genes"]], ","))
  human <- unlist(strsplit(row[["human_genes"]], ","))
  
  zebra <- toupper(trimws(zebra))
  human <- toupper(trimws(human))
  
  # Search fo human orthologs for each zebrafish gene
  zebra2human <- orthologs$human[match(zebra, orthologs$zebra)]
  
  # Intersect human genes and zebrafish orthologs
  common <- intersect(human, zebra2human)
  
  if (length(common) > 0) paste(common, collapse = ",") else NA
})
annotated_hits$ortholog_match[annotated_hits$ortholog_match == "NA"] <- NA

# List of orthologs near the hits (for both species)
orthologs_char <- as.character(annotated_hits$ortholog_match)

# Separate by commas
all_orthologs <- unlist(strsplit(orthologs_char, ","))

# Remove whitspaces around the gene names
all_orthologs <- trimws(all_orthologs)

# Remove duplicated and empty entries
all_orthologs <- unique(all_orthologs[all_orthologs != ""])

all_orthologs <- sort(all_orthologs)



##Subset of those genes which have also been reported in MACs
MACs_genes <- read.table("/data/scMultiome/paper/files/enhancer_sequence_identity/MACs_genes.tsv", quote="\"", comment.char="")

head(MACs_genes)
#       V1
#1   ABCB6
#2    ACTB
#3   ACTG1
#4 ALDH1A3
#5   ATOH7
#6    BCOR

all_orthologs_up <- toupper(all_orthologs)
MACs_genes_up <- toupper(MACs_genes$V1)

# Intersection
common_genes <- intersect(all_orthologs_up, MACs_genes_up)

  
# load zfin genes and descriptions  
zfin_genes <- read.delim("/data/scMultiome/paper/files/enhancer_sequence_identity/zfin_genes.gff3", header=FALSE)  
  
tf_keywords <- c("transcription factor", "homeobox", "basic helix-loop-helix",
                 "forkhead", "paired box", "zinc finger", "ETS")

# Filer and extract
tf_rows <- zfin_genes[grepl(paste(tf_keywords, collapse="|"), zfin_genes$V9, ignore.case = TRUE), ]
tf_genes <- str_extract(tf_rows$V9, "Name=([A-Za-z0-9]+)")
tf_genes <- gsub("Name=", "", tf_genes)
tf_genes <- unique(tf_genes)



# Transform into character
orthologs$zebra <- as.character(orthologs$zebra)
orthologs$human <- as.character(orthologs$human)
tf_genes <- as.character(tf_genes)
all_orthologs_up <- as.character(all_orthologs_up)

# Filter TFs in the list of conserved hits
tf_zebra_in_hits <- intersect(tf_genes, tolower(orthologs$zebra))  
tf_human_in_hits <- orthologs$human[tolower(orthologs$zebra) %in% tf_zebra_in_hits]

# Intersect with the universe of orthlogs near the hits
tf_human_in_hits <- intersect(tf_human_in_hits, all_orthologs_up)
all_orthologs_up <- as.character(all_orthologs_up)
tf_human_in_hits <- as.character(tf_human_in_hits)
MACs_genes_up <- as.character(MACs_genes_up)

# Create table
orthologs_table <- data.frame(
  gene = all_orthologs_up,
  TF = all_orthologs_up %in% tf_human_in_hits,
  MACs = all_orthologs_up %in% MACs_genes_up
)

# Order (optional): first MACs and then TF
orthologs_table <- orthologs_table[order(orthologs_table$MACs, orthologs_table$TF, decreasing = TRUE), ]

#
colnames(annotated_hits) <- c("Zebrafish_Query_Coord", "Human_chromosome", "%ident", "align_len", "mismatch", "gapopen", "Query_start", "Query_end", "Human_start", "Human_end", "evalue", "bitscore", "Query_len", "Human_len", "strand", "Zebrafish_genes", "Human_genes", "ortholog_match")
annotated_hits_filtered <- annotated_hits %>%
  filter(!is.na(ortholog_match))
annotated_hits_filtered <- annotated_hits_filtered[, c(
  "Human_chromosome",  
  "Human_start",       
  "Human_end",        
  "strand",            
  "Zebrafish_Query_Coord",
  "Human_len", "Query_len", "align_len", "%ident", "mismatch", "gapopen",
  "evalue", "bitscore", "Zebrafish_genes", "Human_genes", "ortholog_match"
)]

human_assembly <- "GRCh38_Ensembl"
zebra_assembly <- "GRCz11_NCBI"


annotated_hits_filtered$Human_assembly <- human_assembly
annotated_hits_filtered$Zebrafish_assembly <- zebra_assembly

# Create table of occurrence (TF positive? MAC positive?)

library(dplyr)
library(ggplot2)
library(tidyr)

df <- data.frame(gene = all_orthologs) %>%
  mutate(
    MACs = gene %in% MACs_genes_up,
    TF = gene %in% tf_human_in_hits
  ) %>%
  pivot_longer(cols = c(MACs, TF), names_to = "Attribute", values_to = "Present")

ggplot(df, aes(x = Attribute, y = reorder(gene, gene), fill = Present)) +
  geom_tile(color = "white", width = 0.5, height = 0.8) +  # ajustar altura de filas
  scale_fill_manual(values = c("TRUE" = "red", "FALSE" = "grey90")) +
  theme_minimal() +
  theme(
    axis.title = element_blank(),
    axis.text.x = element_text(size = 10),
    axis.text.y = element_text(size = 4),   # aún más pequeño
    panel.grid = element_blank()
  ) +
  labs(
    title = "Ortholog genes near hits",
    subtitle = "Red: present, grey: absent"
  )

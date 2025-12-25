#!/bin/bash

blastn \
  -task blastn \
  -query enhancers_retina.fa \
  -db human_genome \
  -evalue 1e-1 \
  -word_size 9 \
  -reward 1 -penalty -2 -gapopen 2 -gapextend 2 \
  -soft_masking true -dust yes \
  -max_target_seqs 20 -culling_limit 4 -max_hsps 10 \
  -strand both \
  -num_threads 8 \
  -outfmt "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen" \
  -out retina_hits_DOCRs_vs_human.tsv



#awk '$4>=100' retina_hits_DOCRs_vs_human.tsv | awk -v OFS='\t' '{strand = ($9 < $10 ? "+" : "-"); print $0, strand}' > retina_hits_DOCRs_vs_human_len100.tsv
awk '$4>=50' retina_hits_DOCRs_vs_human.tsv | awk -v OFS='\t' '{strand = ($9 < $10 ? "+" : "-"); print $0, strand}' > retina_hits_DOCRs_vs_human_len50.tsv

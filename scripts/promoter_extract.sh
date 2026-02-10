#!/usr/bin/env bash
set -euo pipefail

# args:
# 1 gene.bed (chr start end id score strand) 0-based
# 2 chr.length (chr len)
# 3 genome.fa
# 4 promoter_len
# 5 out.bed
# 6 out.fa

GENE_BED="$1"
CHRLEN="$2"
GENOME="$3"
L="$4"
OUTBED="$5"
OUTFA="$6"

# need bedtools
command -v bedtools >/dev/null 2>&1 || { echo "ERROR: bedtools not found"; exit 1; }

awk -v L="$L" 'BEGIN{FS=OFS="\t"}
  NR==FNR{len[$1]=$2; next}
  {
    chr=$1; start=$2; end=$3; id=$4; score=$5; strand=$6;
    if(!(chr in len)) next;

    if(strand=="+"){
      pstart=start-L; if(pstart<0) pstart=0;
      pend=start;
    } else if(strand=="-"){
      pstart=end;
      pend=end+L; if(pend>len[chr]) pend=len[chr];
    } else {
      next
    }

    if(pend>pstart){
      print chr, pstart, pend, id, score, strand
    }
  }' "$CHRLEN" "$GENE_BED" > "$OUTBED"

# -name 会用第4列作为 fasta header；-s 根据链取反向互补
bedtools getfasta -fi "$GENOME" -bed "$OUTBED" -s -name > "$OUTFA"

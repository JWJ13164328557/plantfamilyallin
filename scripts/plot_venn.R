args <- commandArgs(trailingOnly = TRUE)
getArg <- function(flag) {
  i <- match(flag, args)
  if (is.na(i) || i == length(args)) stop(paste("Missing", flag))
  args[i + 1]
}
venn_tsv <- getArg("--venn_tsv")
family <- getArg("--family")
out <- getArg("--out")

suppressPackageStartupMessages({
  library(readr)
  library(ggVennDiagram)
  library(ggplot2)   # <--- 关键：ggtitle() 在这里
})

df <- read_tsv(venn_tsv, show_col_types = FALSE)
A <- df$gene[df$BLAST == 1]
B <- df$gene[df$PFAM == 1]
C <- df$gene[df$HMM  == 1]

lst <- list(BLAST = A, PFAM = B, HMM = C)
p <- ggVennDiagram(lst, label_alpha = 0) +
  ggtitle(paste0(family, " Family Identification (Venn)"))

ggsave(filename = out, plot = p, width = 7, height = 6)

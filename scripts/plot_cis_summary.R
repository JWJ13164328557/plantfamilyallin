args <- commandArgs(trailingOnly = TRUE)
getArg <- function(flag) {
  i <- match(flag, args)
  if (is.na(i) || i == length(args)) stop(paste("Missing", flag))
  args[i + 1]
}
in_tsv <- getArg("--in_tsv")
family <- getArg("--family")
out <- getArg("--out")

suppressPackageStartupMessages({
  library(readr); library(ggplot2)
})
df <- read_tsv(in_tsv, show_col_types = FALSE)
if (nrow(df) == 0) {
  pdf(out); plot.new(); text(0.5,0.5,"No FIMO hits"); dev.off()
} else {
  p <- ggplot(df, aes(x=reorder(motif_id, count), y=count)) +
    geom_col() + coord_flip() + theme_bw() +
    labs(title=paste0(family, " Cis-elements (FIMO)"), x="Motif", y="Count")
  ggsave(out, p, width=7, height=6)
}

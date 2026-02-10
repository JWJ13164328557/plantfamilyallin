#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(ggplot2)
})

opt_list <- list(
  make_option("--family", type="character"),
  make_option("--syntenic", type="character"),
  make_option("--out", type="character")
)
opt <- parse_args(OptionParser(option_list=opt_list))

fam <- fread(opt$family)
syn <- fread(opt$syntenic)

fam[, Ks := as.numeric(Ks)]
syn[, Ks := as.numeric(Ks)]
fam <- fam[is.finite(Ks)]
syn <- syn[is.finite(Ks)]

fam_median <- if (nrow(fam) > 0) median(fam$Ks) else NA_real_

p <- ggplot() +
  geom_density(data=syn, aes(x=Ks), linewidth=0.9, alpha=0.35, fill="grey70", color="grey30") +
  geom_density(data=fam, aes(x=Ks), linewidth=1.1, color="#1b9e77") +
  theme_bw(base_size=12) +
  labs(x="Ks", y="Density",
       title="Ks distribution: syntenic background vs family paralogs") +
  theme(
    plot.title=element_text(hjust=0.5, face="bold"),
    panel.grid.minor=element_blank()
  )

if (is.finite(fam_median)) {
  p <- p + geom_vline(xintercept=fam_median, linetype="dashed", linewidth=0.8, color="#1b9e77") +
    annotate("text", x=fam_median, y=Inf, vjust=1.5,
             label=paste0("family median=", sprintf("%.3f", fam_median)),
             color="#1b9e77", size=3.4)
}

ggsave(opt$out, p, width=8.2, height=5.2)

#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(ggplot2)
})

opt_list <- list(
  make_option("--family_kaks", type="character"),
  make_option("--syntenic_kaks", type="character"),
  make_option("--out_pdf", type="character"),
  make_option("--family_name", type="character", default="Family"),
  make_option("--xmax", type="double", default=5)
)
opt <- parse_args(OptionParser(option_list=opt_list))

fam <- fread(opt$family_kaks)
syn <- fread(opt$syntenic_kaks)

fam[, Ks := as.numeric(Ks)]
syn[, Ks := as.numeric(Ks)]
fam <- fam[is.finite(Ks)]
syn <- syn[is.finite(Ks)]

fam[, group := opt$family_name]
syn[, group := "Genome syntenic anchors"]

dt <- rbindlist(list(syn[, .(Ks, group)], fam[, .(Ks, group)]), use.names=TRUE)

p <- ggplot(dt, aes(x=Ks, color=group, fill=group)) +
  geom_density(alpha=0.20, linewidth=1.0, adjust=1.1) +
  coord_cartesian(xlim=c(0, opt$xmax)) +
  theme_bw(base_size = 12) +
  labs(x="Ks", y="Density", title="Ks distribution: family paralogs vs genome-wide syntenic anchors") +
  theme(
    plot.title = element_text(hjust=0.5, face="bold"),
    legend.title = element_blank(),
    panel.grid.minor = element_blank()
  )

ggsave(opt$out_pdf, p, width=8.5, height=5.3)

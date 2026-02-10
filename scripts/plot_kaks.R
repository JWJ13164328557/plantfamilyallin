#!/usr/bin/env Rscript
suppressPackageStartupMessages({
  library(optparse)
  library(ggplot2)
  library(data.table)
  library(RColorBrewer)
})

opt_list <- list(
  make_option(c("--kaks_tsv"), type="character"),
  make_option(c("--out_ks"), type="character"),
  make_option(c("--out_w"), type="character"),
  make_option(c("--out_scatter"), type="character")
)
opt <- parse_args(OptionParser(option_list=opt_list))
dt <- fread(opt$kaks_tsv)

dt[, Ks := as.numeric(Ks)]
dt[, Ka := as.numeric(Ka)]
dt[, KaKs := as.numeric(KaKs)]
dt <- dt[is.finite(Ks) & is.finite(Ka) & is.finite(KaKs)]

# 设置调色板
color_palette <- brewer.pal(3, "Set2")

# 1) Ks 分布图优化
p1 <- ggplot(dt, aes(x=Ks, fill=type)) +
  geom_histogram(bins=60, alpha=0.7, position="identity", color="black") +
  scale_fill_manual(values=color_palette) +
  theme_bw(base_size = 14) +
  labs(x="Ks", y="Count", title="Ks Distribution") +
  theme(plot.title=element_text(hjust=0.5),
        axis.title=element_text(size=12),
        axis.text=element_text(size=10),
        legend.title=element_text(size=12),
        legend.text=element_text(size=10)) +
  theme(panel.grid.major = element_line(size = 0.5, linetype = 'dotted', color = 'gray')) +
  theme(panel.grid.minor = element_blank())
ggsave(opt$out_ks, p1, width=8.5, height=5.5)

# 2) Ka/Ks 分布图优化
p2 <- ggplot(dt, aes(x=type, y=KaKs, fill=type)) +
  geom_violin(trim=TRUE, alpha=0.8, scale="area") +
  geom_boxplot(width=0.15, outlier.size=0.3, alpha=0.5, color="black") +
  scale_fill_manual(values=color_palette) +
  theme_bw(base_size = 14) +
  labs(x="", y="Ka/Ks (ω)", title="Ka/Ks Distribution") +
  theme(plot.title=element_text(hjust=0.5),
        axis.title=element_text(size=12),
        axis.text.x=element_text(angle=25, hjust=1, size=10),
        axis.text.y=element_text(size=10),
        legend.title=element_text(size=12),
        legend.text=element_text(size=10)) +
  theme(panel.grid.major = element_line(size = 0.5, linetype = 'dotted', color = 'gray')) +
  theme(panel.grid.minor = element_blank())
ggsave(opt$out_w, p2, width=8.5, height=5.5)

# 3) Ka vs Ks 散点图优化
p3 <- ggplot(dt, aes(x=Ks, y=Ka, color=type)) +
  geom_point(alpha=0.7, size=2) +
  scale_color_manual(values=color_palette) +
  theme_bw(base_size = 14) +
  labs(x="Ks", y="Ka", title="Ka vs Ks") +
  theme(plot.title=element_text(hjust=0.5),
        axis.title=element_text(size=12),
        axis.text=element_text(size=10),
        legend.title=element_text(size=12),
        legend.text=element_text(size=10)) +
  theme(panel.grid.major = element_line(size = 0.5, linetype = 'dotted', color = 'gray')) +
  theme(panel.grid.minor = element_blank())
ggsave(opt$out_scatter, p3, width=7.5, height=5.8)

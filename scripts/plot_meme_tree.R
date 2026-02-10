#!/usr/bin/env Rscript

# ---------------------------
# args
# ---------------------------
args <- commandArgs(trailingOnly = TRUE)

getArg <- function(flag) {
  i <- match(flag, args)
  if (is.na(i) || i == length(args)) stop(paste("Missing", flag))
  args[i + 1]
}
hasArg <- function(flag) flag %in% args

treefile   <- getArg("--tree")
motif_tsv  <- getArg("--motif_tsv")
len_tsv    <- getArg("--len_tsv")
out        <- getArg("--out")
family     <- getArg("--family")

domain_tsv <- if (hasArg("--domain_tsv")) getArg("--domain_tsv") else ""
gene_tsv   <- if (hasArg("--gene_tsv"))   getArg("--gene_tsv")   else ""

suppressPackageStartupMessages({
  library(ape)
  library(ggtree)
  library(ggplot2)
  library(dplyr)
  library(readr)
  library(patchwork)
})

# ---------------------------
# helpers
# ---------------------------
read_tsv_maybe_header <- function(path, colnames_vec) {
  if (!file.exists(path) || file.info(path)$size == 0) {
    return(tibble())
  }
  first <- readLines(path, n = 1, warn = FALSE)
  # 有表头：第一行出现 seq_id/feature/domain/start/end 之一即可认为有表头
  has_header <- any(grepl("\\bseq_id\\b|\\bdomain\\b|\\bfeature\\b|\\bstart\\b|\\bend\\b", first, ignore.case = TRUE))
  if (has_header) {
    read_tsv(path, show_col_types = FALSE, comment = "#")
  } else {
    read_tsv(path, col_names = colnames_vec, show_col_types = FALSE, comment = "#")
  }
}

stable_levels_numeric_suffix <- function(x) {
  # x like Motif1 Motif2 ... or MEME-1 ...
  tibble(v = unique(x)) |>
    mutate(idx = suppressWarnings(as.integer(gsub("[^0-9]", "", v)))) |>
    arrange(is.na(idx), idx, v) |>
    pull(v)
}

# ---------------------------
# read tree + tip y
# ---------------------------
tr <- read.tree(treefile) |> ladderize(right = TRUE)

p_tree0 <- ggtree(tr, size = 0.35) +
  geom_tiplab(size = 2, align = FALSE, offset = 0.01) +
  theme_tree2()

tip_pos <- p_tree0$data |>
  filter(isTip) |>
  select(label, y)

y_rng <- range(tip_pos$y, na.rm = TRUE)
# ---------------------------
# adaptive tip label offset (based on tree x-range)
# ---------------------------
x_rng <- range(p_tree0$data$x, na.rm = TRUE)
x_off <- diff(x_rng) * 0.01   # 取树宽的 1% 作为 label offset
if (!is.finite(x_off) || x_off <= 0) x_off <- 0.005  # 兜底

# ---------------------------
# read length
# ---------------------------
lens_raw <- read_tsv(len_tsv, show_col_types = FALSE)
stopifnot(all(c("seq_id", "length") %in% colnames(lens_raw)))

lens <- lens_raw |>
  transmute(
    label    = as.character(seq_id),
    prot_len = as.numeric(length)
  ) |>
  inner_join(tip_pos, by = "label")

xmax_len <- max(lens$prot_len, na.rm = TRUE)

# ---------------------------
# read motifs
# ---------------------------
motif_raw <- read_tsv(motif_tsv, show_col_types = FALSE)
need <- c("seq_id", "motif", "start", "end")
stopifnot(all(need %in% colnames(motif_raw)))

motif <- motif_raw |>
  transmute(
    label = as.character(seq_id),
    motif = as.character(motif),
    start = suppressWarnings(as.numeric(start)),
    end   = suppressWarnings(as.numeric(end))
  ) |>
  inner_join(tip_pos, by = "label") |>
  mutate(
    xmin = pmin(start, end, na.rm = TRUE),
    xmax = pmax(start, end, na.rm = TRUE)
  ) |>
  filter(!is.na(xmin), !is.na(xmax))

# motif palette stable
pal_motif <- NULL
if (nrow(motif) > 0) {
  motif_levels <- stable_levels_numeric_suffix(motif$motif)
  motif$motif <- factor(motif$motif, levels = motif_levels)
  # ---- publication-quality motif palette (10, colorblind-friendly) ----
motif_base10 <- c(
  "#E69F00", # orange
  "#56B4E9", # sky blue
  "#009E73", # bluish green
  "#F0E442", # yellow
  "#0072B2", # blue
  "#D55E00", # vermillion
  "#CC79A7", # reddish purple
  "#999999", # grey
  "#332288", # indigo (Tol-like extension)
  "#88CCEE"  # light cyan (Tol-like extension)
)

# 如果 motif_levels > 10，就自动补色（很少见，但做个兜底）
if (length(motif_levels) <= 10) {
  pal_motif <- setNames(motif_base10[seq_along(motif_levels)], motif_levels)
} else {
  pal_motif <- setNames(
    c(motif_base10, grDevices::hcl.colors(length(motif_levels) - 10, "Okabe-Ito")),
    motif_levels
  )
}

}

# ---------------------------
# read domain.tsv (protein aa)
# ---------------------------
domain_df <- tibble()
pal_domain <- NULL

if (nchar(domain_tsv) > 0 && file.exists(domain_tsv) && file.info(domain_tsv)$size > 0) {
  dom_raw <- read_tsv_maybe_header(domain_tsv, c("seq_id", "domain", "start", "end"))
  # 兜底：有些文件可能列数>4
  if (!all(c("seq_id", "domain", "start", "end") %in% colnames(dom_raw))) {
    stop("domain.tsv must have 4 columns: seq_id domain start end (header optional)")
  }

  domain_df <- dom_raw |>
    transmute(
      label  = as.character(seq_id),
      domain = as.character(domain),
      start  = suppressWarnings(as.numeric(start)),
      end    = suppressWarnings(as.numeric(end))
    ) |>
    inner_join(tip_pos, by = "label") |>
    mutate(
      xmin = pmin(start, end, na.rm = TRUE),
      xmax = pmax(start, end, na.rm = TRUE)
    ) |>
    filter(!is.na(xmin), !is.na(xmax))

  if (nrow(domain_df) > 0) {
    dom_levels <- stable_levels_numeric_suffix(domain_df$domain)
    domain_df$domain <- factor(domain_df$domain, levels = dom_levels)
    # ---- publication-quality domain palette (Tol muted style) ----
tol_muted <- c(
  "#332288", "#88CCEE", "#44AA99", "#117733",
  "#999933", "#DDCC77", "#CC6677", "#882255", "#AA4499"
)

# 若 domain 很多就循环使用（一般 Pfam domain 不会特别多）
pal_domain <- setNames(tol_muted[ (seq_along(dom_levels) - 1) %% length(tol_muted) + 1 ], dom_levels)

  }
}

# ---------------------------
# read gene_structure.tsv (bp)
# + convert to per-gene relative coordinates
# + build intron connectors from exon blocks
# ---------------------------
gene_df <- tibble()
intron_df <- tibble()
pal_gene <- NULL

if (nchar(gene_tsv) > 0 && file.exists(gene_tsv) && file.info(gene_tsv)$size > 0) {
  g_raw <- read_tsv_maybe_header(gene_tsv, c("seq_id", "feature", "start", "end"))
  if (!all(c("seq_id", "feature", "start", "end") %in% colnames(g_raw))) {
    stop("gene_structure.tsv must have 4 columns: seq_id feature start end (header optional)")
  }

  gene_df <- g_raw |>
    transmute(
      label   = as.character(seq_id),
      feature = as.character(feature),
      start   = suppressWarnings(as.numeric(start)),
      end     = suppressWarnings(as.numeric(end))
    ) |>
    inner_join(tip_pos, by = "label") |>
    mutate(
      xmin = pmin(start, end, na.rm = TRUE),
      xmax = pmax(start, end, na.rm = TRUE),
      feature = case_when(
        grepl("utr", tolower(feature)) ~ "UTR",
        tolower(feature) == "cds"      ~ "CDS",
        tolower(feature) == "exon"     ~ "exon",
        TRUE ~ feature
      )
    ) |>
    filter(!is.na(xmin), !is.na(xmax)) |>
    group_by(label) |>
    mutate(
      g0   = min(xmin, na.rm = TRUE),
      xmin = xmin - g0,
      xmax = xmax - g0
    ) |>
    ungroup()

  # intron connectors: use exon blocks only
  exon_df <- gene_df |>
    filter(feature == "exon") |>
    group_by(label, y) |>
    arrange(xmin, .by_group = TRUE) |>
    summarise(
      x1 = head(xmax, -1),
      x2 = tail(xmin, -1),
      .groups = "drop"
    ) |>
    filter(length(x1) > 0) |>
    tidyr::unnest(cols = c(x1, x2))

  intron_df <- exon_df

  # gene palette (固定颜色 + 有图例)
  pal_gene <- c(
  UTR  = "#D9D9D9",
  exon = "#A6CEE3",
  CDS  = "#1F78B4"
)


  gene_df$feature <- factor(gene_df$feature, levels = c("UTR", "exon", "CDS"))
}

message(
  "tips=", length(tr$tip.label),
  " lens=", nrow(lens),
  " motif=", nrow(motif),
  " domain_n=", nrow(domain_df),
  " gene_n=", nrow(gene_df)
)

# ---------------------------
# Panel1: tree
# ---------------------------
p_tree <- ggtree(tr, size = 0.35) +
  geom_tiplab(size = 2, align = FALSE, offset = x_off* 0.02) +
  ggtitle("Phylogenetic Tree") +
  theme_tree2() +
  coord_cartesian(ylim = y_rng, clip = "off") +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title = element_blank(),
    plot.margin = margin(5, 0, 5, 2)
  )

# ---------------------------
# Panel2: domain track (protein aa)
# ---------------------------
p_domain <- ggplot() +
  ggtitle("Domain Analysis") +
  scale_y_continuous(limits = y_rng, expand = c(0, 0)) +
  theme_minimal(base_size = 10) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = margin(5, 5, 5, 5)
  )

if (nrow(domain_df) > 0) {
  p_domain <- p_domain +
    geom_segment(
      data = lens,
      aes(x = 0, xend = prot_len, y = y, yend = y),
      linewidth = 0.9, color = "grey85", lineend = "round"
    ) +
    geom_rect(
      data = domain_df,
      aes(xmin = xmin, xmax = xmax, ymin = y - 0.28, ymax = y + 0.28, fill = domain),
      color = "grey20", linewidth = 0.12
    ) +
    scale_fill_manual(values = pal_domain, drop = FALSE) +
    scale_x_continuous(limits = c(0, xmax_len), expand = c(0, 0))
} else {
  p_domain <- p_domain + theme_void()  # 没数据就干净占位
}

# ---------------------------
# Panel3: motif track (protein aa)
# ---------------------------
p_motif <- ggplot() +
  ggtitle("Motif Analysis") +
  scale_y_continuous(limits = y_rng, expand = c(0, 0)) +
  theme_minimal(base_size = 10) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = margin(5, 5, 5, 5)
  )

p_motif <- p_motif +
  geom_segment(
    data = lens,
    aes(x = 0, xend = prot_len, y = y, yend = y),
    linewidth = 0.9, color = "grey85", lineend = "round"
  ) +
  scale_x_continuous(limits = c(0, xmax_len), expand = c(0, 0))

if (nrow(motif) > 0) {
  p_motif <- p_motif +
    geom_rect(
      data = motif,
      aes(xmin = xmin, xmax = xmax, ymin = y - 0.28, ymax = y + 0.28, fill = motif),
      color = "grey20", linewidth = 0.12
    ) +
    scale_fill_manual(values = pal_motif, drop = FALSE)
} else {
  # 没 motif 也保留基线
  p_motif <- p_motif
}

# ---------------------------
# Panel4: gene structure (bp, per-gene relative)
# ---------------------------
p_gene <- ggplot() +
  ggtitle("Gene Structure") +
  scale_y_continuous(limits = y_rng, expand = c(0, 0)) +
  theme_minimal(base_size = 10) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    axis.title = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin = margin(5, 5, 5, 5)
  )

if (nrow(gene_df) > 0) {
  xmax_gene <- max(gene_df$xmax, na.rm = TRUE)

  p_gene <- p_gene +
    # baseline
    geom_segment(
      data = tip_pos,
      aes(x = 0, xend = xmax_gene, y = y, yend = y),
      linewidth = 0.5, color = "grey92", lineend = "round"
    ) +
    # intron connectors
    geom_segment(
      data = intron_df,
      aes(x = x1, xend = x2, y = y, yend = y),
      linewidth = 0.35, color = "grey60"
    ) +
    # blocks
    geom_rect(
      data = gene_df,
      aes(xmin = xmin, xmax = xmax, ymin = y - 0.28, ymax = y + 0.28, fill = feature),
      color = NA
    ) +
    scale_fill_manual(values = pal_gene, drop = FALSE) +
    scale_x_continuous(limits = c(0, xmax_gene), expand = c(0, 0))
} else {
  p_gene <- p_gene + theme_void()
}

# ---------------------------
# combine + collect legends
# ---------------------------
# combine + collect legends (add spacer between tree and domain)
# ---------------------------
final <- (p_tree + plot_spacer() + p_domain + p_motif + p_gene) +
  plot_layout(ncol = 5, widths = c(4.2, 0.00, 1.10, 1.55, 1.55), guides = "collect") &
  theme(
    legend.position = "right",
    legend.title = element_text(size = 9, face = "bold"),
    legend.text  = element_text(size = 8),
    legend.key.height = unit(4.5, "mm"),
    legend.key.width  = unit(4.5, "mm")
  )

ggsave(
  out, final,
  width = 18,
  height = max(7, length(tr$tip.label) * 0.12),
  limitsize = FALSE
)


message("Done: ", out)

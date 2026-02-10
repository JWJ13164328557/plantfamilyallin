suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
  library(dplyr)
  library(readr)
  library(tibble)
})

args <- commandArgs(trailingOnly = TRUE)
getArg <- function(flag) {
  i <- match(flag, args)
  if (is.na(i) || i == length(args)) stop(paste("Missing", flag))
  args[i + 1]
}

chr_file <- getArg("--chrlen")
bed_file <- getArg("--bed")
family   <- getArg("--family")
out_pdf  <- getArg("--out")

# ---------------------------
# 读入
chr_df <- read_tsv(chr_file, col_names = FALSE, show_col_types = FALSE) |>
  setNames(c("chr", "len")) |>
  mutate(
    chr = as.character(.data$chr),
    len = as.numeric(.data$len)
  )

bed_df <- read_tsv(bed_file, col_names = FALSE, show_col_types = FALSE) |>
  setNames(c("chr", "start", "end", "id", "score", "strand")) |>
  mutate(
    chr   = as.character(.data$chr),
    start = as.numeric(.data$start),
    end   = as.numeric(.data$end),
    pos   = (.data$start + .data$end) / 2
  )

# 染色体顺序按 chr.length 文件顺序
chr_levels <- chr_df$chr
chr_df$chr <- factor(chr_df$chr, levels = chr_levels)
bed_df$chr <- factor(bed_df$chr, levels = chr_levels)

# ---------------------------
# 统一用 Mb（更像示例图）
chr_df <- chr_df |> mutate(lenMb = .data$len / 1e6)
bed_df <- bed_df |> mutate(posMb = .data$pos / 1e6) |> arrange(.data$chr, .data$posMb)

ymax <- max(chr_df$lenMb, na.rm = TRUE)

# ---------------------------
# 分面布局：两行
nchr <- nrow(chr_df)
nrow_facets <- 2
ncol_facets <- ceiling(nchr / nrow_facets)

# 每行最左侧 panel 的 chr（用于画标尺，只在这两个 panel 画）
left_chr_1 <- chr_levels[1]
left_chr_2 <- if (nchr >= (ncol_facets + 1)) chr_levels[ncol_facets + 1] else NA_character_
left_chrs  <- na.omit(c(left_chr_1, left_chr_2))

# 标尺刻度（尽量接近示例：约 10 段）
tick_step <- ymax / 10
# 把步长“规整”成常见整数（比如 13 Mb 这种也会自然出现）
tick_step <- max(1, round(tick_step))
ticks <- seq(0, ceiling(ymax / tick_step) * tick_step, by = tick_step)

scalebar <- tibble(
  chr = factor(rep(left_chrs, each = length(ticks)), levels = chr_levels),
  y = rep(ticks, times = length(left_chrs))
)

# ---------------------------
# 几个可调的审美参数（你可以按期刊风格微调）
chr_x      <- 1.00  # 染色体主轴 x
tick_x0    <- 0.82  # 基因短横线起点
tick_x1    <- 1.08  # 基因短横线终点
label_x    <- 1.16  # 基因标签放置的 x
scale_x    <- 0.55  # 标尺竖线 x（只在左侧 panel 画）
scale_w    <- 0.05  # 标尺刻度横线长度
chr_lwd    <- 6.0   # 染色体线宽（圆头）
gene_lwd   <- 0.35  # 基因短横线线宽
seg_lwd    <- 0.25  # 引线线宽
label_size <- 2.6   # 标签字号（ggplot size 是 mm-ish，相对值）
strip_size <- 12    # 染色体标题字号

# 输出尺寸
w <- max(10, min(18, 2.2 * ncol_facets))
h <- 7.5
grDevices::cairo_pdf(out_pdf, width = w, height = h, family = "sans")

p <- ggplot() +
  # 染色体“圆头竖线”
  geom_segment(
    data = chr_df,
    aes(x = chr_x, xend = chr_x, y = 0, yend = lenMb),
    linewidth = chr_lwd,
    lineend = "round",
    color = "#1B7A3A"  # 深绿（接近示例）
  ) +
  # 基因位置短横线
  geom_segment(
    data = bed_df,
    aes(x = tick_x0, xend = tick_x1, y = posMb, yend = posMb),
    linewidth = gene_lwd,
    color = "black"
  ) +
  # 基因标签（沿 y 方向避让）
  geom_text_repel(
    data = bed_df,
    aes(x = label_x, y = posMb, label = id),
    direction = "y",
    hjust = 0,
    size = label_size,
    min.segment.length = 0,
    segment.color = "black",
    segment.size = seg_lwd,
    box.padding = 0.18,
    point.padding = 0.05,
    force = 0.8,
    max.overlaps = Inf
  ) +
  # 只在每行最左侧 panel 画 Mb 标尺：竖线
  geom_segment(
    data = distinct(scalebar, chr),
    aes(x = scale_x, xend = scale_x, y = min(ticks), yend = max(ticks)),
    linewidth = 0.6,
    color = "black"
  ) +
  # 标尺刻度线
  geom_segment(
    data = scalebar,
    aes(x = scale_x, xend = scale_x - scale_w, y = y, yend = y),
    linewidth = 0.6,
    color = "black"
  ) +
  # 标尺文字（示例是 “13 Mb” 这种）
  geom_text(
    data = scalebar,
    aes(x = scale_x - scale_w - 0.02, y = y, label = paste0(y, " Mb")),
    hjust = 1,
    vjust = 0.5,
    size = 3.2,
    color = "black"
  ) +
  facet_wrap(~chr, nrow = nrow_facets) +
  coord_cartesian(ylim = c(0, ymax), clip = "off") +
  scale_x_continuous(limits = c(0.35, 1.70), breaks = NULL, expand = c(0, 0)) +
  scale_y_continuous(breaks = NULL, expand = c(0.02, 0)) +
  labs(title = NULL) +
  theme_void(base_size = 12) +
  theme(
    strip.text = element_text(face = "bold", size = strip_size, margin = margin(b = 6)),
    panel.spacing = unit(1.2, "lines"),
    plot.margin = margin(t = 8, r = 30, b = 8, l = 20)
  )

print(p)
dev.off()
message("Done: ", out_pdf)

#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
getArg <- function(flag, default = NULL) {
  i <- match(flag, args)
  if (is.na(i)) return(default)
  if (i == length(args)) stop(paste("Missing", flag))
  args[i + 1]
}

treefile   <- getArg("--tree")
out        <- getArg("--out")
family     <- getArg("--family", "Family")

# ---- parameters ----
k_clade     <- as.integer(getArg("--k", "6"))                 # 默认 6 组
open_angle  <- as.numeric(getArg("--open_angle", "0"))        # 0=完整圆；想留缺口可设 8~15
show_label  <- tolower(getArg("--show_label", "true")) %in% c("1","true","t","yes","y")

# 线条/字号（论文风格）
tree_lwd    <- as.numeric(getArg("--tree_lwd", "0.35"))
tip_size    <- as.numeric(getArg("--tip_size", "0.90"))
label_size  <- as.numeric(getArg("--label_size", "1.80"))
hilight_alpha <- as.numeric(getArg("--hilight_alpha", "0.18"))

suppressPackageStartupMessages({
  library(ape)
  library(ggplot2)
  library(ggtree)
  library(phangorn)
  library(dplyr)
  library(tibble)
  library(ggtreeExtra)
})

# ---------------------------
# paper-friendly palette (soft, distinct)
# ---------------------------
paper_group_palette <- function(k) {
  base <- c(
    "#66C2A5", # teal
    "#FC8D62", # orange
    "#8DA0CB", # blue-purple
    "#E78AC3", # pink
    "#A6D854", # green
    "#FFD92F", # yellow
    "#E5C494", # tan
    "#B3B3B3"  # grey
  )
  if (k <= length(base)) return(base[1:k])
  extra <- grDevices::hcl.colors(k - length(base), palette = "Dynamic")
  c(base, extra)
}

# ---------------------------
# tree + rooting
# ---------------------------
tr <- read.tree(treefile)

# 没有根就中点定根（fan/fan-like 画法更稳定）
if (!is.rooted(tr)) {
  tr <- phangorn::midpoint(tr)
}

ntip <- length(tr$tip.label)

# ---------------------------
# build child list + root
# ---------------------------
edge <- tr$edge
parents <- edge[, 1]
children <- edge[, 2]

root <- setdiff(unique(parents), unique(children))
if (length(root) != 1) {
  # 兜底：ape 通常根在 Ntip+1（但不绝对）
  root <- ntip + 1
}

child_list <- split(children, parents)

# ---------------------------
# memoized tipset(node): return tip indices under node
# ---------------------------
.tip_cache <- new.env(parent = emptyenv())

get_tip_idx <- function(node) {
  key <- as.character(node)
  if (exists(key, envir = .tip_cache, inherits = FALSE)) {
    return(get(key, envir = .tip_cache, inherits = FALSE))
  }
  res <- NULL
  if (node <= ntip) {
    res <- node
  } else {
    kids <- child_list[[key]]
    if (is.null(kids) || length(kids) == 0) {
      res <- integer(0)
    } else {
      res <- sort(unique(unlist(lapply(kids, get_tip_idx))))
    }
  }
  assign(key, res, envir = .tip_cache)
  res
}

tip_count <- function(node) length(get_tip_idx(node))

# ---------------------------
# Partition the tree into exactly K monophyletic clades
# Strategy: start from root as 1 group, iteratively split the largest internal clade
# (binary tree: each split increases group count by 1)
# ---------------------------
partition_into_k_clades <- function(k) {
  groups <- c(root)

  # 防止死循环
  guard <- 0
  while (length(groups) < k) {
    guard <- guard + 1
    if (guard > 100000) stop("Partition loop guard triggered (unexpected tree structure).")

    # 可分裂的候选：内部节点且有>=2个孩子
    candidates <- groups[groups > ntip]
    candidates <- candidates[sapply(candidates, function(nd) {
      kids <- child_list[[as.character(nd)]]
      !is.null(kids) && length(kids) >= 2
    })]

    if (length(candidates) == 0) break

    # 选“tips最多”的那个 clade 去 split
    sizes <- sapply(candidates, tip_count)
    nd_to_split <- candidates[which.max(sizes)]
    kids <- child_list[[as.character(nd_to_split)]]

    # 如果是多叉节点，split 会一下增加很多组，可能超过 k
    # 我们优先找二叉 split；如果多叉会超，就尝试换下一个候选
    if (length(kids) != 2 && (length(groups) - 1 + length(kids) > k)) {
      # 尝试找另一个不会超的
      ok <- FALSE
      ord <- order(sizes, decreasing = TRUE)
      for (j in ord) {
        cand <- candidates[j]
        kids2 <- child_list[[as.character(cand)]]
        if (length(kids2) == 2) {
          nd_to_split <- cand
          kids <- kids2
          ok <- TRUE
          break
        }
      }
      # 如果找不到二叉，只能接受多叉 split（组数可能略超过 k）
      if (!ok) {
        # 接受多叉，但随后截断到 k（最后几组可能非常小）
        # 注意：多叉树情况下“严格等于k且每组单系”不总能保证
      }
    }

    # replace nd_to_split with its children
    groups <- c(groups[groups != nd_to_split], kids)

    # 多叉导致超过 k：截断到 k（尽量保留大的 clade）
    if (length(groups) > k) {
      # 按 tips 数排序，保留前 k 个
      gs <- sapply(groups, tip_count)
      groups <- groups[order(gs, decreasing = TRUE)][1:k]
      break
    }
  }

  groups
}

group_nodes <- partition_into_k_clades(k_clade)

# 如果因为树结构导致不足 K，做提示但继续画
if (length(group_nodes) < k_clade) {
  warning(sprintf("Only %d groups could be formed (tree structure may be too shallow/multifurcating).", length(group_nodes)))
  k_clade <- length(group_nodes)
}

group_labels <- paste0("Group", seq_len(k_clade))
pal_groups <- paper_group_palette(k_clade)
names(pal_groups) <- group_labels

# tip -> group mapping (each tip belongs to exactly one group clade)
tip_group <- rep(NA_character_, ntip)
for (i in seq_along(group_nodes)) {
  nd <- group_nodes[i]
  tips_i <- get_tip_idx(nd)
  tip_group[tips_i] <- group_labels[i]
}
# 兜底：如果有 NA（理论上不该），全部归到最后一组
if (any(is.na(tip_group))) {
  tip_group[is.na(tip_group)] <- group_labels[length(group_labels)]
}

dat_group <- tibble(
  label = tr$tip.label,
  major_group = tip_group
)

# ---------------------------
# fan tree (circular branches)
# branch.length="none" => 更接近你示例那种“圆形均匀”的内部结构
# ---------------------------
p <- ggtree(
  tr,
  layout = "fan",
  open.angle = open_angle,
  branch.length = "none",
  size = tree_lwd
)

# tree radius (for offsets)
xmax_tree <- max(p$data$x, na.rm = TRUE)

# ---------------------------
# clade background: geom_hilight (SAFE; no central artifact)
# ---------------------------
for (i in seq_along(group_nodes)) {
  p <- p + ggtree::geom_hilight(
    node = group_nodes[i],
    fill = pal_groups[group_labels[i]],
    alpha = hilight_alpha
  )
}

# redraw tree on top (clean black)
p <- p + geom_tree(size = tree_lwd, color = "black")

# ---------------------------
# tip points: color by group
# (use explicit data join to avoid ggtree column conflicts)
# ---------------------------
tip_df <- p$data %>%
  filter(isTip) %>%
  select(label, x, y, angle) %>%
  left_join(dat_group, by = "label")

p <- p + geom_point(
  data = tip_df,
  aes(x = x, y = y, color = major_group),
  size = tip_size,
  alpha = 0.95
)

# ---------------------------
# outer ring: same colors as groups
# ---------------------------
# 放到树外一点
ring_offset <- xmax_tree * 0.004
ring_width  <- xmax_tree * 0.004

if ("pwidth" %in% names(formals(ggtreeExtra::geom_fruit))) {
  p <- p + geom_fruit(
    data = dat_group,
    geom = geom_tile,
    mapping = aes(y = label, x = 1, fill = major_group),
    offset = ring_offset,
    pwidth = ring_width
  )
} else {
  p <- p + geom_fruit(
    data = dat_group,
    geom = geom_tile,
    mapping = aes(y = label, x = 1, fill = major_group),
    offset = ring_offset,
    width = ring_width
  )
}

# ---------------------------
# tip labels (gene names)
# ---------------------------
if (show_label) {
  # 放在 ring 外侧
  lab_off <- ring_offset + ring_width + xmax_tree * 0.0025
  p <- p + geom_tiplab2(
    size = label_size,
    offset = lab_off,
    align = FALSE
  )
}

# ---------------------------
# scales + theme
# ---------------------------
p <- p +
  scale_fill_manual(values = pal_groups, breaks = group_labels) +
  scale_color_manual(values = pal_groups, breaks = group_labels) +
  ggtitle(paste0(family, " Phylogenetic Tree")) +
  theme_void() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    legend.position = "right",
    legend.title = element_blank(),
    legend.text  = element_text(size = 10)
  )
fig_w <- if (show_label) 24 else 12
fig_h <- if (show_label) 24 else 12
ggsave(out, p, width = fig_w, height = fig_h, limitsize = FALSE)
message("Done: ", out)

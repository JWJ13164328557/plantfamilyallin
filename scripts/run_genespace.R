#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(GENESPACE)
  library(data.table)
})

# ---- Force a safe read_bed() that returns a data.table with an 'id' column ----
read_bed <- function(fpath) {
  bd <- data.table::fread(fpath, header = FALSE, sep = "\t", data.table = TRUE)
  if (ncol(bd) < 4) stop("GENESPACE bed must have >=4 columns (chr,start,end,id): ", fpath)
  # normalize names and create id column
  data.table::setnames(bd, paste0("V", seq_len(ncol(bd))))
  bd[, id := V4]
  bd
}


opt_list <- list(
  make_option(c("--wd"), type="character", help="GENESPACE working directory (must contain bed/ and peptide/)"),
  make_option(c("--path2mcscanx"), type="character", help="Directory containing MCScanX executable"),
  make_option(c("--genomes"), type="character", help="Comma-separated genome IDs, e.g. SL,AT,PEP"),
  make_option(c("--ploidy"), type="integer", default=1, help="Ploidy"),
  make_option(c("--overwrite"), action="store_true", default=FALSE, help="Overwrite existing results"),
  make_option(c("--ref"), type="character", default=NULL, help="Reference genome ID for riparian plot"),
  make_option(c("--out_rds"), type="character", help="Output RDS for gsParam"),
  make_option(c("--out_pdf"), type="character", help="Output PDF for riparian plot"),
  make_option(c("--skip_parse"), action="store_true", default=TRUE,
              help="Skip parse_annotations() and assume wd/bed + wd/peptide already prepared (recommended)")
)

opt <- parse_args(OptionParser(option_list=opt_list))

req <- c("wd","path2mcscanx","genomes","out_rds","out_pdf")
missing_args <- req[ sapply(req, function(x) is.null(opt[[x]]) || opt[[x]]=="") ]
if (length(missing_args) > 0) stop(paste("Missing required args:", paste(missing_args, collapse=", ")))

wd <- opt$wd
path2mcscanx <- opt$path2mcscanx
genomes2run <- trimws(strsplit(opt$genomes, ",")[[1]])

cat("== GENESPACE inputs ==\n")
cat("wd:", wd, "\n")
cat("path2mcscanx:", path2mcscanx, "\n")
cat("genomes:", paste(genomes2run, collapse=", "), "\n")
cat("skip_parse:", opt$skip_parse, "\n")
cat("overwrite:", opt$overwrite, "\n")
cat("ploidy:", opt$ploidy, "\n\n")

if (!dir.exists(path2mcscanx)) stop(paste("path2mcscanx dir not found:", path2mcscanx))

# Check MCScanX executable
mcscanx_bin <- file.path(path2mcscanx, "MCScanX")
if (!file.exists(mcscanx_bin)) {
  stop(paste("MCScanX executable not found at:", mcscanx_bin,
             "\nNOTE: --path2mcscanx must be a DIRECTORY that contains MCScanX"))
}

dir.create(wd, recursive=TRUE, showWarnings=FALSE)
dir.create(file.path(wd, "bed"), recursive=TRUE, showWarnings=FALSE)
dir.create(file.path(wd, "peptide"), recursive=TRUE, showWarnings=FALSE)

# If skipping parse, verify prepared inputs exist
if (opt$skip_parse) {
  missing_files <- c()
  for (g in genomes2run) {
    bed <- file.path(wd, "bed", paste0(g, ".bed"))
    pep <- file.path(wd, "peptide", paste0(g, ".fa"))
    if (!file.exists(bed)) missing_files <- c(missing_files, bed)
    if (!file.exists(pep)) missing_files <- c(missing_files, pep)
  }
  if (length(missing_files) > 0) {
    stop(paste("Missing prepared GENESPACE inputs (expected in wd/bed and wd/peptide):\n",
               paste(missing_files, collapse="\n")))
  }
} else {
  stop("skip_parse=FALSE is not implemented in this streamlined script. Use parse_annotations interactively if needed.")
}

cat("== Step1: init_genespace ==\n")
gpar <- init_genespace(
  wd = wd,
  ploidy = opt$ploidy,
  path2mcscanx = path2mcscanx
)

cat("== Step2: run_genespace ==\n")
out <- run_genespace(gpar, overwrite = opt$overwrite)
cat("run_genespace done.\n\n")

saveRDS(out, opt$out_rds)
cat("Saved gsParam RDS:", opt$out_rds, "\n")

refGenome <- opt$ref
if (is.null(refGenome) || !(refGenome %in% genomes2run)) {
  refGenome <- genomes2run[[1]]
  cat("ref not provided or invalid; using:", refGenome, "\n")
}

cat("== Step3: plot_riparian ==\n")

# ---- Fix: force refChr/chr columns to character (avoid data.table join type mismatch) ----
force_chr_cols_character <- function(x) {
  chr_cols <- c("refChr", "chr")
  if (data.table::is.data.table(x)) {
    for (cc in chr_cols) {
      if (cc %in% names(x)) data.table::set(x, j = cc, value = as.character(x[[cc]]))
    }
    return(x)
  }
  if (is.data.frame(x)) {
    for (cc in chr_cols) {
      if (cc %in% names(x)) x[[cc]] <- as.character(x[[cc]])
    }
    return(x)
  }
  if (is.list(x)) {
    for (nm in names(x)) x[[nm]] <- force_chr_cols_character(x[[nm]])
    return(x)
  }
  x
}

# 1) fix in-memory gsParam
out <- force_chr_cols_character(out)

# 2) fix on-disk tables (plot_riparian sometimes reloads these)
rp_files <- c(
  file.path(wd, "riparian", "refPhasedBlkCoords.txt"),
  file.path(wd, "results",  "blkCoords.txt"),
  file.path(wd, "results",  "combBed.txt")
)

for (fp in rp_files) {
  if (file.exists(fp)) {
    dt <- data.table::fread(fp, sep = "\t", data.table = TRUE)
    dt <- force_chr_cols_character(dt)
    data.table::fwrite(dt, fp, sep = "\t")
  }
}

# 3) patch fread inside plot_riparian (last line of defense)
.dt_ns <- asNamespace("data.table")
.fread_orig <- get("fread", envir = .dt_ns)
unlockBinding("fread", .dt_ns)
assign("fread", function(...) {
  dt <- .fread_orig(...)
  if (is.data.frame(dt)) {
    if ("refChr" %in% names(dt)) dt$refChr <- as.character(dt$refChr)
    if ("chr"    %in% names(dt)) dt$chr    <- as.character(dt$chr)
  }
  dt
}, envir = .dt_ns)
lockBinding("fread", .dt_ns)

pdf(opt$out_pdf, width=10, height=7)
plot_riparian(gsParam = out, refGenome = refGenome, forceRecalcBlocks = FALSE)
dev.off()

cat("Saved riparian PDF:", opt$out_pdf, "\n")
cat("\nALL DONE.\n")

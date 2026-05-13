# =============================================================================
# GREAT Enrichment Analysis: Zoonomia cCRE Conservation Groups
# =============================================================================

library(rGREAT)
library(GenomicRanges)
library(rtracklayer)

# -----------------------------------------------------------------------------
# Configuration
# -----------------------------------------------------------------------------

BASE     <- "/project/home/p201120/ryan/cCRE_pipeline/zoonomia_out"
GENOME   <- "hg38"
OUT_DIR  <- file.path(BASE, "great_results")

GROUP_FILES <- c(
  G1 = "G1_final.sorted.bed",
  G2 = "G2_final.sorted.bed",
  G3 = "G3_final.sorted.bed"
)
BG_FILE <- "colon_epithelial_ccres_bg.sorted.bed"

setwd(BASE)
dir.create(OUT_DIR, showWarnings = FALSE)

# -----------------------------------------------------------------------------
# Input validation
# -----------------------------------------------------------------------------

all_files <- c(GROUP_FILES, BG_FILE)
missing   <- all_files[!file.exists(all_files)]
if (length(missing) > 0) {
  stop("Missing input files:\n  ", paste(missing, collapse = "\n  "))
}

# Quick check that chromosomes are chr-prefixed (required for hg38)
check_chr_prefix <- function(path) {
  gr <- import(path)
  if (!any(grepl("^chr", as.character(seqnames(gr))))) {
    warning(path, " does not appear to have chr-prefixed chromosomes — check compatibility with hg38")
  }
  gr
}

# -----------------------------------------------------------------------------
# Load data
# -----------------------------------------------------------------------------

message("Loading background regions...")
bg_gr <- check_chr_prefix(BG_FILE)
message(sprintf("  Background: %d regions", length(bg_gr)))

group_grs <- lapply(GROUP_FILES, function(f) {
  gr <- check_chr_prefix(f)
  message(sprintf("  %s: %d regions", f, length(gr)))
  gr
})

# -----------------------------------------------------------------------------
# Run GREAT
# -----------------------------------------------------------------------------

run_great <- function(name, gr, bg, genome, out_dir) {
  message(sprintf("\n[%s] Submitting GREAT job (%d regions)...", name, length(gr)))

  job <- tryCatch(
    submitGreatJob(gr, genome = genome, bg = bg),
    error = function(e) {
      message(sprintf("  ERROR submitting %s: %s", name, conditionMessage(e)))
      return(NULL)
    }
  )

  if (is.null(job)) return(NULL)

  message(sprintf("[%s] Fetching enrichment tables...", name))
  tables <- getEnrichmentTables(job, availableOntologies(job))

  # Save each ontology as a TSV
  for (ont in names(tables)) {
    clean_name <- gsub("[^A-Za-z0-9_]", "_", ont)
    out_path   <- file.path(out_dir, sprintf("%s_%s.tsv", name, clean_name))
    write.table(tables[[ont]], out_path, sep = "\t", quote = FALSE, row.names = FALSE)
  }

  message(sprintf("[%s] Done. Results written to %s/", name, out_dir))
  tables
}

results <- mapply(
  run_great,
  name    = names(group_grs),
  gr      = group_grs,
  MoreArgs = list(bg = bg_gr, genome = GENOME, out_dir = OUT_DIR),
  SIMPLIFY = FALSE
)

# Save full results object
saveRDS(results, file.path(OUT_DIR, "great_results.rds"))
message("\nFull results saved to great_results.rds")

# -----------------------------------------------------------------------------
# Quick summary: top 10 GO Biological Process hits per group
# -----------------------------------------------------------------------------

message("\n========== Top 10 GO Biological Process hits ==========")
for (grp in names(results)) {
  if (is.null(results[[grp]])) next
  go_bp <- results[[grp]][["GO Biological Process"]]
  if (is.null(go_bp) || nrow(go_bp) == 0) {
    message(sprintf("\n[%s] No GO BP results.", grp))
    next
  }
  message(sprintf("\n--- %s ---", grp))
  print(head(go_bp[order(go_bp$Hyper_Adjp_BH), ], 10))
}

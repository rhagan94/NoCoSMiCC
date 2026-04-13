#### Script 1 - ArchR_Analysis.R ####

# Script for ArchR project creation using ENCODE and HuBMAP scATAC datasets, identification of epithelial cells and peak calling
# Adapted from Hickey et al. 2023 (doi:10.1038/s41586-023-05915-x)

#######################################
# Setup
#######################################
# Load packages
library(ArchR)
library(Seurat)
library(BSgenome.Hsapiens.UCSC.hg38)
library(parallel)
'%notin%' <- Negate('%in%')

# Load genome annotations
addArchRGenome("hg38")

# Set threads
addArchRThreads(threads = 64)

# Set working directory
setwd("/project/home/p201120/ryan/cCRE_pipeline/")

# Define arrow/frag file directories
arrowdir <- "/project/home/p201120/ryan/cCRE_pipeline/outputs/ArrowFiles/"
fragDir <- "/project/home/p201120/ryan/cCRE_pipeline/files/scATAC/"

#######################################
# Define paths to fragment files
#######################################

inputFiles1K <- c(
  "B006-A-001"    = paste0(fragDir, "B006-A-001_atac_fragments.tsv.gz"),
  "B006-A-101"    = paste0(fragDir, "B006-A-101_atac_fragments.tsv.gz"),
  "B006-A-201-R2" = paste0(fragDir, "B006-A-201-R2_atac_fragments.tsv.gz"),
  "ENCSR916RYB"   = paste0(fragDir, "ENCSR916RYB.fragments.tsv.gz"),
  "ENCSR904WIW"   = paste0(fragDir, "ENCSR904WIW.fragments.tsv.gz"),
  "ENCSR830FPR"   = paste0(fragDir, "ENCSR830FPR.fragments.tsv.gz")
)

inputFiles1p5K <- c(
  "B001-A-302" = paste0(fragDir, "B001-A-302_atac_fragments.tsv.gz"),
  "B001-A-401" = paste0(fragDir, "B001-A-401_atac_fragments.tsv.gz"),
  "B001-A-406" = paste0(fragDir, "B001-A-406_atac_fragments.tsv.gz"),
  "B001-A-501" = paste0(fragDir, "B001-A-501_atac_fragments.tsv.gz")
)

inputFiles2K <- c(
  "B004-A-004" = paste0(fragDir, "B004-A-004_atac_fragments.tsv.gz"),
  "B008-A-001" = paste0(fragDir, "B008-A-001_atac_fragments.tsv.gz"),
  "ENCSR997YNO" = paste0(fragDir, "ENCSR997YNO.fragments.tsv.gz"),
  "ENCSR007QIO" = paste0(fragDir, "ENCSR007QIO.fragments.tsv.gz"),
  "ENCSR349XKD" = paste0(fragDir, "ENCSR349XKD.fragments.tsv.gz"),
  "ENCSR434SXE" = paste0(fragDir, "ENCSR434SXE.fragments.tsv.gz"),
  "ENCSR388NCA" = paste0(fragDir, "ENCSR388NCA.fragments.tsv.gz"),
  "ENCSR506YMX" = paste0(fragDir, "ENCSR506YMX.fragments.tsv.gz"),
  "ENCSR367GKP" = paste0(fragDir, "ENCSR367GKP.fragments.tsv.gz")
)

inputFiles3K <- c(
  "B001-A-301" = paste0(fragDir, "B001-A-301_atac_fragments.tsv.gz"),
  "B005-A-001" = paste0(fragDir, "B005-A-001_atac_fragments.tsv.gz"),
  "B005-A-002" = paste0(fragDir, "B005-A-002_atac_fragments.tsv.gz"),
  "B005-A-101" = paste0(fragDir, "B005-A-101_atac_fragments.tsv.gz"),
  "B005-A-201" = paste0(fragDir, "B005-A-201_atac_fragments.tsv.gz"),
  "B006-A-002" = paste0(fragDir, "B006-A-002_atac_fragments.tsv.gz"),
  "B006-A-201" = paste0(fragDir, "B006-A-201_atac_fragments.tsv.gz"),
  "B008-A-002" = paste0(fragDir, "B008-A-002_atac_fragments.tsv.gz"),
  "B008-A-101" = paste0(fragDir, "B008-A-101_atac_fragments.tsv.gz"),
  "B008-A-201" = paste0(fragDir, "B008-A-201_atac_fragments.tsv.gz"),
  "B010-A-001" = paste0(fragDir, "B010-A-001_atac_fragments.tsv.gz"),
  "B010-A-002" = paste0(fragDir, "B010-A-002_atac_fragments.tsv.gz"),
  "B010-A-101" = paste0(fragDir, "B010-A-101_atac_fragments.tsv.gz"),
  "B011-A-001" = paste0(fragDir, "B011-A-001_atac_fragments.tsv.gz"),
  "B011-A-002" = paste0(fragDir, "B011-A-002_atac_fragments.tsv.gz"),
  "B011-A-201" = paste0(fragDir, "B011-A-201_atac_fragments.tsv.gz"),
  "B012-A-001" = paste0(fragDir, "B012-A-001_atac_fragments.tsv.gz"),
  "B012-A-002" = paste0(fragDir, "B012-A-002_atac_fragments.tsv.gz"),
  "B012-A-101" = paste0(fragDir, "B012-A-101_atac_fragments.tsv.gz")
)

inputFiles4K <- c(
  "B004-A-004-R2" = paste0(fragDir, "B004-A-004-R2_atac_fragments.tsv.gz"),
  "B004-A-204"    = paste0(fragDir, "B004-A-204_atac_fragments.tsv.gz"),
  "B009-A-101"    = paste0(fragDir, "B009-A-101_atac_fragments.tsv.gz"),
  "B010-A-201"    = paste0(fragDir, "B010-A-201_atac_fragments.tsv.gz"),
  "B011-A-101"    = paste0(fragDir, "B011-A-101_atac_fragments.tsv.gz")
)

inputFiles5K <- c(
  "B004-A-008" = paste0(fragDir, "B004-A-008_atac_fragments.tsv.gz"),
  "B009-A-001" = paste0(fragDir, "B009-A-001_atac_fragments.tsv.gz")
)

inputFiles6K <- c(
  "B012-A-201" = paste0(fragDir, "B012-A-201_atac_fragments.tsv.gz")
)

#######################################
# Create arrow files
#######################################

ArrowFiles1K <- createArrowFiles(
  inputFiles   = inputFiles1K,
  sampleNames  = names(inputFiles1K),
  minFrags     = 1000,
  minTSS       = 5
)

ArrowFiles1p5K <- createArrowFiles(
  inputFiles   = inputFiles1p5K,
  sampleNames  = names(inputFiles1p5K),
  minFrags     = 1500,
  minTSS       = 5
)

ArrowFiles2K <- createArrowFiles(
  inputFiles   = inputFiles2K,
  sampleNames  = names(inputFiles2K),
  minFrags     = 2000,
  minTSS       = 5
)

ArrowFiles3K <- createArrowFiles(
  inputFiles   = inputFiles3K,
  sampleNames  = names(inputFiles3K),
  minFrags     = 3000,
  minTSS       = 5
)

ArrowFiles4K <- createArrowFiles(
  inputFiles   = inputFiles4K,
  sampleNames  = names(inputFiles4K),
  minFrags     = 4000,
  minTSS       = 5
)

ArrowFiles5K <- createArrowFiles(
  inputFiles   = inputFiles5K,
  sampleNames  = names(inputFiles5K),
  minFrags     = 5000,
  minTSS       = 5
)

ArrowFiles6K <- createArrowFiles(
  inputFiles   = inputFiles6K,
  sampleNames  = names(inputFiles6K),
  minFrags     = 6000,
  minTSS       = 5
)

# Define ENCODE arrows
ENCODE_Arrows <- c(
                  paste0(arrowdir, "ENCSR916RYB.arrow"), 
                  paste0(arrowdir,"ENCSR904WIW.arrow"), 
                  paste0(arrowdir, "ENCSR830FPR.arrow"),
                  paste0(arrowdir, "ENCSR997YNO.arrow"), 
                  paste0(arrowdir, "ENCSR007QIO.arrow"), 
                  paste0(arrowdir, "ENCSR349XKD.arrow"), 
                  paste0(arrowdir, "ENCSR434SXE.arrow"),
                  paste0(arrowdir, "ENCSR388NCA.arrow"), 
                  paste0(arrowdir, "ENCSR506YMX.arrow"), 
                  paste0(arrowdir, "ENCSR367GKP.arrow"))

# Compute doublet scores for ENCODE samples
doubletScores <- addDoubletScores(ENCODE_Arrows, k = 10, knnMethod = "UMAP", LSIMethod = 1, threads = 1)

#######################################
# Create ArchR project for ENCODE samples
#######################################

ENCODEproj <- ArchRProject(
  ArrowFiles        = ENCODE_Arrows,
  outputDirectory   = "/project/home/p201120/ryan/cCRE_pipeline/files/scATAC/ArchR_output"
)

saveArchRProject(
  ArchRProj       = ENCODEproj,
  outputDirectory = "/project/home/p201120/ryan/cCRE_pipeline/files/scATAC/ArchR_output",
  load            = FALSE,
  overwrite       = FALSE
)

#######################################
# Filter doublets
#######################################

ENCODEproj <- loadArchRProject("/project/home/p201120/ryan/cCRE_pipeline/files/scATAC/ArchR_output")
ENCODEproj <- filterDoublets(ENCODEproj)
#write.table(rownames(getCellColData(ENCODEproj)), "initial_post_filter_ENCODE_cells.txt")

#######################################
# Plot quality control metrics
#######################################

p1 <- plotGroups(
  ArchRProj = ENCODEproj, 
  groupBy = "Sample", 
  colorBy = "cellColData", 
  name = "TSSEnrichment",
  plotAs = "violin",
  alpha = 0.8,
  addBoxPlot = TRUE
)

p2 <- plotGroups(
  ArchRProj = ENCODEproj, 
  groupBy = "Sample", 
  colorBy = "cellColData", 
  name = "log10(nFrags)",
  plotAs = "ridges"
)

plotPDF(p1, p2,
  name = "ENCODEproj_QC",
  ArchRProj = ENCODEproj,
  addDOC = FALSE,
  width = 5, height = 5
)

#######################################
# Dim reduction, clustering and epithelial cluster identification
#######################################

# Define marker genes
markerSets <- list(
  General_Epithelial = c("EPCAM", "KRT8", "KRT18", "KRT19", "CDH1"),
  T_cells            = c("CD3D", "CD3E", "CD3G", "CD8A", "CD8B",
                         "TBX21", "IL7R", "CD4", "CD2"),
  B_cells            = c("PAX5", "MS4A1", "CD19", "IGLL5", "VPREB3"),
  Myeloid            = c("CD14", "CD68", "CSF1R", "LYZ", "ITGAM"),
  NK_cells           = c("NCAM1", "KLRD1", "NKG7", "GNLY"),
  Fibroblasts        = c("COL1A1", "COL1A2", "COL6A1", "COL6A2",
                         "FAP", "CBLN2", "SPOCK1", "ACSS3"),
  Endothelial        = c("PECAM1", "VWF", "CDH5", "CLDN5"),
  SmoothMuscle       = c("ACTA2", "MYH11", "TAGLN")
)

allSamples <- unique(ENCODEproj$Sample)

for(samp in allSamples) {

  cat("\n=============================\n")
  cat("Processing sample:", samp, "\n")
  cat("=============================\n")

  outDir <- paste0("ArchR_", samp)
  dir.create(outDir, showWarnings = FALSE)

  sampCells <- getCellNames(ENCODEproj)[ENCODEproj$Sample == samp]
  cat("Total cells:", length(sampCells), "\n")

  sampProj <- subsetArchRProject(
    ArchRProj = ENCODEproj,
    cells = sampCells,
    outputDirectory = outDir,
    force = TRUE
  )

  sampProj <- addIterativeLSI(
    ArchRProj = sampProj,
    useMatrix = "TileMatrix",
    name = "IterativeLSI",
    iterations = 2,
    clusterParams = list(resolution = c(0.2), sampleCells = 10000, n.start = 10),
    varFeatures = 15000,
    dimsToUse = 1:15,
    force = TRUE
  )

  sampProj <- addClusters(
    input = sampProj,
    reducedDims = "IterativeLSI",
    method = "Seurat",
    name = "Clusters",
    resolution = 0.4,
    force = TRUE
  )

  print(table(sampProj$Clusters))

  sampProj <- addUMAP(
    ArchRProj = sampProj,
    reducedDims = "IterativeLSI",
    name = "UMAP",
    nNeighbors = 30,
    minDist = 0.5,
    metric = "cosine",
    force = TRUE
  )

  # Overview UMAPs
  sample_umap <- plotEmbedding(
    ArchRProj = sampProj,
    colorBy = "cellColData",
    name = "Sample",
    embedding = "UMAP"
  )

  cluster_umap <- plotEmbedding(
    ArchRProj = sampProj,
    colorBy = "cellColData",
    name = "Clusters",
    embedding = "UMAP"
  )

  # Read GeneScoreMatrix once
  gsMatrix <- getMatrixFromProject(sampProj, useMatrix = "GeneScoreMatrix")
  allMarkers <- unlist(markerSets)
  availMarkers <- allMarkers[allMarkers %in% rowData(gsMatrix)$name]

  gsRows <- which(rowData(gsMatrix)$name %in% availMarkers)
  gsMat  <- assay(gsMatrix)[gsRows, ]
  rownames(gsMat) <- rowData(gsMatrix)$name[gsRows]
  rm(gsMatrix); gc()

  clusters <- sampProj$Clusters

  meanScores <- sapply(sort(unique(clusters)), function(cl) {
    cells <- which(clusters == cl)
    if(length(cells) == 1) gsMat[, cells]
    else rowMeans(gsMat[, cells])
  })
  colnames(meanScores) <- sort(unique(clusters))

  # Heatmap ordered by marker set
  orderedGenes <- unlist(markerSets)
  orderedGenes <- orderedGenes[orderedGenes %in% rownames(meanScores)]

  heatMatScaled <- t(scale(t(meanScores[orderedGenes, ])))
  heatMatScaled[is.nan(heatMatScaled)] <- 0

  geneSetAnnot <- data.frame(
    CellType = rep(names(markerSets), lengths(markerSets)),
    row.names = unlist(markerSets)
  )[orderedGenes, , drop = FALSE]

  p_heatmap <- pheatmap::pheatmap(
    heatMatScaled,
    annotation_row = geneSetAnnot,
    cluster_rows = FALSE,
    cluster_cols = TRUE,
    show_rownames = TRUE,
    fontsize_row = 7,
    fontsize_col = 9,
    main = paste0(samp, " | Gene Activity Scores per Cluster"),
    silent = TRUE
  )

  # Cluster scores table
  setScores <- sapply(names(markerSets), function(setName) {
    genes <- markerSets[[setName]][markerSets[[setName]] %in% rownames(gsMat)]
    if(length(genes) == 0) return(rep(NA, ncol(meanScores)))
    if(length(genes) == 1) return(meanScores[genes, ])
    colMeans(meanScores[genes, ])
  })

  setScoresDF <- round(as.data.frame(setScores), 3)
  setScoresDF$Cluster <- rownames(setScoresDF)

  cat("\nCluster scores for", samp, ":\n")
  print(setScoresDF)

  write.csv(setScoresDF,
            file = file.path(outDir, paste0("ClusterScores_", samp, ".csv")),
            row.names = TRUE)

  # Save PDF: page 1 = UMAPs, page 2 = heatmap
  pdf(file.path(outDir, paste0("MarkerAnalysis_", samp, ".pdf")),
      width = 15, height = 10)

  # Page 1: UMAPs
  print(cowplot::plot_grid(
    sample_umap, cluster_umap,
    ncol = 2,
    labels = c("Sample", "Clusters"),
    label_size = 12
  ))

  # Page 2: Heatmap on its own page
  grid::grid.newpage()
  grid::grid.draw(p_heatmap$gtable)

  dev.off()
  cat("Saved PDF:", file.path(outDir, paste0("MarkerAnalysis_", samp, ".pdf")), "\n")

  saveArchRProject(sampProj)

  rm(sampProj, gsMat, meanScores, setScores, setScoresDF,
     p_heatmap, heatMatScaled, geneSetAnnot, orderedGenes)
  gc()
}

#######################################
# Define epithelial clusters
#######################################

ENCODE_epi_Clusters <- list(
  "ENCSR997YNO" = c("C1", "C2"),
  "ENCSR830FPR" = c("C1", "C2", "C3", "C4"),
  "ENCSR349XKD" = c("C1", "C2"),
  "ENCSR434SXE" = c("C1"),
  "ENCSR506YMX" = c("C2"),
  "ENCSR904WIW" = c("C1", "C2", "C3", "C4", "C5","C6")
  # Samples with no epithelial clusters are omitted
)

allEpithelialCells <- c()

for(samp in names(ENCODE_epi_Clusters)) {
  
  cat("\nSubsetting epithelial cells from:", samp, "\n")
  
  # Load the already-processed sample project
  sampProj <- loadArchRProject(paste0("ArchR_", samp))
  
  # Get cells from epithelial clusters
  epiClusters <- ENCODE_epi_Clusters[[samp]]
  epiCells <- getCellNames(sampProj)[sampProj$Clusters %in% epiClusters]
  
  cat("Epithelial cells:", length(epiCells), "\n")
  allEpithelialCells <- c(allEpithelialCells, epiCells)
}

cat("\nTotal epithelial cells across ENCODE samples:", length(allEpithelialCells), "\n")

# Subset ENCODE project to epithelial cells
ENCODEProjEpi <- subsetArchRProject(
  ArchRProj = ENCODEproj,
  cells = allEpithelialCells,
  outputDirectory = "ArchR_Epithelial_ENCODE",
  force = TRUE
)

# Verify sample composition
table(ENCODEProjEpi$Sample)

#######################################
# Create combined epithelial project
#######################################

# Define HuBMAP arrows
HuBMAP_Arrows <- c(
                  paste0(arrowdir, "/B006-A-001.arrow"), paste0(arrowdir,"/B006-A-101.arrow"), paste0(arrowdir,"/B006-A-201-R2.arrow"),
                  paste0(arrowdir,"/B001-A-302.arrow"), paste0(arrowdir,"/B001-A-401.arrow"), paste0(arrowdir,"/B001-A-406.arrow"), 
                  paste0(arrowdir, "/B001-A-501.arrow"), paste0(arrowdir, "/B004-A-004.arrow"), paste0(arrowdir, "/B008-A-001.arrow"), 
                  paste0(arrowdir, "/B001-A-301.arrow"), paste0(arrowdir, "/B005-A-001.arrow"), paste0(arrowdir, "/B005-A-002.arrow"), 
                  paste0(arrowdir, "/B005-A-101.arrow"), paste0(arrowdir, "/B005-A-201.arrow"), paste0(arrowdir, "/B006-A-002.arrow"), 
                  paste0(arrowdir, "/B006-A-201.arrow"), paste0(arrowdir, "/B008-A-002.arrow"), paste0(arrowdir, "/B008-A-101.arrow"), 
                  paste0(arrowdir, "/B008-A-201.arrow"), paste0(arrowdir, "/B010-A-001.arrow"), paste0(arrowdir, "/B010-A-002.arrow"), 
                  paste0(arrowdir, "/B010-A-101.arrow"), paste0(arrowdir, "/B011-A-001.arrow"), paste0(arrowdir, "/B011-A-002.arrow"), 
                  paste0(arrowdir, "/B011-A-201.arrow"), paste0(arrowdir, "/B012-A-001.arrow"), paste0(arrowdir, "/B012-A-002.arrow"), 
                  paste0(arrowdir, "/B012-A-101.arrow"), paste0(arrowdir, "/B004-A-004-R2.arrow"), paste0(arrowdir, "/B004-A-204.arrow"), 
                  paste0(arrowdir, "/B009-A-101.arrow"), paste0(arrowdir, "/B010-A-201.arrow"), paste0(arrowdir, "/B011-A-101.arrow"), 
                  paste0(arrowdir, "/B004-A-008.arrow"), paste0(arrowdir, "/B009-A-001.arrow"), paste0(arrowdir, "/B012-A-201.arrow"))

# Create ArchR project for HuBMAP samples                                                               
HuBMAPproj <- ArchRProject(
  ArrowFiles        = HuBMAP_Arrows,
  outputDirectory   = "/project/home/p201120/ryan/cCRE_pipeline/outputs/HuBMAPProj"
)

# Read in the cell barcodes
HuBMAP_multiome_colon_epi_barcodes <- read.table(
  "/mnt/tier2/project/p201120/ryan/cCRE_pipeline/files/colon_epithelial_peak_matrix_cells.tsv",
  header = FALSE,
  stringsAsFactors = FALSE,
  comment.char = ""
)$V1

HuBMAP_Nonmultiome_colon_epi_barcodes <- read.table(
  "/mnt/tier2/project/p201120/ryan/cCRE_pipeline/files/non_multiome_colon_epithelial_peak_matrix_cells.tsv",
  header = FALSE,
  stringsAsFactors = FALSE,
  comment.char = ""
)$V1

# Combine HuBMAP multiome and non-multiome epithelial barcodes
allHuBMAPEpithelialCells <- union(HuBMAP_multiome_colon_epi_barcodes, HuBMAP_Nonmultiome_colon_epi_barcodes)

# Retain only barcodes present in the HuBMAP ArchR project
allHuBMAPEpithelialCells <- allHuBMAPEpithelialCells[allHuBMAPEpithelialCells %in% getCellNames(HuBMAPproj)]

# Subset HuBMAP project to epithelial cells
HuBMAPProjEpi <- subsetArchRProject(
  ArchRProj = HuBMAPproj,
  cells = allHuBMAPEpithelialCells,
  outputDirectory = "ArchR_Epithelial_HuBMAP",
  force = TRUE
)

# Get all arrow files from both projects
allArrows <- c(
  getArrowFiles(ENCODEProjEpi),
  getArrowFiles(HuBMAPProjEpi)
)

# Get all epithelial cell names from both projects
allEpiCells <- c(
  getCellNames(ENCODEProjEpi),
  getCellNames(HuBMAPProjEpi)
)

# Create combined epithelial project from all arrows
CombinedEpiProj <- ArchRProject(
  ArrowFiles      = allArrows,
  outputDirectory = "ArchR_Epithelial_Combined"
)

# Subset to only the epithelial cells from both projects
CombinedEpiProj <- subsetArchRProject(
  ArchRProj       = CombinedEpiProj,
  cells           = allEpiCells[allEpiCells %in% getCellNames(CombinedEpiProj)],
  outputDirectory = "ArchR_Epithelial_Combined",
  force           = TRUE
)

table(CombinedEpiProj$Sample)

#######################################
# Run LSI and Harmony on epithelial cells
#######################################

CombinedEpiProj <- addIterativeLSI(
  ArchRProj = CombinedEpiProj,
  useMatrix = "TileMatrix",
  name = "IterativeLSI",
  iterations = 3,
  sampleCellsPre = 25000,
  dimsToUse = 1:25,
  varFeatures = 15000,
  force = TRUE
)

CombinedEpiProj <- addHarmony(
  ArchRProj = CombinedEpiProj,
  reducedDims = "IterativeLSI",
  name = "Harmony",
  groupBy = "Sample",
  force = TRUE
)

CombinedEpiProj <- addClusters(
  input = CombinedEpiProj,
  reducedDims = "Harmony",
  resolution = 1.0,
  force = TRUE
)

CombinedEpiProj <- addUMAP(
  ArchRProj = CombinedEpiProj,
  reducedDims = "Harmony",
  name = "UMAP_Harmony",
  nNeighbors = 30,
  minDist = 0.5,
  metric = "cosine",
  force = TRUE
)

saveArchRProject(CombinedEpiProj)

#######################################
# Plot combined epithelial UMAPs
#######################################

p1 <- plotEmbedding(
  ArchRProj = CombinedEpiProj,
  colorBy = "cellColData",
  name = "Sample",
  embedding = "UMAP_Harmony"
)

p2 <- plotEmbedding(
  ArchRProj = CombinedEpiProj,
  colorBy = "GeneScoreMatrix",
  name = c("EPCAM", "KRT8", "KRT18", "MUC2"),
  embedding = "UMAP_Harmony"
)

p3 <- plotEmbedding(
  ArchRProj = CombinedEpiProj,
  colorBy = "GeneScoreMatrix",
  name = c("CD3E", "CD14", "COL1A1"),
  embedding = "UMAP_Harmony"
)

plotPDF(p1, p2, p3,
  name = "Epithelial_Combined_QC",
  ArchRProj = CombinedEpiProj,
  addDOC = FALSE,
  width = 5, height = 5
)

saveArchRProject(CombinedEpiProj)
#CombinedEpiProj <- loadArchRProject("/mnt/tier2/project/p201120/ryan/cCRE_pipeline/outputs/ArrowFiles/ArchR_Epithelial_Combined/")
#######################################
# Call peaks using MACS2
#######################################

CombinedEpiProj <- addGroupCoverages(
  ArchRProj = CombinedEpiProj,
  groupBy = "Sample",
  force = TRUE
)

pathToMacs2 <- "/mnt/tier2/project/p201120/ryan/envs/envs/macs2_env/bin/macs2"
CombinedEpiProj <- addReproduciblePeakSet(
  ArchRProj = CombinedEpiProj,
  groupBy = "Sample",
  pathToMacs2 = pathToMacs2,
  force = TRUE
)

# Add peak matrix
CombinedEpiProj <- addPeakMatrix(CombinedEpiProj)

#######################################
# Convert .rds to bed
#######################################

rds_dir <- "/mnt/tier2/project/p201120/ryan/cCRE_pipeline/outputs/ArrowFiles/ArchR_Epithelial_Combined/PeakCalls/Sample/"
out_dir  <- "/mnt/tier2/project/p201120/ryan/cCRE_pipeline/outputs/ArrowFiles/ArchR_Epithelial_Combined/PeakCalls/BED/"

rds_files <- list.files(rds_dir, pattern = "\\.gr\\.rds$", full.names = TRUE)

for (f in rds_files) {
  gr <- readRDS(f)
  sample_name <- gsub("-reproduciblePeaks\\.gr\\.rds$", "", basename(f))
  export(gr, con = file.path(out_dir, paste0(sample_name, ".bed")), format = "BED")
}







macs2_peaks <- getPeakSet(proj_colon_harmony_macs2)
macs2_peaks

# export peaks to directory
rtracklayer::export(macs2_peaks, "/Users/ryanhagan/NoCoSMiCC/ArchR_analysis/MACS2_peaks/Sigmoid_ENCSR388NCA.bed", "bed")


saveArchRProject(CombinedEpiProj)



# Create a unified peak set
projArchR <- addReproduciblePeakSet(
  ArchRProj = projArchR,
  groupBy = "Cluster", 
  pathToMacs2 = pathToMacs2
)

#######################################
# Generate sample BigWigs
#######################################

getGroupBW(
  ArchRProj = CombinedEpiProj,
  groupBy = "Sample",
  normMethod = "ReadsInTSS",
  tileSize = 25,
  maxCells = 5000
)

getOutputDirectory(ArchRProj = CombinedEpiProj)



#######################################
# Non-epithelial peaks for comparison
#######################################

# Define non-epithelial clusters
ENCODE_Nonepi_Clusters <-  list(
  "ENCSR997YNO" = c("C1", "C2"),
  "ENCSR830FPR" = c("C1", "C2", "C3", "C4"),
  "ENCSR349XKD" = c("C1", "C2"),
  "ENCSR434SXE" = c("C1"),
  "ENCSR506YMX" = c("C2"),
  "ENCSR904WIW" = c("C1", "C2", "C3", "C4", "C5","C6"),
  # Samples with no epithelial clusters are omitted
)

NonEpithelialCells <- c()

for(samp in names(ENCODE_Nonepi_Clusters)) {
  
  cat("\nSubsetting non-epithelial cells from:", samp, "\n")

  sampProj <- loadArchRProject(paste0("ArchR_", samp))
  
  NonepiClusters <- ENCODE_Nonepi_Clusters[[samp]]
  NonepiCells <- getCellNames(sampProj)[sampProj$Clusters %in% epiClusters]
  
  cat("Non-Epithelial cells:", length(NonepiCells), "\n")
  NonEpithelialCells <- c(NonEpithelialCells, NonepiCells)
}

cat("\nTotal non-epithelial cells across ENCODE samples:", length(NonEpithelialCells), "\n")







































# Iterative LSI dimensionality reduction
ENCODEproj <- addIterativeLSI(
  ArchRProj = ENCODEproj,
  useMatrix = "TileMatrix",
  name = "IterativeLSI",
  iterations = 3,
  sampleCellsPre = 25000,
  dimsToUse = 1:25,
  varFeatures = 15000,
  force = TRUE
)

# Harmony batch correction across donors
ENCODEproj <- addHarmony(
  ArchRProj = ENCODEproj,
  reducedDims = "IterativeLSI",
  name = "Harmony",
  groupBy = "Sample",
  force = TRUE
)

# Clustering
ENCODEproj <- addClusters(
  input = ENCODEproj,
  reducedDims = "Harmony",
  name = "Clusters",
  resolution = 1.5,
  nOutlier = 20,
  seed = 1,
  sampleCells = 40000,
  maxClusters = 40,
  force = TRUE
)

# UMAP visualisation
ENCODEproj <- addUMAP(
  ArchRProj = ENCODEproj,
  reducedDims = "Harmony",
  name = "UMAP_Harmony",
  nNeighbors = 30,
  minDist = 0.5,
  metric = "cosine",
  force = TRUE
)

# Save the project
saveArchRProject(ArchRProj = ENCODEproj)

# Plot clusters
p1 <- plotEmbedding(
  ArchRProj = ENCODEproj,
  colorBy = "cellColData",
  name = "Clusters",
  embedding = "UMAP_Harmony"
)

# Plot epithelial markers
p2 <- plotEmbedding(
  ArchRProj = ENCODEproj,
  colorBy = "GeneScoreMatrix",
  name = c("EPCAM", "MUC2"),
  embedding = "UMAP_Harmony",
  imputeWeights = getImputeWeights(ENCODEproj)
)

# Plot immune markers
p3 <- plotEmbedding(
  ArchRProj = ENCODEproj,
  colorBy = "GeneScoreMatrix",
  name = c("CD3E", "CD14", "CD8A"),
  embedding = "UMAP_Harmony",
  imputeWeights = getImputeWeights(ENCODEproj)
)

# Plot by sample to check Harmony integration
p4 <- plotEmbedding(
  ArchRProj = ENCODEproj,
  colorBy = "cellColData",
  name = "Sample",
  embedding = "UMAP_Harmony"
)

plotPDF(p1, p2, p3, p4,
  name = "UMAP_Markers_Clusters.pdf",
  ArchRProj = ENCODEproj,
  addDOC = FALSE,
  width = 5, height = 5
)






































#####################################
## Iterative LSI dimensionality reduction
ENCODEproj <- addIterativeLSI(
  ArchRProj = ENCODEproj,
  useMatrix = "TileMatrix", 
  name = "IterativeLSI", 
  iterations = 2, 
  clusterParams = list( #See Seurat::FindClusters
    resolution = c(0.2), 
    sampleCells = 10000, 
    n.start = 10
  ), 
  varFeatures = 15000, 
  dimsToUse = 1:15,
  force = TRUE
)
ENCODEproj

# Clustering
ENCODEproj <- addClusters(
  input = ENCODEproj,
  reducedDims = "IterativeLSI",
  method = "Seurat",
  name = "Clusters",
  resolution = 0.5,
  force = TRUE
)
ENCODEproj

table(ENCODEproj$Clusters)

# UMAP projection
ENCODEproj <- addUMAP(
  ArchRProj = ENCODEproj, 
  reducedDims = "IterativeLSI", 
  name = "UMAP", 
  nNeighbors = 30, 
  minDist = 0.5, 
  metric = "cosine",
  force = TRUE
)

# Plot UMAP by sample and clusters
sample_umap <- plotEmbedding(ArchRProj = ENCODEproj, 
                    colorBy = "cellColData", 
                    name = "Sample", embedding = "UMAP")

cluster_umap <- plotEmbedding(ArchRProj = ENCODEproj, 
                    colorBy = "cellColData", 
                    name = "Clusters", embedding = "UMAP")

ggAlignPlots(sample_umap, cluster_umap, type = "h")






# Subset to epithelial cells using published barcodes






saveArchRProject(
  ArchRProj       = proj,
  outputDirectory = "/project/home/p201120/ryan/cCRE_pipeline/files/scATAC/ArchR_output",
  load            = FALSE,
  overwrite       = FALSE
)







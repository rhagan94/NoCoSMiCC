#### Script 1 - Process-sc-fragments.R ####

## This script enables the identification of epithelial cells from ENCODE colon scATAC-seq datasets

## Author: Dr Ryan Hagan, RCSI

## Date: April 2026

# Load libraries
library(ArchR)
library(stringr)
library(Signac)

# set the WD
setwd("/Users/ryanhagan/NoCoSMiCC/ArchR_analysis")

# Set threads to half total number of cores available
addArchRThreads(threads = 3) 

# Set the genome to human
addArchRGenome("hg38")

# Obtain the list of input files 
inputFiles <- list.files(path = "/Users/ryanhagan/NoCoSMiCC/raw_data/scatac_data/fragment_files",
                         pattern = ".tsv.gz", full.names = TRUE)

# Set the input file names
names(inputFiles) <- str_sub(inputFiles, start = 64, end = -18)
names(inputFiles)

# Create the arrow files
ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,
  sampleNames = names(inputFiles),
  filterTSS = 8, # adjust as necessary
  filterFrags = 1000,
  addTileMat = TRUE,
  addGeneScoreMat = FALSE,
  subThreading = FALSE,
  force = TRUE
)
ArrowFiles

# Set the ArrowFiles to a single sample to process individually
# ArrowFiles <- "/Users/ryanhagan/NoCosMiCC/ArchR_analysis/HBM425.arrow"

# Compute doublet scores
doubScores <- addDoubletScores(
  input = ArrowFiles,
  k = 10,
  knnMethod = "UMAP",
  LSIMethod = 1
)

# Create the ArchR project
projColon1 <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  outputDirectory = "/Users/ryanhagan/NoCosMiCC/ArchR_analysis",
  copyArrows = FALSE
)
projColon1

# Plot quality control metrics
p1 <- plotGroups(
  ArchRProj = projColon1, 
  groupBy = "Sample", 
  colorBy = "cellColData", 
  name = "TSSEnrichment",
  plotAs = "violin",
  alpha = 0.8,
  addBoxPlot = TRUE
)
p1

p2 <- plotGroups(
  ArchRProj = projColon1, 
  groupBy = "Sample", 
  colorBy = "cellColData", 
  name = "log10(nFrags)",
  plotAs = "ridges",
  pal = c("#272E6A")
)
p2

# Filter doublets
proj_Colon_Filt <- filterDoublets(projColon1)
proj_Colon_Filt

# Iterative LSI dimensionality reduction
proj_Colon_Filt <- addIterativeLSI(
  ArchRProj = proj_Colon_Filt,
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
proj_Colon_Filt

# Clustering
proj_Colon_Filt <- addClusters(
  input = proj_Colon_Filt,
  reducedDims = "IterativeLSI",
  method = "Seurat",
  name = "Clusters",
  resolution = 0.4,
  force = TRUE
)
proj_Colon_Filt

table(proj_Colon_Filt$Clusters)

# UMAP projection
proj_Colon_Filt <- addUMAP(
  ArchRProj = proj_Colon_Filt, 
  reducedDims = "IterativeLSI", 
  name = "UMAP", 
  nNeighbors = 30, 
  minDist = 0.5, 
  metric = "cosine",
  force = TRUE
)

# Plot UMAP by sample and clusters
sample_umap <- plotEmbedding(ArchRProj = proj_Colon_Filt, 
                    colorBy = "cellColData", 
                    name = "Sample", embedding = "UMAP")

cluster_umap <- plotEmbedding(ArchRProj = proj_Colon_Filt, 
                    colorBy = "cellColData", 
                    name = "Clusters", embedding = "UMAP")

ggAlignPlots(sample_umap, cluster_umap, type = "h")

# Add group coverages to the ArchR project based on Clusters
#proj_Colon_Filt <- addGroupCoverages(ArchRProj = proj_Colon_Filt, groupBy = "Clusters")

# Assess coverage using the ArchR browser - can be used to assess accessibility around marker genes (e.g. EPCAM and CDH1 for epithelial cells)
ArchRBrowser(proj_Colon_Filt)

# Add the gene score matrix
proj_Colon_Filt <- addGeneScoreMatrix(proj_Colon_Filt)

# Add impute weights
proj_Colon_Filt <- addImputeWeights(proj_Colon_Filt)

# Obtain marker genes for each cluster
markersGS <- getMarkerFeatures(
    ArchRProj = proj_Colon_Filt, 
    useMatrix = "GeneScoreMatrix", 
    groupBy = "Clusters",
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon"
)
markersGS

markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.01 & Log2FC >= 0.5")

#write.csv(markerList$C1$name, "cluster1_markers.csv")

# overlap DA genes with marker gene list(s)

epithelial_markers <- c(
#"MUC2", "MUC5AC", "CLCA1", "MUC13", "FCGBP", "ZG16", "ATOH1", "SPDEF", "REP15" # Goblet cells
"CDH1", "EPCAM", "KRT14", "MUC1", "CD24", "CEACAM1", "KRT3", "ANPEP" # general epithelial cells
) 
epithelial_markers

#intersect(x = markerList$C7$name, y = epithelial_markers)

p <- plotEmbedding(
    ArchRProj = proj_Colon_Filt, 
    colorBy = "GeneScoreMatrix", 
    name = epithelial_markers, 
    embedding = "UMAP",
    quantCut = c(0.01, 0.95),
    imputeWeights = getImputeWeights(proj_Colon_Filt)
)
p

p2 <- lapply(p, function(x){
    x + guides(color = FALSE, fill = FALSE) + 
    theme_ArchR(baseSize = 6.5) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    theme(
        axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(), 
        axis.text.y=element_blank(), 
        axis.ticks.y=element_blank()
    )
})
do.call(cowplot::plot_grid, c(list(ncol = 3),p2))

# Subset the project to retain only Epithelial cells 

epithelial_idx <- BiocGenerics::which(proj_Colon_Filt$Clusters %in% c("C1", "C2", "C3"))
epithelial_cells <- proj_Colon_Filt$cellNames[epithelial_idx]
proj_epithelial <- proj_Colon_Filt[epithelial_cells, ]
proj_epithelial

# Create a separate subset of the non-epithelial cells
# non_epithelial_idx <- BiocGenerics::which(proj_Colon_Filt$Clusters %in% c("C5", "C6", "C7", "C8")
# non_epithelial_cells <- proj_Colon_Filt$cellNames[non_epithelial_idx]
# proj_non_epithelial <- proj_Colon_Filt[non_epithelial_cells, ]
# proj_non_epithelial

# Re-run the dimensionality reduction and clustering steps
# Iterative LSI dimensionality reduction
proj_epithelial <- addIterativeLSI(
  ArchRProj = proj_epithelial,
  useMatrix = "TileMatrix", 
  name = "IterativeLSI", 
  iterations = 2, 
  clusterParams = list(
    resolution = c(0.2), 
    sampleCells = 5000, 
    n.start = 10
  ), 
  varFeatures = 25000, 
  dimsToUse = 1:30,
  force = TRUE
)
proj_epithelial

# clustering
proj_epithelial <- addClusters(
  input = proj_epithelial,
  reducedDims = "IterativeLSI",
  method = "Seurat",
  name = "Clusters",
  resolution = 0.8,
  force = TRUE
)

# UMAP projection
proj_epithelial <- addUMAP(
  ArchRProj = proj_epithelial, 
  reducedDims = "IterativeLSI", 
  name = "UMAP", 
  nNeighbors = 30, 
  minDist = 0.5, 
  metric = "cosine",
  force = TRUE
)
proj_epithelial

umap <- plotEmbedding(ArchRProj = proj_epithelial, 
                    colorBy = "cellColData", 
                    name = "Clusters", embedding = "UMAP", baseSize = 8)
umap

### Signac analysis ###

## Sample 1 - example shown - repeat for all samples

# Obtain the cell barcodes from ArchR project
archr_barcodes <- getCellNames(ArchRProj = proj_epithelial)
name <- stringr::str_sub(archr_barcodes[1], start = 1, end=7)
barcodes <- gsub(name, "", archr_barcodes)

# Create the filtered fragment file
fpath <- "/Users/ryanhagan/NoCoSMiCC/raw_data/scatac_data/fragment_files/HuBMAP/HBM338.SSPB.265_frags.sort.bed.gz"
cells <- barcodes
outfile <- "/Users/ryanhagan/NoCoSMiCC/raw_data/scatac_data/fragment_files/HuBMAP/filtered_HBM338.SSPB.265_frags.sort.bed.gz"
FilterCells(
  fragments = fpath,
  cells = cells,
  outfile = outfile
)

# Define Macs2 path
macs2="/Users/ryanhagan/miniconda3/envs/macs2_env/bin/macs2"

# call peaks on the filtered fragment list
Filtered_HBM338_peaks <- CallPeaks(object = outfile,
    macs2.path = macs2)
Filtered_HBM338_peaks

# Plot peak widths
ggplot(as.data.frame(Filtered_HBM425_peaks@ranges@width), aes(x = Filtered_HBM425_peaks@ranges@width)) + 
  geom_histogram(color="#3B9AB2", fill="#3B9AB2", bins = 100, alpha = 0.5) +
  xlim(0,3000) +
  labs(x = "Width (bp)", y = "Frequency")+
  geom_vline(aes(xintercept = mean(Filtered_HBM425_peaks@ranges@width)),col='coral',size=0.8, linetype = "dashed")+
  theme_bw() +
  theme(axis.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"),
        plot.title = element_text(face = "bold")) +
  ggtitle("Epithelial peaks")

# Export the peaks
rtracklayer::export(Filtered_HBM338_peaks, "/Users/ryanhagan/NoCoSMiCC/ArchR_analysis/MACS2_peaks/HBM338_epithelial.bed", "bed")

# Create a BigWig file 
getGroupBW(proj_epithelial,
  groupBy = "Sample",
  maxCells = 10000)

# Save the project
saveArchRProject(ArchRProj = proj_epithelial, outputDirectory = "/Users/ryanhagan/NoCosMiCC/ArchR_analysis", load = FALSE)

# To load the ArchR project for later use:
# loadArchRProject(path = "/Users/ryanhagan/NoCosMiCC/ArchR_analysis")

# END OF SCRIPT




## New function for calling epithelial cells
# ============================================================
# Epithelial Cluster Identification and Extraction Pipeline
# ============================================================

# Define marker lists

epithelialMarkers <- c(
  # General epithelial
  "EPCAM", "KRT8", "KRT18", "CDH1",
  # Colonocyte/absorptive
  "CA1", "SLC26A2", "SLC26A3", "RAB6B", "FABP1",
  # Best4+ enterocytes
  "BEST4", "CA7", "OTOP2", "OTOP3",
  # Goblet
  "MUC2", "TFF1", "FCGBP", "FFAR4", "CLCA1", "RETNLB", "SPINK4", "AGR2",
  # Immature goblet
  "KLK1", "ITLN1", "WFDC2", "LRRC26",
  # Stem
  "SMOC2", "LGR5", "ASCL2", "SOX9", "RGMB",
  # Cycling TA
  "TICRR", "CDC25C",
  # Enteroendocrine
  "CHGA", "NEUROD1", "SCGN", "FEV", "PYY",
  # Tuft
  "TRPM5", "GNG13", "DCLK1", "LRMP",
  # Paneth
  "DEFA5", "DEFA6", "LYZ", "REG3A",
  # M cells
  "GP2", "TNFRSF11A"
)

nonEpithelialMarkers <- c(
  # Pan-immune
  "PTPRC",
  # T cells
  "CD3D", "CD3E", "CD3G", "CD8A", "CD8B", "CD4", "CD2", "IL7R",
  # Tregs
  "FOXP3", "CTLA4", "BATF",
  # Exhausted T
  "PDCD1", "HAVCR2", "LAG3", "TOX",
  # Th17
  "IL17A", "KLRB1",
  # NK
  "KLRF1", "NCAM1", "SH2D1B",
  # B cells
  "PAX5", "MS4A1", "CD19",
  # Plasma
  "IGLL1", "MZB1", "JCHAIN", "SSR4",
  # GC B cells
  "FCRLA", "CD180",
  # Mast
  "TPSAB1", "CMA1", "GATA2", "HDC",
  # Monocyte/macrophage
  "CD14", "CD68", "FCGR1A", "LILRB2", "FOLR2",
  # DC
  "CLEC9A", "ITGAX", "LY75",
  # Cycling immune
  "MKI67", "TOP2A",
  # Fibroblasts
  "COL1A1", "COL1A2", "COL6A1", "FAP", "DCN", "CBLN2",
  # Inflammatory fibroblasts
  "MMP3", "MMP1", "CHI3L1",
  # Myofibroblasts
  "MYH11", "TAGLN", "ACTA2", "DES",
  # Endothelial
  "VWF", "CDH5", "PLVAP", "PECAM1",
  # Post-capillary venules
  "SELP", "MADCAM1",
  # Pericytes
  "RGS5", "MCAM", "NOTCH3", "KCNJ8",
  # Schwann/neural
  "S100B", "SOX10", "MBP", "MAP2",
  # Interstitial cells of Cajal
  "ANO1"
)

# Define function for epithelial scoring

scoreClustersByMarkers <- function(ArchRProj, groupBy = "Clusters",
                                   positiveMarkers, negativeMarkers) {
  
  markerGS <- getMarkerFeatures(
    ArchRProj = ArchRProj,
    useMatrix = "GeneScoreMatrix",
    groupBy = groupBy,
    bias = c("TSSEnrichment", "log10(nFrags)"),
    testMethod = "wilcoxon"
  )
  
  log2fc <- assay(markerGS, "Log2FC")
  rownames(log2fc) <- rowData(markerGS)$name
  
  # Report how many markers were found
  cat("Epithelial markers found:", sum(positiveMarkers %in% rownames(log2fc)),
      "of", length(positiveMarkers), "\n")
  cat("Non-epithelial markers found:", sum(negativeMarkers %in% rownames(log2fc)),
      "of", length(negativeMarkers), "\n")
  
  clusterScores <- sapply(colnames(log2fc), function(cluster) {
    
    posFound <- positiveMarkers[positiveMarkers %in% rownames(log2fc)]
    negFound <- negativeMarkers[negativeMarkers %in% rownames(log2fc)]
    
    posScore <- mean(log2fc[posFound, cluster], na.rm = TRUE)
    negScore <- mean(log2fc[negFound, cluster], na.rm = TRUE)
    
    posScore - negScore
  })
  
  return(sort(clusterScores, decreasing = TRUE))
}

# Run scoring

# Score clusters
clusterScores <- scoreClustersByMarkers(
  ArchRProj = proj_Colon_Filt,
  groupBy = "Clusters",
  positiveMarkers = epithelialMarkers,
  negativeMarkers = nonEpithelialMarkers
)

cat("\nCluster scores:\n")
print(clusterScores)

# Find epithelial clusters
epithelialClusters <- names(clusterScores[clusterScores > 0])

# Subset to epithelial clusters
cat("\nSubsetting to epithelial clusters:", epithelialClusters, "\n")

proj_epithelial <- subsetArchRProject(
  ArchRProj = proj_Colon_Filt,
  cells = getCellNames(proj_Colon_Filt)[proj_Colon_Filt$Clusters %in% epithelialClusters],
  outputDirectory = "ArchR_epithelial",
  dropCells = TRUE,
  force = TRUE
)

cat("\nEpithelial project summary:\n")
print(proj_epithelial)




#############################################################################################
##########################             2026 Update                  #########################
#############################################################################################

# Script for ArchR project creation using ENCODE and HuBMAP scATAC datasets
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
setwd("/project/home/p201120/ryan/cCRE_pipeline/outputs/ArrowFiles/")

# Define arrow/frag file directories
arrowdir <- getwd()
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








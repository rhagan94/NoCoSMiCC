#### Script 1 - Process-sc-fragments.R ####

## This script enables the generation of Chromatin Accessible Regions (CARs) ((aka peaks)) from scATAC-seq samples

## Author: Dr Ryan Hagan, RCSI

## Date: July 2024

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

############################################################################################################################
#..........................................................................................................................#
############################################################################################################################
# Load packages
library(ArchR)
library(Seurat)
library(BSgenome.Hsapiens.UCSC.hg38)
library(parallel)
'%notin%' <- Negate('%in%')

# Load Genome Annotations
addArchRGenome("hg38")

# Set Threads to be used
addArchRThreads(threads = 96)

# Set/Create Working Directory to Folder
setwd("/project/home/p201120/ryan/cCRE_pipeline/outputs/ArrowFiles")

fragDir <- "/project/home/p201120/ryan/cCRE_pipeline/files/scATAC/"

############################################################################################################################
#..........................................................................................................................#
############################################################################################################################
# Paths to ENCODE and HuBMAP fragment files

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

############################################################################################################################
#..........................................................................................................................#
############################################################################################################################
# Create Arrow files

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

# Reassign as character vectors of arrow file paths
ArrowFiles1K   <- paste0(names(inputFiles1K),   ".arrow")
ArrowFiles1p5K <- paste0(names(inputFiles1p5K), ".arrow")
ArrowFiles2K   <- paste0(names(inputFiles2K),   ".arrow")
ArrowFiles3K   <- paste0(names(inputFiles3K),   ".arrow")
ArrowFiles4K   <- paste0(names(inputFiles4K),   ".arrow")
ArrowFiles5K   <- paste0(names(inputFiles5K),   ".arrow")
ArrowFiles6K   <- paste0(names(inputFiles6K),   ".arrow")

# ENCODE samples
ENCODE_Arrows <- c("ENCSR916RYB.arrow", "ENCSR904WIW.arrow", "ENCSR830FPR.arrow",
                   "ENCSR997YNO.arrow", "ENCSR007QIO.arrow", "ENCSR349XKD.arrow", "ENCSR434SXE.arrow",
                   "ENCSR388NCA.arrow", "ENCSR506YMX.arrow", "ENCSR367GKP.arrow")

doubletScores <- addDoubletScores(ENCODE_Arrows,   
                                  k = 10, knnMethod = "UMAP", 
                                  LSIMethod = 1, threads = 1)
                   
############################################################################################################################
#..........................................................................................................................#
############################################################################################################################
# Compute doublet scores for ENCODE samples

doubScores <- addDoubletScores(ArrowFiles1K,   k = 10, knnMethod = "UMAP", LSIMethod = 1, threads = 1)
doubScores <- addDoubletScores(ArrowFiles1p5K, k = 10, knnMethod = "UMAP", LSIMethod = 1, threads = 1)
doubScores <- addDoubletScores(ArrowFiles2K,   k = 10, knnMethod = "UMAP", LSIMethod = 1, threads = 1)
doubScores <- addDoubletScores(ArrowFiles3K,   k = 10, knnMethod = "UMAP", LSIMethod = 1, threads = 1)
doubScores <- addDoubletScores(ArrowFiles4K,   k = 10, knnMethod = "UMAP", LSIMethod = 1, threads = 1)
doubScores <- addDoubletScores(ArrowFiles5K,   k = 10, knnMethod = "UMAP", LSIMethod = 1, threads = 1)
doubScores <- addDoubletScores(ArrowFiles6K,   k = 10, knnMethod = "UMAP", LSIMethod = 1, threads = 1)

############################################################################################################################
#..........................................................................................................................#
############################################################################################################################
# Create ArchR project

set.seed(1)
allArrowFiles <- c(ArrowFiles1K, ArrowFiles1p5K, ArrowFiles2K, ArrowFiles3K, ArrowFiles4K, ArrowFiles5K, ArrowFiles6K)

proj <- ArchRProject(
  ArrowFiles        = allArrowFiles,
  outputDirectory   = "/project/home/p201120/ryan/cCRE_pipeline/files/scATAC/ArchR_output"
)

saveArchRProject(
  ArchRProj       = proj,
  outputDirectory = "/project/home/p201120/ryan/cCRE_pipeline/files/scATAC/ArchR_output",
  load            = FALSE,
  overwrite       = FALSE
)

proj <- loadArchRProject("/project/home/p201120/ryan/cCRE_pipeline/files/scATAC/ArchR_output")
proj <- filterDoublets(proj, filterRatio = 1.2)
write.table(rownames(getCellColData(proj)), "initial_post_filter_atac_cells.txt")






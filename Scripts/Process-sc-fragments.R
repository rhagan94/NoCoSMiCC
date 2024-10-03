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
p1 <- plotEmbedding(ArchRProj = proj_Colon_Filt, 
                    colorBy = "cellColData", 
                    name = "Sample", embedding = "UMAP")

p2 <- plotEmbedding(ArchRProj = proj_Colon_Filt, 
                    colorBy = "cellColData", 
                    name = "Clusters", embedding = "UMAP")

ggAlignPlots(p1, p2, type = "h")

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

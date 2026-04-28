#### Synteny analysis using GENESPACE ####

# Load libraries
library(GENESPACE)
library(dplyr)
library(GenomicRanges)
library(rtracklayer)

# Get data together
genomeRepo <- "/project/home/p201120/ryan/Synteny/Genomes"
wd <- "/project/home/p201120/ryan/Synteny"
path2mcscanx <- "/mnt/tier2/project/p201120/ryan/envs/get/bin"

# Download the Human, Mouse and Dog genome data
urls <- c(
  human ="000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_",
  mouse = "000/001/635/GCF_000001635.27_GRCm39/GCF_000001635.27_GRCm39_",
  dog = "000/002/285/GCF_000002285.5_Dog10K_Boxer_Tasha/GCF_000002285.5_Dog10K_Boxer_Tasha_")

genomes2run <- names(urls)
urls <- file.path("https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF", urls)
translatedCDS <- sprintf("%stranslated_cds.faa.gz", urls)
geneGff <- sprintf("%sgenomic.gff.gz", urls)

names(translatedCDS) <- genomes2run
names(geneGff) <- genomes2run
writeDirs <- file.path(genomeRepo, genomes2run)
names(writeDirs) <- genomes2run
for(i in genomes2run){
  print(i)
  if(!dir.exists(writeDirs[i]))
    dir.create(writeDirs[i])
  download.file(
    url = geneGff[i], 
    destfile = file.path(writeDirs[i], basename(geneGff[i])))
  download.file(
    url = translatedCDS[i], 
    destfile = file.path(writeDirs[i], basename(translatedCDS[i])))
}

# Parse the annotations
genomes2run <- c("human", "mouse", "dog")
parsedPaths <- parse_annotations(
  rawGenomeRepo = genomeRepo,
  genomeDirs = genomes2run,
  genomeIDs = genomes2run,
  presets = "ncbi",
  genespaceWd = wd)

# Initialize the GENESPACE run and set parameters
gpar <- init_genespace(
  wd = wd,
  path2mcscanx = path2mcscanx)

# Conduct GENESPACE run
out <- run_genespace(gpar, overwrite = T)

# Read the syntenic block coordinates into R
syntenic_blocks <- read.csv(file.path(wd, "results/syntenicBlock_coordinates.csv"))

syntenic_blocks_human_mouse <- as.data.frame(syntenic_blocks) %>% 
                        filter(genome1 == 'human', genome2 == 'mouse'| (genome1 == 'mouse' & genome2 == 'human')) %>%
                        dplyr::select(chr1, startBp1, endBp1) %>%
                        dplyr::rename(Chr = chr1,
                               Start = startBp1,
                               End = endBp1)
syntenic_blocks_human_mouse     

syntenic_blocks_human_mouse$Chr <- paste("chr", syntenic_blocks_human_mouse$Chr, sep = "")

syntenic_blocks_human_mouse_GR <- makeGRangesFromDataFrame(syntenic_blocks_human_mouse)
syntenic_blocks_human_mouse_GR

syntenic_blocks_human_dog <- as.data.frame(syntenic_blocks) %>% 
  filter(genome1 == 'human', genome2 == 'dog'| (genome1 == 'dog' & genome2 == 'human')) %>%
  dplyr::select(chr1, startBp1, endBp1) %>%
  dplyr::rename(Chr = chr1,
                Start = startBp1,
                End = endBp1)
syntenic_blocks_human_dog     

syntenic_blocks_human_dog$Chr <- paste("chr", syntenic_blocks_human_dog$Chr, sep = "")

syntenic_blocks_human_dog_GR <- makeGRangesFromDataFrame(syntenic_blocks_human_dog)
syntenic_blocks_human_dog_GR

syntenic_blocks_both_GR <- subsetByOverlaps(syntenic_blocks_human_mouse_GR,
                                             syntenic_blocks_human_dog_GR)

cat("Human-mouse syntenic blocks: ", length(syntenic_blocks_human_mouse_GR), "\n")
cat("Human-dog syntenic blocks: ", length(syntenic_blocks_human_dog_GR), "\n")
cat("Conserved in both mouse and dog: ", length(syntenic_blocks_both_GR), "\n")

# ===========================================================================
# 7. Overlap with colon epithelial cCREs
# ===========================================================================

ccre_gr <- import(
  "/mnt/tier2/project/p201120/ryan/cCRE_pipeline/outputs/colon-epithelial-cCREs.bed",
  format    = "BED",
  extraCols = c(type = "character")
)

# cCREs within syntenic regions
ccres_in_syntenic <- subsetByOverlaps(ccre_gr, syntenic_blocks_both_GR)

# cCREs outside syntenic regions (potentially human-specific)
ccres_not_in_syntenic <- subsetByOverlaps(ccre_gr, syntenic_blocks_both_GR, invert = TRUE)

cat("Total cCREs: ", length(ccre_gr), "\n")
cat("cCREs in syntenic regions: ", length(ccres_in_syntenic), "\n")
cat("cCREs NOT in syntenic regions: ", length(ccres_not_in_syntenic), "\n")
cat("% in syntenic regions:", round(length(ccres_in_syntenic) / length(ccre_gr) * 100, 1), "%\n")

export(ccres_in_syntenic,
  file.path(wd, "results/cCREs_in_syntenic_regions.bed"),
  format = "BED")

export(ccres_not_in_syntenic,
  file.path(wd, "results/cCREs_NOT_in_syntenic_regions.bed"),
  format = "BED")












############################################################################################


# assess peak widths 
blocks_width <- as.data.frame(syntenic_blocks_both_GR@ranges@width)
blocks_width

library(ggplot2)
p <- ggplot(blocks_width, aes(x = syntenic_blocks_both_GR@ranges@width)) + 
  geom_histogram(color="#bd7ebe", fill="#bd7ebe", bins = 100, alpha = 0.5) +
 # xlim(0,2000) +
  labs(x = "Width (bp)", y = "Frequency")+
  geom_vline(aes(xintercept = median(syntenic_blocks_human_GR@ranges@width)),col='#b2e061',size=0.8, linetype = "dashed")+
  theme_bw() +
  theme(axis.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"),
        plot.title = element_text(face = "bold")) +
  ggtitle("Human-mouse-dog syntenic regions")

ggsave(
  filename = "/mnt/tier2/project/p201120/ryan/Synteny/results/syntenic_regions.png",
  plot = p,
  width = 8,
  height = 6,
  units = "in",
  dpi = 300,
  bg = "white"
)

# export synteny blocks as a bed file
rtracklayer::export.bed(syntenic_blocks_human_GR, "/Users/ryanhagan/NoCoSMiCC/Synteny/human_dog_syntenic_blocks.bed")

# Dog
# assess peak widths 
blocks_width <- as.data.frame(syntenic_blocks_human_dog_GR@ranges@width)
blocks_width

library(ggplot2)
ggplot(blocks_width, aes(x = syntenic_blocks_human_dog_GR@ranges@width)) + 
  geom_histogram(color="#bd7ebe", fill="#bd7ebe", bins = 100, alpha = 0.5) +
  # xlim(0,2000) +
  labs(x = "Width (bp)", y = "Frequency")+
  geom_vline(aes(xintercept = mean(syntenic_blocks_human_dog_GR@ranges@width)),col='#b2e061',size=0.8, linetype = "dashed")+
  theme_bw() +
  theme(axis.title.x = element_text(face = "bold"),
        axis.title.y = element_text(face = "bold"),
        plot.title = element_text(face = "bold")) +
  ggtitle("Human-dog syntenic regions")

# export synteny blocks as a bed file
rtracklayer::export.bed(syntenic_blocks_human_dog_GR, "/Users/ryanhagan/NoCoSMiCC/Synteny/human_dog_syntenic_blocks.bed")






# Compare syntenic and non-syntenic cCREs
all_cCREs <- ChIPQC:::GetGRanges("/project/home/p201120/ryan/cCRE_pipeline/outputs/colon-epithelial-cCREs.bed", simple = TRUE)
all_cCREs

syntenic_cCREs <- ChIPQC:::GetGRanges("results/cCREs_in_syntenic_regions.bed", simple = TRUE)
syntenic_cCREs

nonsyntenic_cCREs <- ChIPQC:::GetGRanges("results/cCREs_NOT_in_syntenic_regions.bed", simple = TRUE)
nonsyntenic_cCREs

# plot genomic locations
atac.annotation = annotatePeak(all_cCREs, tssRegion=c(-3000, 3000),TxDb=txdb, annoDb="org.Hs.eg.db")
atac.annotation
plotAnnoPie(atac.annotation,  col = c("#fd7f6f", "#7eb0d5", "#b2e061", "#bd7ebe", "#ffb55a", "#ffee65", "#beb9db", "#fdcce5", "#8bd3c7", "#d7e3fc", "#ffc09f"))

atac.annotation = annotatePeak(syntenic_cCREs, tssRegion=c(-3000, 3000),TxDb=txdb, annoDb="org.Hs.eg.db")
atac.annotation
plotAnnoPie(atac.annotation,  col = c("#fd7f6f", "#7eb0d5", "#b2e061", "#bd7ebe", "#ffb55a", "#ffee65", "#beb9db", "#fdcce5", "#8bd3c7", "#d7e3fc", "#ffc09f"))

atac.annotation = annotatePeak(nonsyntenic_cCREs, tssRegion=c(-3000, 3000),TxDb=txdb, annoDb="org.Hs.eg.db")
atac.annotation
plotAnnoPie(atac.annotation,  col = c("#fd7f6f", "#7eb0d5", "#b2e061", "#bd7ebe", "#ffb55a", "#ffee65", "#beb9db", "#fdcce5", "#8bd3c7", "#d7e3fc", "#ffc09f"))


# Read the common syntenic regions into R
common_syntenic_regions <- rtracklayer::import.bed("/Users/ryanhagan/NoCoSMiCC/Synteny/common_syntenic_blocks.bed")
common_syntenic_regions

# Visualize

# convert cCREs into df
common_syntenic_df <- as.data.frame(common_syntenic_regions) %>% 
  dplyr::select(seqnames, start, end) %>%
  dplyr::rename(chr = seqnames)
common_syntenic_df

# add genome info
common_syntenic_df$genome <- "human"

# format Chr name
common_syntenic_df$chr <- gsub("chr", "", common_syntenic_df$chr)

# create ROI df
roi <- data.frame(
  genome = common_syntenic_df$genome,
  chr = common_syntenic_df$chr,
  start = common_syntenic_df$start,
  end = common_syntenic_df$end)

ggthemes <- ggplot2::theme(
  panel.background = ggplot2::element_rect(fill = "white"))

ripd <- plot_riparian(
  gsParam = out, 
  useRegions = FALSE, 
  refGenome = "human",
  genomeIDs = c("dog", "mouse", "human"),
  highlightBed = head(roi, n=15),
  backgroundColor = NULL,
  addThemes = ggthemes,
  chrFill = "#17B5C5",
  useOrder = FALSE,
  braidAlpha = 0.8,
  scaleGapSize = .5)




# Find unique cCREs
# Read the unlifted BED file for dog
bed_data_dog_unlifted <- read.table("/Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/LiftOver/hg38-cCREs-sc-dog-unlifted.bed", 
                                    header = FALSE, 
                                    sep = "\t", 
                                    stringsAsFactors = FALSE)

# View the first few rows to ensure it is read correctly
head(bed_data_dog_unlifted)

# Find unique values in the 4th column (Column 4)
unique_column4_values_dog_unlifted <- unique(bed_data_dog_unlifted$V4)

# Print unique values
print(unique_column4_values_dog_unlifted)



### Plot a cCRE in a lifted/syntenic region
# Load the necessary libraries
library(Gviz)
library(GenomicRanges)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)
library(org.Hs.eg.db)

# Read the BED file
bed_file <- "/Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/sc-only-rCARs/MaxZ/hg38-cCREs-sc-only-lifted.bed"
cCREs <- read.table(bed_file, sep="\t", stringsAsFactors=FALSE,
                    col.names=c("chrom", "start", "end", "name", "type"))

# Since your file doesn't explicitly list "dELS", let's assume you're looking 
# for distal enhancer-like signatures, which might be labeled differently
# For now, let's look for "CA-only" entries, as these are often enhancers
# (Adjust this criteria based on your knowledge of the data)
enhancer_indices <- grep("dELS", cCREs$type)

if(length(enhancer_indices) > 0) {
  # Take the first match
  cCRE_to_plot <- cCREs[enhancer_indices[44], ]
} else {
  # If no CA-only, just take the first entry
  cCRE_to_plot <- cCREs[44, ]
}

# Print the selected cCRE for reference
print(cCRE_to_plot)

# Define the genomic region with larger padding to see nearby genes
chr <- cCRE_to_plot$chrom
start_pos <- cCRE_to_plot$start - 20000  # 50kb padding
end_pos <- cCRE_to_plot$end + 20000      # 50kb padding

# Create GRanges object for the selected cCRE
cCRE_gr <- GRanges(
  seqnames = chr,
  ranges = IRanges(start = cCRE_to_plot$start, end = cCRE_to_plot$end),
  strand = "*",
  name = cCRE_to_plot$name,
  type = cCRE_to_plot$type
)

# Create tracks for visualization
gtrack <- GenomeAxisTrack()
itrack <- IdeogramTrack(genome = "hg38", chromosome = chr)

# Create a track for the selected cCRE
dtrack <- AnnotationTrack(
  cCRE_gr,
  name = paste0(cCRE_to_plot$name, " (", cCRE_to_plot$type, ")"),
  fill = "darkred",
  col = "black"
)

# Get gene models from UCSC
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
genes_track <- GeneRegionTrack(
  txdb,
  chromosome = chr,
  start = start_pos,
  end = end_pos,
  name = "Genes",
  showId = TRUE
)

# Add gene symbols
symbols <- mapIds(org.Hs.eg.db, 
                  keys = genes_track@range$gene, 
                  column = "SYMBOL", 
                  keytype = "ENTREZID")
genes_track@range$symbol <- symbols
displayPars(genes_track) <- list(
  transcriptAnnotation = "symbol",
  showId = TRUE
)

# Plot all tracks
plotTracks(
  list(itrack, gtrack, dtrack, genes_track),
  from = start_pos,
  to = end_pos,
  chromosome = chr,
  sizes = c(1, 2, 2, 5)
)

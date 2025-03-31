### Script for plotting the cCRE regions of interest - including syntenic blocks and genes


ccre_bed = read.table("/Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/LiftOver/sc-only-lifted-cCREs.bed")
ccre_gr <- makeGRangesFromDataFrame(ccre_bed, keep.extra.columns=TRUE,
                                    start.field="V2",
                                    end.field = "V3",
                                    seqnames.field="V1")
ccre_gr

# Load the gene annotation GTF file
gtf <- import("/Users/ryanhagan/NoCoSMiCC/TAD_Analysis/Files/gencode.v46.annotation.gtf.gz")

# Convert the GTF to GRanges format and keep only protein-coding genes
gtf_genes <- gtf[gtf$type == "gene" & gtf$gene_type == "protein_coding"]

genes_gr <- as(gtf_genes, "GRanges")

nearest_genes <- nearest(ccre_gr, genes_gr)
nearest_genes_info <- genes_gr[nearest_genes]
length(unique(nearest_genes_info$gene_id))

library(org.Hs.eg.db)

cleaned_gene_ids <- gsub("\\..*", "", nearest_genes_info$gene_id)
gene_symbols <- mapIds(org.Hs.eg.db, 
                       keys =cleaned_gene_ids, 
                       column = "SYMBOL", 
                       keytype = "ENSEMBL", 
                       multiVals = "first")

# Combine the results with cCRE data
ccre_nearest_genes <- cbind(as.data.frame(ccre_bed), gene_symbols)
write.csv(ccre_nearest_genes, "/Users/ryanhagan/NoCoSMiCC/HiC_Analysis/constrained_ccres_nearest_genes.csv")
head(ccre_nearest_genes)
unique(ccre_nearest_genes$gene_symbols)
write.csv(ccre_nearest_genes$gene_symbols, "/Users/ryanhagan/NoCoSMiCC/HiC_Analysis/constrained_ccres_nearest_genes.csv")

library(clusterProfiler)

# Perform GO enrichment analysis
ego <- enrichGO(gene = cleaned_gene_ids,
                OrgDb = org.Hs.eg.db,
                keyType = "ENSEMBL",
                ont = "BP",         # Biological Process
                pAdjustMethod = "BH",
                qvalueCutoff = 0.05)

# Visualize GO enrichment results
barplot(ego, showCategory = 15, order = TRUE)
dotplot(ego, showCategory = 20)

head(ego@result$Description, n=20)









# Perform Gene Ontology enrichment (GO BP for biological processes)
ego <- enrichGO(gene = gene_symbols,
                OrgDb = org.Hs.eg.db,
                keyType = "SYMBOL",
                ont = "BP",  # BP = Biological Process, adjust as needed
                pAdjustMethod = "BH",
                pvalueCutoff = 0.05)

# View the top enriched GO terms
head(ego)

barplot(ego, showCategory = 20) 

dotplot(ego, showCategory=15)

plotGOgraph(ego)






library(karyoploteR)

kp <- plotKaryotype(genome = "hg38",chromosomes="chr19")
kpPlotRegions(kp, data=c("chr19:38899326-38900000", "chr19:39369197-39390502"), col="#AACCFF", r0=0, r1=0.25)

kp <- plotKaryotype(genome = "mm10",chromosomes="chr19")
kpPlotRegions(kp, data=c("chr19:38899326-38900000", "chr19:39369197-39390502"), col=c("#AACCFF","green"), r0=0, r1=0.25)



# Load the libraries
library(GenomicRanges)
library(rtracklayer)
library(biomaRt)

# Load cCRE and syntenic blocks data from BED files
ccre_file <- "/Users/ryanhagan/NoCoSMiCC/Synteny/sc-lifted-syntenic-compI-G2-cCREs.bed"
syntenic_blocks_file <- "/Users/ryanhagan/NoCoSMiCC/Synteny/common_syntenic_blocks.bed"

# Import BED files as GRanges objects
ccre_bed = read.table(ccre_file)
ccre_gr <- makeGRangesFromDataFrame(unique(ccre_bed), keep.extra.columns=TRUE,
                                    start.field="V2",
                                    end.field = "V3",
                                    seqnames.field="V1")
ccre_gr

syntenic_blocks_bed = read.table(syntenic_blocks_file)
syntenic_blocks_gr <- makeGRangesFromDataFrame(syntenic_blocks_bed, keep.extra.columns=TRUE,
                                               start.field="V2",
                                               end.field = "V3",
                                               seqnames.field="V1")
syntenic_blocks_gr

# Load the gene annotation GTF file
gtf <- import("/Users/ryanhagan/NoCoSMiCC/TAD_Analysis/Files/gencode.v46.annotation.gtf.gz")

# Convert the GTF to GRanges format and keep only protein-coding genes
gtf_genes <- gtf[gtf$type == "gene" & gtf$gene_type == "protein_coding"]

genes_gr <- as(gtf_genes, "GRanges")
genes_gr

nearest_genes <- nearest(ccre_gr, genes_gr)
nearest_genes_info <- genes_gr[nearest_genes]
length(unique(nearest_genes_info$gene_name))

nearest_genes_info_df <- as.data.frame(nearest_genes_info)

# Combine nearest gene information with cCRE data
ccre_nearest_genes <- cbind(as.data.frame(ccre_gr), nearest_genes_info_df)

# Preview the result
head(ccre_nearest_genes)

# Syntenic blocks already loaded from BED
syntenic_gr <- GRanges(seqnames = seqnames(syntenic_blocks_gr),
                       ranges = IRanges(start = start(syntenic_blocks_gr), end = end(syntenic_blocks_gr)),
                       id = mcols(syntenic_blocks_gr)$name)
syntenic_gr <- unique(syntenic_gr)
# Nearest genes
genes_gr <- GRanges(seqnames = nearest_genes_info_df$seqnames,
                    ranges = IRanges(start = nearest_genes_info_df$start, end = nearest_genes_info_df$end),
                    gene = nearest_genes_info_df$gene)
genes_gr <- unique(genes_gr)
# Load Gviz library
library(Gviz)

genome_axis <- GenomeAxisTrack()

# Create Track for cCREs
ccre_track <- AnnotationTrack(ccre_gr, name = "cCREs", col = "#C2B280", fill = "#C2B280",
                              group = ccre_gr$V4, groupAnnotation = "group", just.group = "above",
                              showId = TRUE)

# Create Track for Nearest Genes (with labels)
gene_track <- AnnotationTrack(genes_gr, name = "Nearest Genes", col = "#00308F", fill = "#00308F",
                              group = genes_gr$gene_name, groupAnnotation = "group",
                              just.group = "above")

# Create Track for Syntenic Blocks
synteny_track <- AnnotationTrack(syntenic_gr, name = "Syntenic Blocks", col = "#88c4b4", fill = "#88c4b4")

# Define the plot range
plot_range <- GRanges(seqnames = "chr5", ranges = IRanges(start = 705246, end = 1295246))

# Plot the tracks with labeled genes
plotTracks(list(genome_axis, ccre_track, gene_track, synteny_track), 
           from = start(plot_range), to = end(plot_range),
           chromosome = as.character(seqnames(plot_range)))

plot_range <- GRanges(seqnames = "chr2", ranges = IRanges(start = 5553147, end = 13595068))

# Plot the tracks with labeled genes
plotTracks(list(genome_axis, ccre_track, gene_track, synteny_track), 
           from = start(plot_range), to = end(plot_range),
           chromosome = as.character(seqnames(plot_range)),
           rotation.item = 90)


### Updated version

# Load the libraries
library(GenomicRanges)
library(rtracklayer)
library(biomaRt)

gtf_genes <- gtf[gtf$type == "gene" & gtf$gene_type == "protein_coding"]

genes_gr <- as(gtf_genes, "GRanges")

# Find overlaps between cCREs and genes
overlaps <- findOverlaps(ccre_gr, genes_gr, type = "any", select = "all")

# Extract the indices of cCREs that overlap with genes
ccre_overlap_indices <- queryHits(overlaps)

# Get the indices of cCREs that do NOT overlap with genes
non_overlap_indices <- setdiff(seq_along(ccre_gr), unique(ccre_overlap_indices))

# Subset the cCREs based on the non-overlapping indices
nc_ccre_gr <- ccre_gr[non_overlap_indices]
nc_ccre_gr

# Load syntenic blocks data
syntenic_blocks_file <- "/Users/ryanhagan/NoCoSMiCC/Synteny/common_syntenic_blocks.bed"

syntenic_blocks_bed = read.table(syntenic_blocks_file)
syntenic_blocks_gr <- makeGRangesFromDataFrame(syntenic_blocks_bed, keep.extra.columns=TRUE,
                                               start.field="V2",
                                               end.field = "V3",
                                               seqnames.field="V1")
syntenic_blocks_gr


nearest_genes <- nearest(nc_ccre_gr, genes_gr)
nearest_genes_info <- genes_gr[nearest_genes]
length(unique(nearest_genes_info$gene_id))

nearest_genes_info_df <- as.data.frame(nearest_genes_info)

# Combine nearest gene information with cCRE data
ccre_nearest_genes <- cbind(as.data.frame(nc_ccre_gr), nearest_genes_info_df)

# Preview the result
head(ccre_nearest_genes)

# Syntenic blocks Granges
syntenic_gr <- GRanges(seqnames = seqnames(syntenic_blocks_gr),
                       ranges = IRanges(start = start(syntenic_blocks_gr), end = end(syntenic_blocks_gr)),
                       id = mcols(syntenic_blocks_gr)$name)
syntenic_gr <- unique(syntenic_gr)
# Genes Granges
#genes_gr <- GRanges(seqnames = gtf_genes@seqnames,
#                   ranges = IRanges(start = gtf_genes@start, end = gtf_genes@end),
#                  gene = gtf_genes$genename)
#genes_gr <- unique(genes_gr)
# Load Gviz library
library(Gviz)

genome_axis <- GenomeAxisTrack()

# subset the ccres by type
CA_only <- nc_ccre_gr[nc_ccre_gr$V5 == "CA-only" ]
PLS <- nc_ccre_gr[nc_ccre_gr$V5 == "PLS" | nc_ccre_gr$V5 == "PLS,CTCF-bound" ]
dELS <- nc_ccre_gr[nc_ccre_gr$V5 == "dELS" | nc_ccre_gr$V5 == "dELS,CTCF-bound" ]
pELS <- nc_ccre_gr[nc_ccre_gr$V5 == "pELS" | nc_ccre_gr$V5 == "pELS,CTCF-bound" ]
DNase_H3K4me3 <- nc_ccre_gr[nc_ccre_gr$V5 == "DNase-H3K4me3" | nc_ccre_gr$V5 == "DNase-H3K4me3,CTCF-bound" ]
CTCF <- nc_ccre_gr[nc_ccre_gr$V5 == "CTCF-only,CTCF-bound" ]

# Create Track for cCREs
ca_only_track <- AnnotationTrack(CA_only, name = "CA only", col = "#88c4b4", fill = "#88c4b4",
                                 group = CA_only$V4, groupAnnotation = "group", just.group = "above",
                                 showId = FALSE)
PLS_track <- AnnotationTrack(PLS, name = "PLS", col = "#88c4b4", fill = "#88c4b4",
                             group = PLS$V4, groupAnnotation = "group", just.group = "above",
                             showId = FALSE, shape = "ellipse")
pELS_track <- AnnotationTrack(pELS, name = "pELS", col = "#88c4b4", fill = "#88c4b4",
                              group = pELS$V4, groupAnnotation = "group", just.group = "above",
                              showId = FALSE)
dELS_track <- AnnotationTrack(dELS, name = "dELS", col = "#88c4b4", fill = "#88c4b4",
                              group = dELS$V4, groupAnnotation = "group", just.group = "above",
                              showId = FALSE)
dnase_h3k4me3_track <- AnnotationTrack(DNase_H3K4me3, name = "Dnase_H3K4me3", col = "#88c4b4", fill = "#C2B280",
                                       group = DNase_H3K4me3$V4, groupAnnotation = "group", just.group = "above",
                                       showId = FALSE)
CTCF_only_track <- AnnotationTrack(CTCF, name = "CTCF", col = "#88c4b4", fill = "#88c4b4",
                                   group = CTCF$V4, groupAnnotation = "group", just.group = "above",
                                   showId = FALSE)

# Create Track for Nearest Genes (with labels)
gene_track <- AnnotationTrack(gtf_genes, name = "Genes", col = "#87CEEB", fill = "#87CEEB",
                              group = gtf_genes$gene_name, groupAnnotation = "group",
                              just.group = "above", rotation.item = 90, track.height = 0.5)

# Create Track for Syntenic Blocks
synteny_track <- AnnotationTrack(syntenic_gr, name = "Syntenic Blocks", col = "lightgrey", fill = "lightgrey", track.height = 0.5)

# Define the plot range
plot_range <- GRanges(seqnames = "chr1", ranges = IRanges(start = 8878686, end = 8981000))

# Plot the tracks with labeled genes
# Plot the tracks with labeled genes and genome coordinates at the bottom
plotTracks(list(ca_only_track, PLS_track, pELS_track, dELS_track, dnase_h3k4me3_track, CTCF_only_track,
                gene_track, synteny_track, genome_axis),  # Place genome_axis at the end
           from = start(plot_range), to = end(plot_range),
           chromosome = as.character(seqnames(plot_range)),
           rotation.item = 90)

plot_range <- GRanges(seqnames = "chr2", ranges = IRanges(start = 8070000, end = 8990000))

# Plot the tracks with labeled genes
plotTracks(list(genome_axis, PLS_track, dELS_track, gene_track, synteny_track), 
           from = start(plot_range), to = end(plot_range),
           chromosome = as.character(seqnames(plot_range)),
           rotation.item = 90)


# Repeat for the nonG1 cCREs

other_ccre_file <- "/Users/ryanhagan/NoCoSMiCC/Zoonomia/liftover_compA_synt_nonG1_cCREs_sorted.bed"

# Import BED files as GRanges objects
other_ccre_bed = read.table(other_ccre_file)
other_ccre_gr <- makeGRangesFromDataFrame(unique(other_ccre_bed), keep.extra.columns=TRUE,
                                    start.field="V2",
                                    end.field = "V3",
                                    seqnames.field="V1")
other_ccre_gr

# Find overlaps between cCREs and genes
overlaps <- findOverlaps(other_ccre_gr, genes_gr, type = "any", select = "all")

# Extract the indices of cCREs that overlap with genes
ccre_overlap_indices <- queryHits(overlaps)

# Get the indices of cCREs that do NOT overlap with genes
non_overlap_indices <- setdiff(seq_along(other_ccre_gr), unique(ccre_overlap_indices))

# Subset the cCREs based on the non-overlapping indices
nc_other_ccre_gr <- other_ccre_gr[non_overlap_indices]
nc_other_ccre_gr

genome_axis <- GenomeAxisTrack()

# subset the ccres by type
other_CA_only <- nc_other_ccre_gr[nc_other_ccre_gr$V5 == "CA-only" ]
other_PLS <- nc_other_ccre_gr[nc_other_ccre_gr$V5 == "PLS" | nc_other_ccre_gr$V5 == "PLS,CTCF-bound" ]
other_dELS <- nc_other_ccre_gr[nc_other_ccre_gr$V5 == "dELS" | nc_other_ccre_gr$V5 == "dELS,CTCF-bound" ]
other_pELS <- nc_other_ccre_gr[nc_other_ccre_gr$V5 == "pELS" | nc_other_ccre_gr$V5 == "pELS,CTCF-bound" ]
other_DNase_H3K4me3 <- nc_other_ccre_gr[nc_other_ccre_gr$V5 == "DNase-H3K4me3" | nc_other_ccre_gr$V5 == "DNase-H3K4me3,CTCF-bound" ]
other_CTCF <- nc_other_ccre_gr[nc_other_ccre_gr$V5 == "CTCF-only,CTCF-bound" ]

# Create Track for the nonG1 cCREs
other_ca_only_track <- AnnotationTrack(other_CA_only, name = "CA only", col = "#e0115f", fill = "#e0115f",
                                 group = other_CA_only$V4, groupAnnotation = "group", just.group = "above",
                                 showId = FALSE)
other_PLS_track <- AnnotationTrack(other_PLS, name = "PLS", col = "#e0115f", fill = "#e0115f",
                             group = other_PLS$V4, groupAnnotation = "group", just.group = "above",
                             showId = FALSE)
other_pELS_track <- AnnotationTrack(other_pELS, name = "pELS", col = "#e0115f", fill = "#e0115f",
                              group = other_pELS$V4, groupAnnotation = "group", just.group = "above",
                              showId = FALSE)
other_dELS_track <- AnnotationTrack(other_dELS, name = "dELS", col = "#e0115f", fill = "#e0115f",
                              group = other_dELS$V4, groupAnnotation = "group", just.group = "above",
                              showId = FALSE)
other_dnase_h3k4me3_track <- AnnotationTrack(other_DNase_H3K4me3, name = "Dnase_H3K4me3", col = "#e0115f", fill = "#e0115f",
                                       group = other_DNase_H3K4me3$V4, groupAnnotation = "group", just.group = "above",
                                       showId = FALSE)
other_CTCF_only_track <- AnnotationTrack(other_CTCF, name = "CTCF", col = "#e0115f", fill = "#e0115f",
                                   group = other_CTCF$V4, groupAnnotation = "group", just.group = "above",
                                   showId = FALSE)

# Overlay CA-only tracks
ca_combined_track <- OverlayTrack(trackList = list(ca_only_track, other_ca_only_track), name = "CA only")

# Overlay PLS tracks
PLS_combined_track <- OverlayTrack(trackList = list(PLS_track, other_PLS_track), name = "PLS")

# Overlay pELS tracks
pELS_combined_track <- OverlayTrack(trackList = list(pELS_track, other_pELS_track), name = "pELS")

# Overlay dELS tracks
dELS_combined_track <- OverlayTrack(trackList = list(dELS_track, other_dELS_track), name = "dELS")

# Overlay DNase_H3K4me3 tracks
dnase_combined_track <- OverlayTrack(trackList = list(dnase_h3k4me3_track, other_dnase_h3k4me3_track), name = "Dnase H3K4me3")

# Overlay CTCF tracks
CTCF_combined_track <- OverlayTrack(trackList = list(CTCF_only_track, other_CTCF_only_track), name = "CTCF")

# Define the plot range
plot_range <- GRanges(seqnames = "chr8", ranges = IRanges(start = 24722952, end = 25622952))
 

plotTracks(list(ca_combined_track,
                PLS_combined_track,
                pELS_combined_track,
                dELS_combined_track,
                dnase_combined_track,
                CTCF_combined_track,
                gene_track, synteny_track, genome_axis),  
           from = start(plot_range), to = end(plot_range),
           chromosome = as.character(seqnames(plot_range)))


# plot all the human cCREs together
human_g1_ccres_track <- AnnotationTrack(nc_ccre_gr, name = "cCREs", col = "#C2B280", fill = "#C2B280",
                                       showId = FALSE)
human_nong1_ccres_track <- AnnotationTrack(nc_other_ccre_gr, name = "cCREs", col = "#C2B280", fill = "#C2B280",
                                        showId = FALSE)

all_human_cCREs_track <- OverlayTrack(trackList = list(human_g1_ccres_track, human_nong1_ccres_track), name = "cCREs")

plot_range <- GRanges(seqnames = "chr5", ranges = IRanges(start = 6022952, end = 7022952))
plotTracks(list(all_human_cCREs_track,
                gene_track, synteny_track, genome_axis),
           from = start(plot_range), to = end(plot_range),
           chromosome = as.character(seqnames(plot_range)))

# Wnt genes - Wnt11
plot_range <- GRanges(seqnames = "chr11", ranges = IRanges(start = 76116325, end = 76310761))
plotTracks(list(all_human_cCREs_track,
                gene_track, synteny_track, genome_axis),
           from = start(plot_range), to = end(plot_range),
           chromosome = as.character(seqnames(plot_range)))

# Wnt genes - Trabd2a
plot_range <- GRanges(seqnames = "chr2", ranges = IRanges(start = 84621650, end = 84987008))
plotTracks(list(all_human_cCREs_track,
                gene_track, synteny_track, genome_axis),
           from = start(plot_range), to = end(plot_range),
           chromosome = as.character(seqnames(plot_range)))

# Wnt genes - Wnt4
plot_range <- GRanges(seqnames = "chr1", ranges = IRanges(start = 21097313, end = 23149969))


## Now create the equivalent plots in mouse

# Load the mouse cCREs
mouse_cCREs_bed = read.table("/Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/LiftOver/cCRE_mouse_lifted")
mouse_cCREs_gr <- makeGRangesFromDataFrame(mouse_cCREs_bed, keep.extra.columns=TRUE,
                                    start.field="V2",
                                    end.field = "V3",
                                    seqnames.field="V1")
mouse_cCREs_gr


# Load the mouse GTF file
mouse_gtf <- import("/Users/ryanhagan/NoCoSMiCC/Files/gencode.vM10.annotation.gtf")

# Convert the GTF to GRanges format and keep only protein-coding genes
mouse_gtf_genes <- mouse_gtf[mouse_gtf$type == "gene" & mouse_gtf$gene_type == "protein_coding"]

mouse_gene_track <- AnnotationTrack(mouse_gtf_genes, name = "Genes", col = "#87CEEB", fill = "#87CEEB",
                              group = mouse_gtf_genes$gene_name, groupAnnotation = "group",
                              just.group = "above", track.height = 0.5)

genome_axis <- GenomeAxisTrack()

plot_range <- GRanges(seqnames = "chr14", ranges = IRanges(start = 67509377, end = 68599377))

plotTracks(list(mouse_gene_track, genome_axis),
           from = start(plot_range), to = end(plot_range),
           chromosome = as.character(seqnames(plot_range)))


# Write the non-coding ccre files 

nc_ccre_df <- data.frame(seqnames=seqnames(nc_ccre_gr),
                 starts=start(nc_ccre_gr)-1,
                 ends=end(nc_ccre_gr),
                 names=c(rep(".", length(nc_ccre_gr))),
                 scores=c(rep(".", length(nc_ccre_gr))),
                 strands=strand(nc_ccre_gr))
nc_ccre_df
write.table(nc_ccre_df, file="/Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/LiftOver/human_nc_ccre_G1.bed", quote=F, sep="\t", row.names=F, col.names=F)

nc_other_ccre_df <- data.frame(seqnames=seqnames(nc_other_ccre_gr),
                         starts=start(nc_other_ccre_gr)-1,
                         ends=end(nc_other_ccre_gr),
                         names=c(rep(".", length(nc_other_ccre_gr))),
                         scores=c(rep(".", length(nc_other_ccre_gr))),
                         strands=strand(nc_other_ccre_gr))
nc_other_ccre_df
write.table(nc_other_ccre_df, file="/Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/LiftOver/human_nc_ccre_nonG1.bed", quote=F, sep="\t", row.names=F, col.names=F)

# Read in the cCRE mouse liftover co-ordinates 

nc_mouse_ccre_bed = read.table("/Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/LiftOver/nc_cCRE_mouse_G1_lifted.bed")
nc_mouse_ccre_gr <- makeGRangesFromDataFrame(nc_mouse_ccre_bed, keep.extra.columns=TRUE,
                                             start.field="V2",
                                             end.field = "V3",
                                             seqnames.field="V1")
nc_mouse_ccre_gr

# Create Track for cCREs
mouse_ccre_track <- AnnotationTrack(nc_mouse_ccre_gr, name = "Mouse cCRE", col = "#88c4b4", fill = "#88c4b4",
                                 showId = FALSE)

genome_axis <- GenomeAxisTrack()

plot_range <- GRanges(seqnames = "chr14", ranges = IRanges(start = 67509377, end = 68599377))

plotTracks(list(mouse_ccre_track,
                mouse_gene_track, genome_axis),
           from = start(plot_range), to = end(plot_range),
           chromosome = as.character(seqnames(plot_range)))

nc_mouse_nonG1_ccre_bed = read.table("/Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/LiftOver/nc_cCRE_mouse_nonG1_lifted.bed")
nc_mouse_nonG1_ccre_gr <- makeGRangesFromDataFrame(nc_mouse_nonG1_ccre_bed, keep.extra.columns=TRUE,
                                             start.field="V2",
                                             end.field = "V3",
                                             seqnames.field="V1")
nc_mouse_nonG1_ccre_gr

mouse_ccre_nonG1_track <- AnnotationTrack(nc_mouse_nonG1_ccre_gr, name = "Mouse cCRE", col = "#e0115f", fill = "#e0115f",
                                    showId = FALSE)

mouse_combined_track <- OverlayTrack(trackList = list(mouse_ccre_track, mouse_ccre_nonG1_track), name = "Mouse cCRE")


plot_range <- GRanges(seqnames = "chr14", ranges = IRanges(start = 67509377, end = 68599377))

plotTracks(list(mouse_combined_track,
                mouse_gene_track, genome_axis),
           from = start(plot_range), to = end(plot_range),
           chromosome = as.character(seqnames(plot_range)))




## Plot some other regions - an entire syntenic block

# Transcript track
gtf_transcripts <- gtf[gtf$type == "gene"]
transcripts_gr <- as(gtf_transcripts, "GRanges")
transcript_track <- AnnotationTrack(transcripts_gr, name = "Transcripts", col = "#87CEEB", fill = "#87CEEB",
                                    group = transcripts_gr$gene_name, groupAnnotation = "group",
                                    just.group = "above")

# Human
plot_range <- GRanges(seqnames = "chr3", ranges = IRanges(start = 40804741, end = 41110096))
plotTracks(list(ca_combined_track,
                PLS_combined_track,
                pELS_combined_track,
                dELS_combined_track,
                dnase_combined_track,
                CTCF_combined_track,
                gene_track, transcript_track, synteny_track, genome_axis),
           from = start(plot_range), to = end(plot_range),
           chromosome = as.character(seqnames(plot_range)))





# Mouse 
plot_range <- GRanges(seqnames = "chr3", ranges = IRanges(start = 122729124, end = 142395636))

plotTracks(list(mouse_combined_track,
                mouse_gene_track, genome_axis),
           from = start(plot_range), to = end(plot_range),
           chromosome = as.character(seqnames(plot_range)))



### GO enrichment analysis
nearest_genes <- nearest(ccre_gr, genes_gr)
nearest_genes_info <- genes_gr[nearest_genes]
length(unique(nearest_genes_info$gene_id))

library(org.Hs.eg.db)

cleaned_gene_ids <- gsub("\\..*", "", nearest_genes_info$gene_id)
gene_symbols <- mapIds(org.Hs.eg.db, 
                       keys =cleaned_gene_ids, 
                       column = "SYMBOL", 
                       keytype = "ENSEMBL", 
                       multiVals = "first")

# Combine the results with cCRE data
ccre_nearest_genes <- cbind(as.data.frame(ccre_gr), gene_symbols)
#write.csv(ccre_nearest_genes, "/Users/ryanhagan/NoCoSMiCC/HiC_Analysis/constrained_ccres_nearest_genes.csv")
head(ccre_nearest_genes)
unique(ccre_nearest_genes$gene_symbols)
#write.csv(ccre_nearest_genes$gene_symbols, "/Users/ryanhagan/NoCoSMiCC/HiC_Analysis/constrained_ccres_nearest_genes.csv")

library(clusterProfiler)

# Perform GO enrichment analysis
ego <- enrichGO(gene = cleaned_gene_ids,
                OrgDb = org.Hs.eg.db,
                keyType = "ENSEMBL",
                ont = "BP",         # Biological Process
                pAdjustMethod = "BH",
                qvalueCutoff = 0.05)

# Visualize GO enrichment results
#barplot(ego, showCategory = 15)
dotplot(ego, showCategory = 15)

head(ego@result$Description, n=20)
ego_df <- as.data.frame(ego@result)
head(ego_df, n=15)
wnt_genes <- ego_df[15,]$geneID
wnt_genes

# repeat for nonG1

nearest_genes <- nearest(nc_other_ccre_gr, genes_gr)
nearest_genes_info <- genes_gr[nearest_genes]
length(unique(nearest_genes_info$gene_id))

cleaned_gene_ids <- gsub("\\..*", "", nearest_genes_info$gene_id)
gene_symbols <- mapIds(org.Hs.eg.db, 
                       keys =cleaned_gene_ids, 
                       column = "SYMBOL", 
                       keytype = "ENSEMBL", 
                       multiVals = "first")

ccre_nearest_genes <- cbind(as.data.frame(nc_other_ccre_gr), gene_symbols)
head(ccre_nearest_genes)
unique(ccre_nearest_genes$gene_symbols)

# Perform GO enrichment analysis
ego2 <- enrichGO(gene = cleaned_gene_ids,
                OrgDb = org.Hs.eg.db,
                keyType = "ENSEMBL",
                ont = "BP", 
                pAdjustMethod = "BH",
                qvalueCutoff = 0.05)

# Visualize GO enrichment results
barplot(ego2, showCategory = 15)
dotplot(ego2, showCategory = 20)

head(ego2@result$Description, n=20)

# combine the G1 and nonG1 genes
combined_genes <- union(ego@gene, ego2@gene)
length(combined_genes)

ego3 <- enrichGO(gene = combined_genes,
                 OrgDb = org.Hs.eg.db,
                 keyType = "ENSEMBL",
                 ont = "BP", 
                 pAdjustMethod = "BH",
                 qvalueCutoff = 0.05)

dotplot(ego3, showCategory = 20)
ego3_df <- as.data.frame(ego3@result)
head(ego3_df, n = 10)









#################################################################




##### create a hg38 TSS file for classifying cCREs

transcripts <- gtf[gtf$type == "transcript"]

# Extract the TSS coordinates
tss_gr <- transcripts
tss_gr

strand_plus <- as.vector(strand(tss_gr)) == "+"
strand_minus <- as.vector(strand(tss_gr)) == "-"

start(tss_gr)[strand_plus] <- start(tss_gr)[strand_plus]
start(tss_gr)[strand_minus] <- end(tss_gr)[strand_minus]
end(tss_gr) <- start(tss_gr) 

tss_gr <- tss_gr[, c("transcript_id", "gene_id")]
tss_df <- as.data.frame(tss_gr)

tss_df$chr <- as.character(seqnames(tss_gr))
tss_df$start <- start(tss_gr)
tss_df$end <- end(tss_gr)
tss_df$strand <- as.vector(strand(tss_gr))

head(tss_df)

# Rearrange columns to match the desired order
tss_df <- tss_df[, c("chr", "start", "end", "transcript_id", "strand", "gene_id")]

# Optionally add a column for the score (a dot to mimic BED format)
tss_df$score <- "."

# Rearrange again to place the score in the 5th position
tss_df <- tss_df[, c("chr", "start", "end", "transcript_id", "score", "strand", "gene_id")]

# View the first few rows to verify the format
head(tss_df)

# create the 4k bp TSS bed file
extended_tss_gr <- resize(tss_gr, width = width(tss_gr) + 4000, fix = "center")
extended_tss_gr

extended_tss_df <- as.data.frame(extended_tss_gr)

extended_tss_df$chr <- as.character(seqnames(extended_tss_gr))
extended_tss_df$start <- start(extended_tss_gr)
extended_tss_df$end <- end(extended_tss_gr)
extended_tss_df$strand <- as.vector(strand(extended_tss_gr))

head(extended_tss_df)

# Rearrange columns to match the desired order
extended_tss_df <- extended_tss_df[, c("chr", "start", "end", "transcript_id", "strand", "gene_id")]

# Optionally add a column for the score (a dot to mimic BED format)
extended_tss_df$score <- "."

# Rearrange again to place the score in the 5th position
extended_tss_df <- extended_tss_df[, c("chr", "start", "end", "transcript_id", "score", "strand", "gene_id")]

head(extended_tss_df)

# save the TSS file
write.table(tss_df, file = "/Users/ryanhagan/NoCoSMiCC/Files/hg38_tss.bed", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

# save the TSS +4k file
write.table(extended_tss_df, file = "/Users/ryanhagan/NoCoSMiCC/Files/hg38_tss_4k.bed", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

cactus \
--jobStore /path/to/jobStore \
/path/to/seqFile.txt \
/path/to/output/cactus.hal





### Colour plots

# Define a common color for ca, PLS, pELS, dELS, and dnase_CTCF tracks
common_background_color <- "#AA4A44"

# Apply the common background color and set text color to black for specific tracks
displayPars(ca_combined_track) <- list(background.title = common_background_color, col.title = "black")
displayPars(PLS_combined_track) <- list(background.title = common_background_color, col.title = "black")
displayPars(pELS_combined_track) <- list(background.title = common_background_color, col.title = "black")
displayPars(dELS_combined_track) <- list(background.title = common_background_color, col.title = "black")
displayPars(dnase_combined_track) <- list(background.title = common_background_color, col.title = "black")
displayPars(CTCF_combined_track) <- list(background.title = common_background_color, col.title = "black")

# Apply distinct background colors for the other tracks
displayPars(gene_track) <- list(background.title = "lightsteelblue", col.title = "black")
displayPars(synteny_track) <- list(background.title = "lightsteelblue", col.title = "black")
displayPars(genome_axis) <- list(background.title = "lightsteelblue", col.title = "black")

# Define the plot range
plot_range <- GRanges(seqnames = "chr8", ranges = IRanges(start = 24722952, end = 25622952))

# Plot the tracks with background title colors and black text
plotTracks(list(ca_combined_track,
                PLS_combined_track,
                pELS_combined_track,
                dELS_combined_track,
                dnase_combined_track,
                CTCF_combined_track,
                gene_track, synteny_track, genome_axis),
           from = start(plot_range), 
           to = end(plot_range),
           chromosome = as.character(seqnames(plot_range)))



###########################################################################
###########################################################################


##### GO analysis of filtered cCREs #######

# Load the libraries
library(GenomicRanges)
library(rtracklayer)
library(biomaRt)
library(org.Hs.eg.db)
library(clusterProfiler)

# Load cCRE file
ccre_file <- "/Users/ryanhagan/NoCoSMiCC/synteny/bulk-lifted-G1-cCREs.bed"
ccre_bed = read.table(ccre_file)
ccre_gr <- makeGRangesFromDataFrame(ccre_bed, keep.extra.columns=TRUE,
                                    start.field="V2",
                                    end.field = "V3",
                                    seqnames.field="V1")
ccre_gr

# Load the gene annotation GTF file
gtf <- import("/Users/ryanhagan/NoCoSMiCC/TAD_Analysis/Files/gencode.v46.annotation.gtf.gz")

gtf_genes <- gtf[gtf$type == "gene" & gtf$gene_type == "protein_coding"]

genes_gr <- as(gtf_genes, "GRanges")

# Load the GTF and extract exons for protein-coding genes
gtf_exons <- gtf[gtf$type == "exon" & gtf$gene_type == "protein_coding"]

# Convert the filtered exons into a GRanges object
exons_gr <- as(gtf_exons, "GRanges")

# Find overlaps between cCREs and genes
overlaps <- findOverlaps(ccre_gr, gtf_exons, type = "any", select = "all")

# Extract the indices of cCREs that overlap with genes/exons
ccre_overlap_indices <- queryHits(overlaps)

# Get the indices of cCREs that do NOT overlap with genes/exons
non_overlap_indices <- setdiff(seq_along(ccre_gr), unique(ccre_overlap_indices))

# Subset the cCREs based on the non-overlapping indices
nc_ccre_gr <- ccre_gr[non_overlap_indices]
nc_ccre_gr

coding_ccre_gr <- ccre_gr[ccre_overlap_indices]
coding_ccre_gr

### GO enrichment analysis
nearest_genes <- nearest(coding_ccre_gr, genes_gr)
nearest_genes_info <- genes_gr[nearest_genes]
length(unique(nearest_genes_info$gene_id))

cleaned_gene_ids <- gsub("\\..*", "", nearest_genes_info$gene_id)
gene_symbols <- mapIds(org.Hs.eg.db, 
                       keys =cleaned_gene_ids, 
                       column = "SYMBOL", 
                       keytype = "ENSEMBL", 
                       multiVals = "first")

# Combine the results with cCRE data
ccre_nearest_genes <- cbind(as.data.frame(coding_ccre_gr), gene_symbols)
#write.csv(ccre_nearest_genes, "/Users/ryanhagan/NoCoSMiCC/HiC_Analysis/constrained_ccres_nearest_genes.csv")
head(ccre_nearest_genes)
unique(ccre_nearest_genes$gene_symbols)
#write.csv(ccre_nearest_genes$gene_symbols, "/Users/ryanhagan/NoCoSMiCC/HiC_Analysis/constrained_ccres_nearest_genes.csv")

# Perform GO enrichment analysis
ego <- enrichGO(gene = cleaned_gene_ids,
                OrgDb = org.Hs.eg.db,
                keyType = "ENSEMBL",
                ont = "BP",         # Biological Process
                pAdjustMethod = "BH",
                qvalueCutoff = 0.05)

# Visualize GO enrichment results
#barplot(ego, showCategory = 15)
dotplot(ego, showCategory = 15)


ego2 <- enrichGO(gene = unique(ccre_nearest_genes$gene_symbols),
                OrgDb = org.Hs.eg.db,
                keyType = "SYMBOL",
                ont = "BP",         # Biological Process
                pAdjustMethod = "BH",
                qvalueCutoff = 0.05)

# Visualize GO enrichment results
#barplot(ego, showCategory = 15)
dotplot(ego2, showCategory = 15)



awk 'BEGIN{OFS="\t"} {gsub(/.*ColoncCRE/, "", $4); printf "ColoncCRE%s\t%.12f\n", $4, $5}' PhyloP_scores.bed | sort -k1,1 > sorted_PhyloP_scores.bed


echo "cCRE_ID scATAC_score H3K4me3_score H3K27ac_score CTCF_score PhyloP_score"
join hg38-scATAC-maxZ.txt hg38-H3K4me3-maxZ.txt > merged_1.txt
join merged_1.txt hg38-H3K27ac-maxZ.txt > merged_2.txt
join merged_2.txt hg38-CTCF-maxZ.txt > merged_3.txt
join merged_3.txt sorted_PhyloP_scores.bed > final_merged.txt




# syntenic cCREs between human, mouse and dog

# Use the approaches above to plot this - gtf file to get gene labels on tracks
# The liftover files contain the coordinates for lifted cCREs in mouse and dog
# Plot human, mouse and dog seperately


human_ccre_bed = read.table("/Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/sc-only-rCARs/MaxZ/hg38-cCREs-sc-only-filtered.bed")
human_ccre_gr <- makeGRangesFromDataFrame(human_ccre_bed, keep.extra.columns=TRUE,
                                    start.field="V2",
                                    end.field = "V3",
                                    seqnames.field="V1")
human_ccre_gr

# Load the gene annotation GTF file
gtf <- import("/Users/ryanhagan/NoCoSMiCC/TAD_Analysis/Files/gencode.v46.annotation.gtf.gz")
human_gtf_genes <- gtf[gtf$type == "gene" & gtf$gene_type == "protein_coding"]
genes_gr <- as(gtf_genes, "GRanges")

# Mouse
mouse_ccre_bed = read.table("/Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/LiftOver/hg38-cCREs-sc-mouse-lifted.bed")
mouse_ccre_gr <- makeGRangesFromDataFrame(mouse_ccre_bed, keep.extra.columns=TRUE,
                                          start.field="V2",
                                          end.field = "V3",
                                          seqnames.field="V1")
mouse_ccre_gr

gtf <- import("/Users/ryanhagan/NoCoSMiCC/TAD_Analysis/Files/gencode.vM10.annotation.gtf.gz")
mouse_gtf_genes <- gtf[gtf$type == "gene" & gtf$gene_type == "protein_coding"]
mouse_genes_gr <- as(gtf_genes, "GRanges")
mouse_genes_gr

# dog
dog_ccre_bed = read.table("/Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/LiftOver/hg38-cCREs-sc-dog-lifted.bed")
dog_ccre_gr <- makeGRangesFromDataFrame(dog_ccre_bed, keep.extra.columns=TRUE,
                                          start.field="V2",
                                          end.field = "V3",
                                          seqnames.field="V1")
dog_ccre_gr

gtf <- import("/Users/ryanhagan/NoCoSMiCC/TAD_Analysis/Files/Canis_lupus_familiaris.ROS_Cfam_1.0.113.gtf.gz")
dog_gtf_genes <- gtf[gtf$type == "gene" & gtf$gene_biotype == "protein_coding"]
dog_genes_gr <- as(gtf_genes, "GRanges")
dog_genes_gr

# modify seqnames of dog genome
current_seqlevels <- seqlevels(dog_gtf_genes)
new_seqlevels <- ifelse(grepl("^chr", current_seqlevels),
                        current_seqlevels,
                        paste0("chr", current_seqlevels))

# Set the new seqlevels to the GRanges object
seqlevels(dog_gtf_genes) <- new_seqlevels


library(Gviz)

genome_axis <- GenomeAxisTrack()


# Create Track for Nearest Genes (with labels)
human_gene_track <- AnnotationTrack(human_gtf_genes, name = "Genes", col = "#87CEEB", fill = "#87CEEB",
                              group = human_gtf_genes$gene_name, groupAnnotation = "group",
                              just.group = "above", rotation.item = 90, track.height = 0.5)

mouse_gene_track <- AnnotationTrack(mouse_gtf_genes, name = "Genes", col = "#87CEEB", fill = "#87CEEB",
                                    group = mouse_gtf_genes$gene_name, groupAnnotation = "group",
                                    just.group = "above", rotation.item = 90, track.height = 0.5)

dog_gene_track <- AnnotationTrack(dog_gtf_genes, name = "Genes", col = "#87CEEB", fill = "#87CEEB",
                                    group = dog_gtf_genes$gene_name, groupAnnotation = "group",
                                    just.group = "above", rotation.item = 90, track.height = 0.5,
                                  options(ucscChromosomeNames=FALSE))


# subset the ccres by type
human_dELS <- human_ccre_gr[human_ccre_gr$V5 == "dELS" | human_ccre_gr$V5 == "dELS,CTCF-bound" ]
mouse_dELS <- mouse_ccre_gr[mouse_ccre_gr$V5 == "dELS" | mouse_ccre_gr$V5 == "dELS,CTCF-bound" ]
dog_dELS <- dog_ccre_gr[dog_ccre_gr$V5 == "dELS" | dog_ccre_gr$V5 == "dELS,CTCF-bound" ]



# Create Track for cCREs
human_pELS_track <- AnnotationTrack(human_dELS, name = "human_dELS", col = "#88c4b4", fill = "#88c4b4",
                              group = human_dELS$V4, groupAnnotation = "group", just.group = "above",
                              showId = FALSE)

mouse_pELS_track <- AnnotationTrack(mouse_dELS, name = "mouse_dELS", col = "#88c4b4", fill = "#88c4b4",
                                    group = mouse_dELS$V4, groupAnnotation = "group", just.group = "above",
                                    showId = FALSE)

dog_pELS_track <- AnnotationTrack(dog_dELS, name = "dog_dELS", col = "#88c4b4", fill = "#88c4b4",
                                    group = dog_dELS$V4, groupAnnotation = "group", just.group = "above",
                                    showId = FALSE)




overlap_ccres <- intersect(human_ccre_gr$V4, mouse_ccre_gr$V4)
overlap_ccres_all <- intersect(overlap_ccres, dog_ccre_gr$V4)
overlap_ccres_all

human_ccre_filtered <- human_ccre_gr[
  mcols(human_ccre_gr)$V4 %in% overlap_ccres_all & 
    mcols(human_ccre_gr)$V5 == "dELS"
]
head(human_ccre_filtered, n = 1176)

# Define the plot range for human
plot_range <- GRanges(seqnames = "chr2", ranges = IRanges(start = 47302297, end = 47684740))

# Plot
plotTracks(list(human_pELS_track, human_gene_track, genome_axis),  # Place genome_axis at the end
           from = start(plot_range), to = end(plot_range),
           chromosome = as.character(seqnames(plot_range)),
           rotation.item = 90)

# Define the plot range for mouse
plot_range <- GRanges(seqnames = "chr17", ranges = IRanges(start = 87543407, end = 87858555))

# Plot
plotTracks(list(mouse_pELS_track, mouse_gene_track, genome_axis),  # Place genome_axis at the end
           from = start(plot_range), to = end(plot_range),
           chromosome = as.character(seqnames(plot_range)),
           rotation.item = 90)


dog_prkcz_gene <- subset(dog_gtf_genes, gene_name == "EPCAM")
dog_prkcz_gene

# Define the plot range for dog
plot_range <- GRanges(seqnames = "chr10", ranges = IRanges(start = 50200727, end = 50549285))

# Plot
plotTracks(list(dog_pELS_track, dog_gene_track, genome_axis),  # Place genome_axis at the end
           from = start(plot_range), to = end(plot_range),
           chromosome = as.character(seqnames(plot_range)),
           rotation.item = 90)






######
# PDZK1IP1


# human
plot_range <- GRanges(seqnames = "chr1", ranges = IRanges(start = 47049265, end = 48756716))

# Plot
plotTracks(list(human_pELS_track, human_gene_track, genome_axis),  # Place genome_axis at the end
           from = start(plot_range), to = end(plot_range),
           chromosome = as.character(seqnames(plot_range)),
           rotation.item = 90)

# mouse
plot_range <- GRanges(seqnames = "chr4", ranges = IRanges(start = 112945905, end = 116951091))

# Plot
plotTracks(list(mouse_pELS_track, mouse_gene_track, genome_axis),  # Place genome_axis at the end
           from = start(plot_range), to = end(plot_range),
           chromosome = as.character(seqnames(plot_range)),
           rotation.item = 90)

# dog
plot_range <- GRanges(seqnames = "chr15", ranges = IRanges(start = 12929353, end = 13735332))

# Plot
plotTracks(list(dog_pELS_track, dog_gene_track, genome_axis),  # Place genome_axis at the end
           from = start(plot_range), to = end(plot_range),
           chromosome = as.character(seqnames(plot_range)),
           rotation.item = 90)

######
# SSBP2

# human
plot_range <- GRanges(seqnames = "chr5", ranges = IRanges(start = 81000000, end = 83000000))

# Plot
plotTracks(list(human_pELS_track, human_gene_track, genome_axis),  # Place genome_axis at the end
           from = start(plot_range), to = end(plot_range),
           chromosome = as.character(seqnames(plot_range)),
           rotation.item = 90)
# mouse
plot_range <- GRanges(seqnames = "chr13", ranges = IRanges(start = 90608402, end = 92851548))

# Plot
plotTracks(list(mouse_pELS_track, mouse_gene_track, genome_axis),  # Place genome_axis at the end
           from = start(plot_range), to = end(plot_range),
           chromosome = as.character(seqnames(plot_range)),
           rotation.item = 90)

# dog
plot_range <- GRanges(seqnames = "chr3", ranges = IRanges(start = 22530545, end = 28815874))

# Plot
plotTracks(list(dog_pELS_track, dog_gene_track, genome_axis),  # Place genome_axis at the end
           from = start(plot_range), to = end(plot_range),
           chromosome = as.character(seqnames(plot_range)),
           rotation.item = 90)



## Function to plot cCREs

library(GenomicRanges)
library(Gviz)
library(rtracklayer)

# Set a consistent color palette
COLOR_PALETTE <- list(
  genes = "#2E8B57",       # Sea Green for genes
  dELS = "#4169E1",        # Royal Blue for dELS tracks
  axis = "#333333",        # Dark gray for genome axis
  background = "#F5F5F5"   # Light gray background
)

# Enhanced function to create gene tracks
create_enhanced_gene_track <- function(gtf_genes, track_name) {
  AnnotationTrack(
    gtf_genes, 
    name = track_name, 
    col = COLOR_PALETTE$genes, 
    fill = COLOR_PALETTE$genes,
    group = gtf_genes$gene_name, 
    groupAnnotation = "group",
    just.group = "above", 
    rotation.item = 90, 
    track.height = 3,  # Increased track height
    shape = "box",     # More defined box shape
    fontsize = 10,     # Adjusted font size
    fontcolor = "black"
  )
}

# Enhanced function to create dELS tracks
create_enhanced_dELS_track <- function(ccre_gr, track_name) {
  AnnotationTrack(
    ccre_gr, 
    name = track_name, 
    col = COLOR_PALETTE$dELS, 
    fill = COLOR_PALETTE$dELS,
    group = ccre_gr$V4, 
    groupAnnotation = "group", 
    just.group = "above",
    showId = FALSE,
    track.height = 2,  # Adjusted track height
    shape = "box",
    alpha = 0.7        # Slight transparency
  )
}

# Prepare data loading function
prepare_species_data <- function(gtf_path, ccre_bed_path, 
                                 chromosome_conversion = FALSE) {
  # Load GTF
  gtf <- import(gtf_path)
  gtf_genes <- gtf[gtf$type == "gene" & gtf$gene_type == "protein_coding"]
  
  # Optional chromosome name conversion (e.g., for dog genome)
  if (chromosome_conversion) {
    current_seqlevels <- seqlevels(gtf_genes)
    new_seqlevels <- ifelse(grepl("^chr", current_seqlevels),
                            current_seqlevels,
                            paste0("chr", current_seqlevels))
    seqlevels(gtf_genes) <- new_seqlevels
  }
  
  # Load cCREs
  ccre_bed <- read.table(ccre_bed_path)
  ccre_gr <- makeGRangesFromDataFrame(ccre_bed, 
                                      keep.extra.columns = TRUE,
                                      start.field = "V2",
                                      end.field = "V3",
                                      seqnames.field = "V1")
  
  # Subset dELS tracks
  dELS <- ccre_gr[ccre_gr$V5 == "dELS" | ccre_gr$V5 == "dELS,CTCF-bound"]
  
  return(list(
    gtf_genes = gtf_genes,
    dELS = dELS
  ))
}

# Prepare species data with chromosome conversion for dog
human_data <- prepare_species_data(
  "/Users/ryanhagan/NoCoSMiCC/TAD_Analysis/Files/gencode.v46.annotation.gtf.gz",
  "/Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/sc-only-rCARs/MaxZ/hg38-cCREs-sc-only-filtered.bed"
)

mouse_data <- prepare_species_data(
  "/Users/ryanhagan/NoCoSMiCC/TAD_Analysis/Files/gencode.vM10.annotation.gtf.gz",
  "/Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/LiftOver/hg38-cCREs-sc-mouse-lifted.bed"
)

dog_data <- prepare_species_data(
  "/Users/ryanhagan/NoCoSMiCC/TAD_Analysis/Files/Canis_lupus_familiaris.ROS_Cfam_1.0.113.gtf.gz",
  "/Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/LiftOver/hg38-cCREs-sc-dog-lifted.bed",
  chromosome_conversion = TRUE
)

# Function to find a specific gene and get its genomic coordinates
find_gene_coordinates <- function(gtf, gene_name) {
  gene_subset <- gtf[gtf$gene_name == gene_name]
  if(length(gene_subset) == 0) {
    stop(paste("Gene", gene_name, "not found in the GTF file"))
  }
  return(gene_subset)
}

# Function to create a plot range around a gene
create_gene_plot_range <- function(gene_gr, flank_size = 2000) {
  GRanges(
    seqnames = seqnames(gene_gr),
    ranges = IRanges(
      start = start(gene_gr) - flank_size, 
      end = end(gene_gr) + flank_size
    )
  )
}

# Function to plot syntenic regions for a specific species
plot_species_region <- function(species_data, plot_range, species_name) {
  # Filter dELS tracks within the plot range
  dELS_track <- create_enhanced_dELS_track(
    species_data$dELS[overlapsAny(species_data$dELS, plot_range)], 
    paste(species_name, "dELS")
  )
  
  # Create gene track for genes within the plot range
  genes_in_range <- species_data$gtf_genes[overlapsAny(species_data$gtf_genes, plot_range)]
  gene_track <- create_enhanced_gene_track(genes_in_range, paste(species_name, "Genes"))
  
  # Genome axis
  genome_axis <- GenomeAxisTrack(
    fontcolor = COLOR_PALETTE$axis,
    col = COLOR_PALETTE$axis,
    fontsize = 10
  )
  
  # Plot
  plotTracks(
    list(dELS_track, gene_track, genome_axis),
    from = start(plot_range), 
    to = end(plot_range),
    chromosome = as.character(seqnames(plot_range))
  )
}

# Wrapper function to plot regions for multiple species
plot_syntenic_regions <- function(gene_name = "APC", 
                                  species_data_list, 
                                  flank_size = 5000) {
  # Find gene in human (reference species)
  human_gene <- find_gene_coordinates(species_data_list$human$gtf_genes, gene_name)
  
  # Create plot range
  human_plot_range <- create_gene_plot_range(human_gene, flank_size)
  
  # Plot for each species
  plot_species_region(species_data_list$human, human_plot_range, "Human")
  plot_species_region(species_data_list$mouse, human_plot_range, "Mouse")
  plot_species_region(species_data_list$dog, human_plot_range, "Dog")
}

# Example usage
plot_syntenic_regions("APC", 
                      list(
                        human = human_data,
                        mouse = mouse_data,
                        dog = dog_data
                      ),
                      flank_size = 5000  # Easy to modify flank size
)





#######

########

library(GenomicRanges)
library(Gviz)
library(rtracklayer)

# Function to standardize chromosome names
standardize_chromosome_names <- function(gr) {
  # Remove 'chr' prefix if present
  seqlevels(gr) <- gsub("^chr", "", seqlevels(gr))
  
  # Standard chromosome names for human, mouse, and dog
  standard_chroms <- c(as.character(1:38), "X", "Y", "M")
  
  # Filter to keep only standard chromosomes
  gr <- gr[seqnames(gr) %in% standard_chroms]
  
  # Add 'chr' prefix back
  seqlevels(gr) <- paste0("chr", seqlevels(gr))
  
  return(gr)
}

# Set a consistent color palette
COLOR_PALETTE <- list(
  genes = "#2E8B57",       # Sea Green for genes
  dELS = "#4169E1",        # Royal Blue for dELS tracks
  axis = "#333333",        # Dark gray for genome axis
  background = "#F5F5F5"   # Light gray background
)

# Enhanced function to create gene tracks
create_enhanced_gene_track <- function(gtf_genes, track_name) {
  AnnotationTrack(
    gtf_genes,
    name = track_name,
    col = COLOR_PALETTE$genes,
    fill = COLOR_PALETTE$genes,
    group = gtf_genes$gene_name,
    groupAnnotation = "group",
    just.group = "above",
    rotation.item = 90,
    track.height = 3,  # Increased track height
    shape = "box",     # More defined box shape
    fontsize = 10,     # Adjusted font size
    fontcolor = "black"
  )
}

# Enhanced function to create dELS tracks
create_enhanced_dELS_track <- function(ccre_gr, track_name) {
  AnnotationTrack(
    ccre_gr,
    name = track_name,
    col = COLOR_PALETTE$dELS,
    fill = COLOR_PALETTE$dELS,
    group = ccre_gr$V4,
    groupAnnotation = "group",
    just.group = "above",
    showId = FALSE,
    track.height = 2,  # Adjusted track height
    shape = "box",
    alpha = 0.7        # Slight transparency
  )
}

# Prepare data loading function
prepare_species_data <- function(gtf_path, ccre_bed_path) {
  # Load GTF
  gtf <- import(gtf_path)
  gtf_genes <- gtf[gtf$type == "gene" & gtf$gene_type == "protein_coding"]
  
  # Standardize chromosome names
  gtf_genes <- standardize_chromosome_names(gtf_genes)
  
  # Load cCREs
  ccre_bed <- read.table(ccre_bed_path)
  ccre_gr <- makeGRangesFromDataFrame(ccre_bed,
                                      keep.extra.columns = TRUE,
                                      start.field = "V2",
                                      end.field = "V3",
                                      seqnames.field = "V1")
  
  # Standardize chromosome names for cCREs
  ccre_gr <- standardize_chromosome_names(ccre_gr)
  
  # Subset dELS tracks
  dELS <- ccre_gr[ccre_gr$V5 == "dELS" | ccre_gr$V5 == "dELS,CTCF-bound"]
  
  return(list(
    gtf_genes = gtf_genes,
    dELS = dELS
  ))
}

# Prepare species data with specific paths
human_data <- prepare_species_data(
  "/Users/ryanhagan/NoCoSMiCC/TAD_Analysis/Files/gencode.v46.annotation.gtf.gz",
  "/Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/sc-only-rCARs/MaxZ/hg38-cCREs-sc-only-filtered.bed"
)

mouse_data <- prepare_species_data(
  "/Users/ryanhagan/NoCoSMiCC/TAD_Analysis/Files/gencode.vM10.annotation.gtf.gz",
  "/Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/LiftOver/hg38-cCREs-sc-mouse-lifted.bed"
)

dog_data <- prepare_species_data(
  "/Users/ryanhagan/NoCoSMiCC/TAD_Analysis/Files/Canis_lupus_familiaris.ROS_Cfam_1.0.113.gtf.gz",
  "/Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/LiftOver/hg38-cCREs-sc-dog-lifted.bed"
)

# Function to find a specific gene and get its genomic coordinates
find_gene_coordinates <- function(gtf, gene_name) {
  gene_subset <- gtf[gtf$gene_name == gene_name]
  if(length(gene_subset) == 0) {
    stop(paste("Gene", gene_name, "not found in the GTF file"))
  }
  return(gene_subset)
}

# Function to create a plot range around a gene
create_gene_plot_range <- function(gene_gr, flank_size = 2000) {
  GRanges(
    seqnames = seqnames(gene_gr),
    ranges = IRanges(
      start = start(gene_gr) - flank_size,
      end = end(gene_gr) + flank_size
    )
  )
}

# Function to plot syntenic regions for a specific species
plot_species_region <- function(species_data, plot_range, species_name) {
  # Filter dELS tracks within the plot range
  dELS_track <- create_enhanced_dELS_track(
    species_data$dELS[overlapsAny(species_data$dELS, plot_range)],
    paste(species_name, "dELS")
  )
  
  # Create gene track for genes within the plot range
  genes_in_range <- species_data$gtf_genes[overlapsAny(species_data$gtf_genes, plot_range)]
  gene_track <- create_enhanced_gene_track(genes_in_range, paste(species_name, "Genes"))
  
  # Genome axis
  genome_axis <- GenomeAxisTrack(
    fontcolor = COLOR_PALETTE$axis,
    col = COLOR_PALETTE$axis,
    fontsize = 10
  )
  
  # Plot
  plotTracks(
    list(dELS_track, gene_track, genome_axis),
    from = start(plot_range),
    to = end(plot_range),
    chromosome = as.character(seqnames(plot_range))
  )
}

# Wrapper function to plot regions for multiple species
plot_syntenic_regions <- function(gene_name = "APC",
                                  species_data_list,
                                  flank_size = 5000) {
  # Find gene in human (reference species)
  human_gene <- find_gene_coordinates(species_data_list$human$gtf_genes, gene_name)
  
  # Create plot range
  human_plot_range <- create_gene_plot_range(human_gene, flank_size)
  
  # Plot for each species
  plot_species_region(species_data_list$human, human_plot_range, "Human")
  plot_species_region(species_data_list$mouse, human_plot_range, "Mouse")
  plot_species_region(species_data_list$dog, human_plot_range, "Dog")
}

# Function to plot syntenic regions with predefined species data
plot_syntenic_regions_with_predefined_data <- function(
    gene_name = "APC", 
    flank_size = 5000
) {
  # Combine species data
  species_data_list <- list(
    human = human_data,
    mouse = mouse_data,
    dog = dog_data
  )
  
  # Plot syntenic regions
  plot_syntenic_regions(
    gene_name = gene_name, 
    species_data_list = species_data_list, 
    flank_size = flank_size
  )
}

# Example usage
# Plot APC gene with 5000 bp padding
plot_syntenic_regions_with_predefined_data(
  gene_name = "APC",  # Change this to your gene of interest
  flank_size = 50000   # Adjust padding as needed
)


plot_syntenic_regions_with_predefined_data("TP53", flank_size = 50000)







library(GenomicRanges)
library(Gviz)
library(rtracklayer)

# Modified color palette for cCRE types
COLOR_PALETTE <- list(
  cCRE = list(
    dELS = "#4169E1",
    pELS = "#6A5ACD",
    PLS = "#32CD32",
    CTCF = "#FF4500",
    H3K4me3 = "#FFD700",
    DNase = "#40E0D0",
    default = "#B0C4DE"
  ),
  axis = "#333333",
  highlight = "#FF0000"
)

# Modified prepare_species_data to handle full cCRE data
prepare_species_data <- function(gtf_path, ccre_bed_path) {
  # Load GTF (retained for potential gene context)
  gtf <- import(gtf_path)
  gtf_genes <- gtf[gtf$type == "gene" & gtf$gene_type == "protein_coding"]
  gtf_genes <- standardize_chromosome_names(gtf_genes)
  
  # Load full cCRE data
  ccre_bed <- read.table(ccre_bed_path, col.names = c("chr", "start", "end", "ccre_name", "type"))
  ccre_gr <- makeGRangesFromDataFrame(ccre_bed,
                                      keep.extra.columns = TRUE,
                                      start.field = "start",
                                      end.field = "end",
                                      seqnames.field = "chr")
  
  # Standardize chromosome names
  ccre_gr <- standardize_chromosome_names(ccre_gr)
  
  return(list(
    ccre = ccre_gr,
    genes = gtf_genes  # Retained for optional gene context
  ))
}

# Function to find specific cCRE coordinates
find_ccre_coordinates <- function(ccre_gr, ccre_name) {
  ccre_subset <- ccre_gr[ccre_gr$ccre_name == ccre_name]
  if(length(ccre_subset) == 0) {
    stop(paste("cCRE", ccre_name, "not found in the dataset"))
  }
  return(ccre_subset)
}

# Create plot range around cCRE
create_ccre_plot_range <- function(ccre_gr, flank_size = 5000) {
  GRanges(
    seqnames = seqnames(ccre_gr),
    ranges = IRanges(
      start = start(ccre_gr) - flank_size,
      end = end(ccre_gr) + flank_size
    )
  )
}

# Enhanced cCRE track creation with type-specific coloring
create_ccre_track <- function(ccre_gr, track_name, highlight_ccre = NULL) {
  # Assign colors based on cCRE type
  type_colors <- sapply(ccre_gr$type, function(t) {
    COLOR_PALETTE$cCRE[[unlist(strsplit(t, ","))[1]]] %||% COLOR_PALETTE$cCRE$default
  })
  
  AnnotationTrack(
    ccre_gr,
    name = track_name,
    col = type_colors,
    fill = type_colors,
    group = ccre_gr$ccre_name,
    groupAnnotation = "group",
    just.group = "above",
    shape = "box",
    rotation.item = 45,
    cex.group = 0.8,
    background.title = "transparent"
  )
}

# Main plotting function
plot_ccre_species <- function(ccre, flank_size = 5000, species = c("human", "mouse", "dog")) {
  species <- match.arg(species)
  species_data <- switch(species,
                         human = human_data,
                         mouse = mouse_data,
                         dog = dog_data
  )
  
  # Get cCRE coordinates
  target_ccre <- find_ccre_coordinates(species_data$ccre, ccre)
  plot_range <- create_ccre_plot_range(target_ccre, flank_size)
  
  # Create tracks
  genome_axis <- GenomeAxisTrack(
    fontcolor = COLOR_PALETTE$axis,
    col = COLOR_PALETTE$axis,
    fontsize = 10
  )
  
  ccre_track <- create_ccre_track(
    species_data$ccre[overlapsAny(species_data$ccre, plot_range)],
    paste(species, "cCREs")
  )
  
  # Highlight the target cCRE
  highlight_track <- HighlightTrack(
    trackList = list(ccre_track, genome_axis),
    start = start(target_ccre),
    end = end(target_ccre),
    chromosome = as.character(seqnames(target_ccre)),
    col = COLOR_PALETTE$highlight,
    fill = NA,
    inBackground = FALSE
  )
  
  # Plot
  plotTracks(
    highlight_track,
    from = start(plot_range),
    to = end(plot_range),
    chromosome = as.character(seqnames(plot_range)),
    main = paste("cCRE Landscape:", ccre, "in", species),
    sizes = c(2, 0.5)
  )
}

# Example usage:
plot_ccre_species("ColoncCRE2435", flank_size = 10000, species = "human")





##### ggbio visualisation ######

library(ggbio)
library(GenomicRanges)
library(ggplot2)
library(patchwork)

# Define the plotting region for human (chr2)
plot_range <- GRanges(seqnames = "chr2", 
                      ranges = IRanges(start = 47302297, end = 47684740))

# Subset gene annotations and cCREs to the plotting region
human_genes_plot <- subsetByOverlaps(human_gtf_genes, plot_range)
human_ccre_plot  <- subsetByOverlaps(human_ccre_filtered, plot_range)

# Convert the GRanges objects to data frames for plotting
genes_df <- as.data.frame(human_genes_plot)
ccre_df  <- as.data.frame(human_ccre_plot)

# Extract cCRE names from column V4
ccre_df$cCRE_name <- ccre_df$V4  

### Compute the distance to the nearest TSS and annotate nearest gene
human_TSS <- promoters(human_gtf_genes, upstream = 0, downstream = 1)
nearest_hits <- distanceToNearest(human_ccre_plot, human_TSS)

# Extract gene names directly from human_TSS
ccre_df$nearest_gene <- mcols(human_TSS)$gene_name[subjectHits(nearest_hits)]
ccre_df$nearest_tss_distance <- mcols(nearest_hits)$distance

# Create annotation label: "Nearest gene = GENE (distance = X bp)"
ccre_df$annotation_label <- paste0("Nearest gene = ", ccre_df$nearest_gene, 
                                   " (", ccre_df$nearest_tss_distance, " bp)")

### Adjust label positions for cCRE annotations
ccre_df$x_cCRE_label <- (ccre_df$start + ccre_df$end) / 2 - 20 * (ccre_df$end - ccre_df$start)  # Shift left
ccre_df$x_nearest_label <- (ccre_df$start + ccre_df$end) / 2 + 150 * (ccre_df$end - ccre_df$start)  # Shift right
ccre_df$y_nearest_label <- 0.05  

### Gene Plot with gene labels just above the arrows
gene_plot <- ggplot(genes_df, aes(xmin = start, xmax = end, y = gene_name)) +
  geom_gene_arrow(fill = "#800020") +
  # Compute the midpoint for each gene and nudge labels upward
  geom_text(aes(x = (start + end) / 2, label = gene_name), 
            position = position_nudge(y = 0.3),
            size = 3,
            color = "black") +
  theme_genes() +
  scale_x_continuous(limits = c(start(plot_range), end(plot_range))) +
  theme(
    axis.title.x = element_blank(),
    axis.text.x  = element_text(size = 10),
    axis.ticks.x = element_line(),  
    axis.title.y = element_blank(),
    axis.text.y  = element_blank(),
    axis.ticks.y = element_line(),
    plot.title = element_blank()
  )

### cCRE Track Plot with dynamic title using different colors
ccre_plot <- ggplot(ccre_df) + 
  geom_rect(aes(xmin = start, xmax = end, ymin = 0.45, ymax = 0.55), 
            fill = "#4CBB17", color = "#4CBB17") +
  scale_y_continuous(limits = c(0.4, 0.6)) +
  ggtitle(bquote(atop(
    bold(.(ccre_df$cCRE_name)) * 
      bold(" | Nearest gene = ") * 
      .(ccre_df$nearest_gene) * 
      bold(" | Distance = ") * 
      bold(.(ccre_df$nearest_tss_distance)) * " bp"
  ))) +
  scale_x_continuous(limits = c(start(plot_range), end(plot_range))) + 
  theme_classic() + 
  theme(
    axis.title.x = element_blank(),
    axis.text.x  = element_text(size = 10),
    axis.ticks.x = element_line(),
    axis.title.y = element_blank(),
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank(), 
    axis.line.x = element_line(),
    plot.title = element_text(hjust = 0.5, size = 12)
  )

### Combine the Plots with Patchwork
combined_plot <- ccre_plot / gene_plot +
  plot_layout(heights = c(1, 3)) +
  plot_annotation(
    title = "",
    theme = theme(
      plot.title = element_text(hjust = 0.5, size = 14),
      panel.border = element_rect(color = "black", fill = NA, size = 1),  # Add a box
      plot.margin = margin(10, 10, 10, 10)  # Add space around the plot
    )
  )

# Display the combined plot
print(combined_plot)

## add a border

# Add a border around the entire combined plot
combined_plot <- ccre_plot / gene_plot + 
  plot_layout(heights = c(0.5, 2)) + 
  plot_annotation(
    title = "",
    theme = theme(
      plot.title = element_text(hjust = 0.5, size = 14),
      plot.background = element_rect(color = "black", fill = NA, size = 1)  # Border added here
    )
  )

# Display the final plot with a box
print(combined_plot)




### FUNCTION
library(ggbio)
library(GenomicRanges)
library(ggplot2)
library(patchwork)

plot_cCRE_region <- function(ccre_id, flank = 50000, human_gtf_genes, human_ccre_filtered) {
  # Subset the cCRE based on the given ID
  ccre_selected <- human_ccre_filtered[human_ccre_filtered$V4 == ccre_id]
  
  if (length(ccre_selected) == 0) {
    stop("Error: cCRE ID not found in dataset.")
  }
  
  # Define the plot range based on flanking distance
  plot_range <- GRanges(
    seqnames = seqnames(ccre_selected),
    ranges = IRanges(
      start = start(ccre_selected) - flank,
      end = end(ccre_selected) + flank
    )
  )
  
  # Subset gene annotations and cCREs to the plotting region
  human_genes_plot <- subsetByOverlaps(human_gtf_genes, plot_range)
  human_ccre_plot  <- subsetByOverlaps(human_ccre_filtered, plot_range)
  
  # Convert GRanges to data frames
  genes_df <- as.data.frame(human_genes_plot)
  ccre_df  <- as.data.frame(human_ccre_plot)
  
  # Extract cCRE names
  ccre_df$cCRE_name <- ccre_df$V4  
  
  # Identify the selected cCRE for color mapping
  ccre_df$highlight <- ifelse(ccre_df$cCRE_name == ccre_id, "selected", "other")
  
  # Compute the distance to the nearest TSS and annotate nearest gene
  human_TSS <- promoters(human_gtf_genes, upstream = 0, downstream = 1)
  nearest_hits <- distanceToNearest(human_ccre_plot, human_TSS)
  
  # Extract gene names from human_TSS
  ccre_df$nearest_gene <- mcols(human_TSS)$gene_name[subjectHits(nearest_hits)]
  ccre_df$nearest_tss_distance <- mcols(nearest_hits)$distance
  
  # Create annotation labels
  ccre_df$annotation_label <- paste0("Nearest gene = ", ccre_df$nearest_gene, 
                                     " (", ccre_df$nearest_tss_distance, " bp)")
  
  # Adjust label positions for cCRE annotations
  ccre_df$x_cCRE_label <- (ccre_df$start + ccre_df$end) / 2 - 20 * (ccre_df$end - ccre_df$start)
  ccre_df$x_nearest_label <- (ccre_df$start + ccre_df$end) / 2 + 150 * (ccre_df$end - ccre_df$start)
  ccre_df$y_nearest_label <- 0.05  
  
  # Define colors: red for the selected cCRE, green for others
  color_mapping <- c("selected" = "#4CBB17", "other" = "grey")
  
  # Gene Plot
  gene_plot <- ggplot(genes_df, aes(xmin = start, xmax = end, y = gene_name)) +
    geom_gene_arrow(fill = "#800020") +
    geom_text(aes(x = (start + end) / 2, label = gene_name), 
              position = position_nudge(y = 0.3),
              size = 3,
              color = "black") +
    theme_genes() +
    scale_x_continuous(limits = c(start(plot_range), end(plot_range))) +
    theme(
      axis.title.x = element_blank(),
      axis.text.x  = element_text(size = 10),
      axis.ticks.x = element_line(),  
      axis.title.y = element_blank(),
      axis.text.y  = element_blank(),
      axis.ticks.y = element_line(),
      plot.title = element_blank()
    )
  
  # cCRE Track Plot
  ccre_plot <- ggplot(ccre_df) +  
    geom_rect(aes(xmin = start, xmax = end, ymin = 0.45, ymax = 0.55, fill = highlight), 
              color = "black") +  
    scale_fill_manual(values = color_mapping) +  
    scale_y_continuous(limits = c(0.4, 0.6)) +
    ggtitle(bquote(atop(
      bold(.(ccre_df$cCRE_name[ccre_df$cCRE_name == ccre_id])) *  
        bold(" | Nearest gene = ") *  
        .(ccre_df$nearest_gene[ccre_df$cCRE_name == ccre_id]) *  
        bold(" | Distance = ") *  
        bold(.(ccre_df$nearest_tss_distance[ccre_df$cCRE_name == ccre_id])) * " bp"
    ))) +
    scale_x_continuous(limits = c(start(plot_range), end(plot_range))) +  
    theme_classic() +  
    theme(
      axis.title.x = element_blank(),
      axis.text.x  = element_text(size = 10),
      axis.ticks.x = element_line(),
      axis.title.y = element_blank(),
      axis.text.y  = element_blank(),
      axis.ticks.y = element_blank(),
      axis.line.y = element_blank(),  
      axis.line.x = element_line(),
      plot.title = element_text(hjust = 0.5, size = 12),
      legend.position = "none"  # Hide legend since colors are predefined
    )
  
  # Combine Plots
  combined_plot <- ccre_plot / gene_plot +  
    plot_layout(heights = c(0.5, 2)) +  
    plot_annotation(
      title = "",
      theme = theme(
        plot.title = element_text(hjust = 0.5, size = 14),
        plot.background = element_rect(color = "black", fill = NA, size = 1)  
      )
    )
  
  # Display the final plot
  print(combined_plot)
}



# Examples
plot_cCRE_region("ColoncCRE91407", flank = 50000, human_gtf_genes, human_ccre_filtered)


plot_cCRE_region("ColoncCRE20983", flank = 50000, human_gtf_genes, human_ccre_filtered)


plot_cCRE_region("ColoncCRE32335", flank = 35000, human_gtf_genes, human_ccre_filtered)

plot_cCRE_region("ColoncCRE14482", flank = 80000, human_gtf_genes, human_ccre_filtered)


plot_cCRE_region("ColoncCRE2324", flank = 120000, human_gtf_genes, human_ccre_filtered)

plot_cCRE_region("ColoncCRE6532", flank = 320000, human_gtf_genes, human_ccre_filtered)

plot_cCRE_region("ColoncCRE66693", flank = 62000, human_gtf_genes, human_ccre_filtered)

plot_cCRE_region("ColoncCRE73406", flank = 62000, human_gtf_genes, human_ccre_filtered)

plot_cCRE_region("ColoncCRE121299", flank = 178000, human_gtf_genes, human_ccre_filtered)

plot_cCRE_region("ColoncCRE121635", flank = 278000, human_gtf_genes, human_ccre_filtered)

plot_cCRE_region("ColoncCRE24316", flank = 278000, human_gtf_genes, human_ccre_filtered)


# Updated function
plot_cCRE_region <- function(ccre_id, flank = 50000, human_gtf_genes, human_ccre_filtered) { 
  # Subset the cCRE based on the given ID 
  ccre_selected <- human_ccre_filtered[human_ccre_filtered$V4 == ccre_id] 
  
  if (length(ccre_selected) == 0) { 
    stop("Error: cCRE ID not found in dataset.") 
  } 
  
  # Define the plot range based on flanking distance 
  plot_range <- GRanges( 
    seqnames = seqnames(ccre_selected), 
    ranges = IRanges( 
      start = start(ccre_selected) - flank, 
      end = end(ccre_selected) + flank 
    ) 
  ) 
  
  # Subset gene annotations and cCREs to the plotting region 
  human_genes_plot <- subsetByOverlaps(human_gtf_genes, plot_range) 
  human_ccre_plot  <- subsetByOverlaps(human_ccre_filtered, plot_range) 
  
  # Convert GRanges to data frames 
  genes_df <- as.data.frame(human_genes_plot) 
  ccre_df  <- as.data.frame(human_ccre_plot) 
  
  # Extract cCRE names 
  ccre_df$cCRE_name <- ccre_df$V4   
  
  # Identify the selected cCRE for color mapping 
  ccre_df$highlight <- ifelse(ccre_df$cCRE_name == ccre_id, "selected", "other") 
  
  # Compute the distance to the nearest TSS and annotate nearest gene 
  human_TSS <- promoters(human_gtf_genes, upstream = 0, downstream = 1) 
  nearest_hits <- distanceToNearest(human_ccre_plot, human_TSS) 
  
  # Extract gene names from human_TSS 
  ccre_df$nearest_gene <- mcols(human_TSS)$gene_name[subjectHits(nearest_hits)] 
  ccre_df$nearest_tss_distance <- mcols(nearest_hits)$distance 
  
  # Create annotation labels 
  ccre_df$annotation_label <- paste0("Nearest gene = ", ccre_df$nearest_gene, 
                                     " (", ccre_df$nearest_tss_distance, " bp)") 
  
  # Adjust label positions for cCRE annotations 
  ccre_df$x_cCRE_label <- (ccre_df$start + ccre_df$end) / 2 - 20 * (ccre_df$end - ccre_df$start) 
  ccre_df$x_nearest_label <- (ccre_df$start + ccre_df$end) / 2 + 150 * (ccre_df$end - ccre_df$start) 
  ccre_df$y_nearest_label <- 0.05 
  
  # Define colors: red for the selected cCRE, green for others 
  color_mapping <- c("selected" = "#4CBB17", "other" = "grey") 
  
  # Gene Plot: Adjust to show portion of the gene within the plot range
  gene_plot <- ggplot(genes_df, aes(xmin = pmax(start, start(plot_range)), 
                                    xmax = pmin(end, end(plot_range)), 
                                    y = gene_name)) + 
    # Positive strand
    geom_gene_arrow(data = subset(genes_df, strand == "+"), 
                    fill = "#800020") +
    # Negative strand 
    geom_gene_arrow(data = subset(genes_df, strand == "-"), 
                    fill = "#2171b5", forward = FALSE) +
    
    geom_text(aes(x = (pmax(start, start(plot_range)) + pmin(end, end(plot_range))) / 2, 
                  label = gene_name), 
              position = position_nudge(y = 0.3), 
              size = 3, 
              color = "black") + 
    theme_genes() + 
    scale_x_continuous(limits = c(start(plot_range), end(plot_range))) + 
    theme( 
      axis.title.x = element_blank(), 
      axis.text.x  = element_text(size = 10), 
      axis.ticks.x = element_line(),   
      axis.title.y = element_blank(), 
      axis.text.y  = element_blank(), 
      axis.ticks.y = element_line(), 
      plot.title = element_blank() 
    )
  
  # cCRE Track Plot with Highlighting Arrow (only for selected cCRE)
  ccre_plot <- ggplot(ccre_df) +   
    geom_rect(aes(xmin = start, xmax = end, ymin = 0.45, ymax = 0.55, fill = highlight),  
              color = "black") +   
    # Only add arrow for the selected cCRE
    geom_segment(data = subset(ccre_df, cCRE_name == ccre_id), 
                 aes(x = (start + end) / 2, xend = (start + end) / 2, 
                     y = 0.8, yend = 0.65),  # Adjusted y and yend for proper positioning
                 arrow = arrow(type = "closed", length = unit(0.1, "inches")), 
                 color = "red", size = 1) +  # Adjust arrow color and size
    scale_fill_manual(values = color_mapping) +   
    scale_y_continuous(limits = c(0.2, 0.8)) +  # Adjust y-axis limits to fit the downward arrow
    ggtitle(bquote(atop(
      bold(.(ccre_df$cCRE_name[ccre_df$cCRE_name == ccre_id])) ~ ", category = " ~ bold(.(ccre_df$V5[ccre_df$cCRE_name == ccre_id])),
      "Nearest gene = " ~ bold(.(ccre_df$nearest_gene[ccre_df$cCRE_name == ccre_id])) ~ ", distance = " ~ bold(.(ccre_df$nearest_tss_distance[ccre_df$cCRE_name == ccre_id])) ~ " bp"
    ))) + 
    scale_x_continuous(limits = c(start(plot_range), end(plot_range))) +   
    theme_classic() +   
    theme( 
      axis.title.x = element_blank(), 
      axis.text.x  = element_text(size = 10), 
      axis.ticks.x = element_line(), 
      axis.title.y = element_blank(), 
      axis.text.y  = element_blank(), 
      axis.ticks.y = element_blank(), 
      axis.line.y = element_blank(),   
      axis.line.x = element_line(), 
      plot.title = element_text(hjust = 0.5, size = 12), 
      legend.position = "none"  # Hide legend since colors are predefined 
    )
  
  # Combine Plots 
  combined_plot <- ccre_plot / gene_plot +   
    plot_layout(heights = c(0.5, 2)) +   
    plot_annotation( 
      title = "", 
      theme = theme( 
        plot.title = element_text(hjust = 0.5, size = 14), 
        plot.background = element_rect(color = "black", fill = NA, size = 1)   
      ) 
    )
  
  # Display the final plot 
  print(combined_plot)
}



plot_cCRE_region("ColoncCRE24316", flank = 278000, human_gtf_genes, human_ccre_filtered)
plot_cCRE_region("ColoncCRE91407", flank = 368000, human_gtf_genes, human_ccre_filtered)
plot_cCRE_region("ColoncCRE20983", flank = 50000, human_gtf_genes, human_ccre_filtered)
plot_cCRE_region("ColoncCRE32335", flank = 35000, human_gtf_genes, human_ccre_filtered)
plot_cCRE_region("ColoncCRE14482", flank = 80000, human_gtf_genes, human_ccre_filtered)



plot_cCRE_region("ColoncCRE91407", flank = 18900, human_gtf_genes, human_ccre_gr)

plot_cCRE_region("ColoncCRE91407", flank = 3048900, dog_gtf_genes, dog_ccre_gr)

plot_cCRE_region("ColoncCRE3958", flank = 98900, mouse_gtf_genes, mouse_ccre_gr)

plot_cCRE_region("ColoncCRE3958", flank = 98900, dog_gtf_genes, dog_ccre_gr)

plot_cCRE_region("ColoncCRE3958", flank = 98900, human_gtf_genes, human_ccre_gr)

# Promoters

plot_cCRE_region("ColoncCRE37", flank = 18900, human_gtf_genes, human_ccre_gr)
plot_cCRE_region("ColoncCRE37", flank = 18900, mouse_gtf_genes, mouse_ccre_gr)
plot_cCRE_region("ColoncCRE37", flank = 18900, dog_gtf_genes, dog_ccre_gr)

plot_cCRE_region("ColoncCRE189598", flank = 18900, human_gtf_genes, human_ccre_gr)
plot_cCRE_region("ColoncCRE189598", flank = 18900, mouse_gtf_genes, mouse_ccre_gr)
plot_cCRE_region("ColoncCRE189598", flank = 18900, dog_gtf_genes, dog_ccre_gr)










################################################################
##### UPDATED MARCH 2025
################################################################

### Load libraries
library(ggbio)
library(GenomicRanges)
library(ggplot2)
library(patchwork)

# Load the cCREs and annotation files for each species
## Human
human_ccre_bed = read.table("/Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/LiftOver/sc-only-lifted-cCREs-2025.bed")
human_ccre_gr <- makeGRangesFromDataFrame(human_ccre_bed, keep.extra.columns=TRUE,
                                          start.field="V2",
                                          end.field = "V3",
                                          seqnames.field="V1")
human_ccre_gr

human_gtf <- import("/Users/ryanhagan/NoCoSMiCC/Files/gencode.v47.annotation.gtf")
human_gtf_genes <- human_gtf[human_gtf$type == "gene" & human_gtf$gene_type == "protein_coding"]
human_genes_gr <- as(gtf_genes, "GRanges")
human_genes_gr

## Mouse
mouse_ccre_bed = read.table("/Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/LiftOver/sc_cCRE_mouse_lifted.bed")
mouse_ccre_gr <- makeGRangesFromDataFrame(mouse_ccre_bed, keep.extra.columns=TRUE,
                                          start.field="V2",
                                          end.field = "V3",
                                          seqnames.field="V1")
mouse_ccre_gr

mouse_gtf <- import("/Users/ryanhagan/NoCoSMiCC/TAD_Analysis/Files/gencode.vM10.annotation.gtf.gz")
mouse_gtf_genes <- mouse_gtf[mouse_gtf$type == "gene" & mouse_gtf$gene_type == "protein_coding"]
mouse_genes_gr <- as(gtf_genes, "GRanges")
mouse_genes_gr

## Dog
dog_ccre_bed = read.table("/Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/LiftOver/sc_cCRE_dog_lifted.bed")
dog_ccre_gr <- makeGRangesFromDataFrame(dog_ccre_bed, keep.extra.columns=TRUE,
                                        start.field="V2",
                                        end.field = "V3",
                                        seqnames.field="V1")
dog_ccre_gr

dog_gtf <- import("/Users/ryanhagan/NoCoSMiCC/TAD_Analysis/Files/Canis_lupus_familiarisboxer.Dog10K_Boxer_Tasha.113.gtf")
dog_gtf_genes <- dog_gtf[dog_gtf$type == "gene" & dog_gtf$gene_biotype == "protein_coding"]
dog_genes_gr <- as(gtf_genes, "GRanges")
dog_genes_gr

## modify seqnames of dog genome
current_seqlevels <- seqlevels(dog_gtf_genes)
new_seqlevels <- ifelse(grepl("^chr", current_seqlevels),
                        current_seqlevels,
                        paste0("chr", current_seqlevels))

# Set the new seqlevels to the GRanges object
seqlevels(dog_gtf_genes) <- new_seqlevels

## Define the function for plotting
plot_cCRE_region <- function(ccre_id, flank = 50000, human_gtf_genes, human_ccre_filtered, species = species) { 
  # Subset the cCRE based on the given ID 
  ccre_selected <- human_ccre_filtered[human_ccre_filtered$V4 == ccre_id] 
  
  if (length(ccre_selected) == 0) { 
    stop("Error: cCRE ID not found in dataset.") 
  } 
  
  # Define the plot range based on flanking distance 
  plot_range <- GRanges( 
    seqnames = seqnames(ccre_selected), 
    ranges = IRanges( 
      start = start(ccre_selected) - flank, 
      end = end(ccre_selected) + flank 
    ) 
  ) 
  
  # Subset gene annotations and cCREs to the plotting region 
  human_genes_plot <- subsetByOverlaps(human_gtf_genes, plot_range) 
  human_ccre_plot  <- subsetByOverlaps(human_ccre_filtered, plot_range) 
  
  # Convert GRanges to data frames 
  genes_df <- as.data.frame(human_genes_plot) 
  ccre_df  <- as.data.frame(human_ccre_plot) 
  
  # Extract cCRE names 
  ccre_df$cCRE_name <- ccre_df$V4   
  
  # Identify the selected cCRE for color mapping 
  ccre_df$highlight <- ifelse(ccre_df$cCRE_name == ccre_id, "selected", "other") 
  
  # Compute the distance to the nearest TSS and annotate nearest gene 
  human_TSS <- promoters(human_gtf_genes, upstream = 0, downstream = 1) 
  nearest_hits <- distanceToNearest(human_ccre_plot, human_TSS) 
  
  # Extract gene names from human_TSS 
  ccre_df$nearest_gene <- mcols(human_TSS)$gene_name[subjectHits(nearest_hits)] 
  ccre_df$nearest_tss_distance <- mcols(nearest_hits)$distance 
  
  # Create annotation labels 
  ccre_df$annotation_label <- paste0("Nearest gene = ", ccre_df$nearest_gene, 
                                     " (", ccre_df$nearest_tss_distance, " bp)") 
  
  # Adjust label positions for cCRE annotations 
  ccre_df$x_cCRE_label <- (ccre_df$start + ccre_df$end) / 2 - 20 * (ccre_df$end - ccre_df$start) 
  ccre_df$x_nearest_label <- (ccre_df$start + ccre_df$end) / 2 + 150 * (ccre_df$end - ccre_df$start) 
  ccre_df$y_nearest_label <- 0.05 
  
  # Define colors: red for the selected cCRE, green for others 
  color_mapping <- c("selected" = "#4CBB17", "other" = "grey") 
  
  # Gene Plot: Adjust to show portion of the gene within the plot range
  gene_plot <- ggplot(genes_df, aes(xmin = pmax(start, start(plot_range)), 
                                    xmax = pmin(end, end(plot_range)), 
                                    y = gene_name)) + 
    # Positive strand
    geom_gene_arrow(data = subset(genes_df, strand == "+"), 
                    fill = "#800020") +
    # Negative strand 
    geom_gene_arrow(data = subset(genes_df, strand == "-"), 
                    fill = "#2171b5", forward = FALSE) +
    
    geom_text(aes(x = (pmax(start, start(plot_range)) + pmin(end, end(plot_range))) / 2, 
                  label = gene_name), 
              position = position_nudge(y = 0.3), 
              size = 3, 
              color = "black") + 
    theme_genes() + 
    scale_x_continuous(limits = c(start(plot_range), end(plot_range))) + 
    labs(x = paste0("Region: ", as.character(seqnames(plot_range)),":", as.character(ranges(plot_range)))) +
    theme( 
      axis.title.x  = element_text(size = 10, face = "bold", vjust=-0.5),
      axis.text.x  = element_text(size = 10), 
      axis.ticks.x = element_line(),   
      axis.title.y = element_blank(), 
      axis.text.y  = element_blank(), 
      axis.ticks.y = element_line(), 
      plot.title = element_blank() 
    )
  
  # cCRE Track Plot with Highlighting Arrow (only for selected cCRE)
  ccre_plot <- ggplot(ccre_df) +   
    geom_rect(aes(xmin = start, xmax = end, ymin = 0.45, ymax = 0.55, fill = highlight),  
              color = "black") +   
    # Only add arrow for the selected cCRE
    geom_segment(data = subset(ccre_df, cCRE_name == ccre_id), 
                 aes(x = (start + end) / 2, xend = (start + end) / 2, 
                     y = 0.8, yend = 0.65),  # Adjusted y and yend for proper positioning
                 arrow = arrow(type = "closed", length = unit(0.1, "inches")), 
                 color = "red", size = 1) +  # Adjust arrow color and size
    scale_fill_manual(values = color_mapping) +   
    scale_y_continuous(limits = c(0.2, 0.8)) +  # Adjust y-axis limits to fit the downward arrow
    ggtitle(bquote(atop(
      bold(.(ccre_df$cCRE_name[ccre_df$cCRE_name == ccre_id])) ~ bold(", Category = ") ~ .(ccre_df$V5[ccre_df$cCRE_name == ccre_id]),
      bold("Nearest gene = ") ~ .(ccre_df$nearest_gene[ccre_df$cCRE_name == ccre_id]) ~ bold(", Distance = ") ~ bold(.(ccre_df$nearest_tss_distance[ccre_df$cCRE_name == ccre_id])) ~ " bp"
    ))) + 
    scale_x_continuous(limits = c(start(plot_range), end(plot_range))) +
    theme_classic() +   
    theme( 
     # axis.title.x = element_blank(), 
      axis.text.x  = element_text(size = 10), 
      axis.ticks.x = element_line(),
      axis.title.y = element_blank(),
      axis.title.x  = element_blank(),
      axis.text.y  = element_blank(), 
      axis.ticks.y = element_blank(), 
      axis.line.y = element_blank(),   
      axis.line.x = element_line(), 
      plot.title = element_text(hjust = 0.5, size = 12), 
      legend.position = "none"
    )
  
  # Combine Plots 
  combined_plot <- ccre_plot / gene_plot +   
    plot_layout(heights = c(0.5, 2)) +   
    plot_annotation( 
      title = species, 
      theme = theme( 
        plot.title = element_text(hjust = 0.5, size = 18, face = "bold", colour = "#7F3D17"), 
        plot.background = element_rect(color = "black", fill = NA, size = 1)   
      ) 
    )
  
  # Display the final plot 
  print(combined_plot)
}

human_cCRE <- plot_cCRE_region("ColoncCRE42", flank = 6800, human_gtf_genes, human_ccre_gr, species = "Human")
human_cCRE
mouse_cCRE <- plot_cCRE_region("ColoncCRE42", flank = 6800, mouse_gtf_genes, mouse_ccre_gr, species = "Mouse")
mouse_cCRE
dog_cCRE <- plot_cCRE_region("ColoncCRE42", flank = 6800, dog_gtf_genes, dog_ccre_gr, species = "Dog")
dog_cCRE

combined_plot <- wrap_elements(human_cCRE) | wrap_elements(mouse_cCRE) | wrap_elements(dog_cCRE)
print(combined_plot)

# EPCAM enhancer
plot_cCRE_region("ColoncCRE97443", flank = 10000, human_gtf_genes, human_ccre_gr, species = "Human")
plot_cCRE_region("ColoncCRE97443", flank = 40000, mouse_gtf_genes, mouse_ccre_gr, species = "Mouse")
plot_cCRE_region("ColoncCRE97443", flank = 40000, dog_gtf_genes, dog_ccre_gr, species = "Dog")

# PDZK1IP1 enhancer - activated in epithelial cancer cells (paper)
plot_cCRE_region("ColoncCRE5636", flank = 40000, human_gtf_genes, human_ccre_gr, species = "Human")
plot_cCRE_region("ColoncCRE5636", flank = 40000, mouse_gtf_genes, mouse_ccre_gr, species = "Mouse")
plot_cCRE_region("ColoncCRE5636", flank = 40000, dog_gtf_genes, dog_ccre_gr, species = "Dog")









### Plot proportions of cCREs across species with distance/ category conserved

library(ggplot2)
library(dplyr)
library(RColorBrewer)

# extract transcripts for each species
human_gtf_transcripts <- human_gtf[human_gtf$type == "transcript"]
mouse_gtf_transcripts <- mouse_gtf[mouse_gtf$type == "transcript"]
dog_gtf_transcripts <- dog_gtf[dog_gtf$type == "transcript"]
## modify seqnames of dog genome
current_seqlevels <- seqlevels(dog_gtf_transcripts)
new_seqlevels <- ifelse(grepl("^chr", current_seqlevels),
                        current_seqlevels,
                        paste0("chr", current_seqlevels))
seqlevels(dog_gtf_transcripts) <- new_seqlevels

# Function to calculate distances and categorize them
get_distance_categories <- function(ccre_gr, gene_gr, species_name) {
  # Calculate distance to nearest gene
  nearest_hits <- distanceToNearest(ccre_gr, gene_gr)
  distances <- mcols(nearest_hits)$distance
  
  # Create a data frame with cCRE IDs and distance categories
  distance_category <- ifelse(distances == 0, "Distance = 0",
                              ifelse(distances <= 1000, "Distance <= 1000bp", "Distance > 1000bp"))
  
  # Return a data frame with cCRE IDs and their distance categories
  data.frame(
    cCRE_ID = mcols(ccre_gr)$V4,  # Ensure we are using correct column name for cCRE IDs
    Distance_Category = distance_category,
    Species = species_name
  )
}

# Get the distance categories for each species (Human, Mouse, Dog)
human_data <- get_distance_categories(human_ccre_gr, human_gtf_transcripts, "Human")

# Check if there are any human cCREs with distance = 0
human_zero_distances <- human_data %>% filter(Distance_Category == "Distance = 0")

# Filter mouse and dog data based on the same cCRE IDs
mouse_data <- get_distance_categories(mouse_ccre_gr, mouse_gtf_transcripts, "Mouse")
dog_data <- get_distance_categories(dog_ccre_gr, dog_gtf_transcripts, "Dog")

# Filter mouse and dog data to only include the cCREs that were 0 in human
mouse_filtered <- mouse_data %>% filter(cCRE_ID %in% human_zero_distances$cCRE_ID)
dog_filtered <- dog_data %>% filter(cCRE_ID %in% human_zero_distances$cCRE_ID)

# Combine all the data (Human, Mouse, Dog) into a single data frame for plotting
combined_data <- bind_rows(
  data.frame(cCRE_ID = human_zero_distances$cCRE_ID, Distance_Category = "Distance = 0", Species = "Human"),
  mouse_filtered,
  dog_filtered
)

# Reorder the levels of Distance_Category to get the desired stacking order
combined_data$Distance_Category <- factor(combined_data$Distance_Category, 
                                          levels = c("Distance = 0", "Distance <= 1000bp", "Distance > 1000bp"))

# Reorder the levels of Species for the desired order of "Human", "Mouse", "Dog"
combined_data$Species <- factor(combined_data$Species, levels = c("Human", "Mouse", "Dog"))

# Choose a color palette (I am using a brewer palette, but feel free to choose your own)
color_palette <- c("Distance = 0" = "#66c2a5", 
                   "Distance <= 1000bp" = "#fc8d62", 
                   "Distance > 1000bp" = "#8da0cb")

# Plot a stacked bar chart with percentage on y-axis
ggplot(combined_data, aes(x = Species, fill = Distance_Category)) +
  geom_bar(position = "fill", width = 0.7) +  # Use position="fill" to show fractions (percentages)
  scale_fill_manual(values = color_palette) +
  scale_y_continuous(labels = scales::percent) +  # Convert y-axis to percentage
  labs(
    title = "Fraction of cCREs with Distance Categories for Human, Mouse, and Dog",
    x = "Species",
    y = "Percentage",
    fill = "Distance Category"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 14, face = "bold"),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    legend.position = "top"
  )

# Reverse the factor levels for Distance_Category
combined_data$Distance_Category <- factor(combined_data$Distance_Category, 
                                          levels = c("Distance > 1000bp", "Distance <= 1000bp", "Distance = 0"))

# Plot again with the new ordering
ggplot(combined_data, aes(x = Species, fill = Distance_Category)) +
  geom_bar(position = "fill", width = 0.7) +  # Use position="fill" to show fractions (percentages)
  scale_fill_manual(values = color_palette) +
  scale_y_continuous(labels = scales::percent) +  # Convert y-axis to percentage
  labs(
    title = "Fraction of cCREs with Distance Categories for Human, Mouse, and Dog",
    x = "Species",
    y = "Percentage",
    fill = "Distance Category"
  ) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 14, face = "bold"),
    plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
    legend.position = "top"
  )


####

library(rtracklayer)
library(GenomicRanges)
library(ggplot2)
library(dplyr)

# Human
# 1. Subset the GTF to protein-coding transcripts
human_gtf_transcripts <- human_gtf[human_gtf$type == "transcript"]

# 2. Compute TSS positions from the transcript annotations
human_TSS <- promoters(human_gtf_transcripts, upstream = 0, downstream = 1)

# 3. Calculate the distance from each cCRE in human_ccre_gr to its nearest TSS
nearest_hits <- distanceToNearest(human_ccre_gr, human_TSS)

# 4. Convert your cCRE GRanges to a data frame and add the distance information
ccre_df <- as.data.frame(human_ccre_gr)
ccre_df$nearest_tss_distance <- NA
ccre_df$nearest_tss_distance[queryHits(nearest_hits)] <- mcols(nearest_hits)$distance

# 5. Group by cCRE category (stored in column 'V5') and compute the mean distance, ordering by mean distance
mean_distance <- ccre_df %>%
  group_by(V5) %>%
  summarise(mean_distance = mean(nearest_tss_distance, na.rm = TRUE)) %>%
  arrange(mean_distance)

# 6. Plot the mean distances by cCRE category with mean values shown in the middle of each bar
ggplot(mean_distance, aes(x = reorder(V5, mean_distance), y = mean_distance, fill = V5)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = round(mean_distance, 0)), vjust = 1.5, color = "black", size = 5) +
  scale_fill_brewer(palette = "Set3") +
  labs(x = "cCRE Category", y = "Mean distance to nearest TSS (bp)",
       title = "Mean Distance to Nearest TSS by cCRE Category") +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14, face = "bold"),
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5))

# Mouse
# 1. Subset the GTF to transcripts
mouse_gtf_transcripts <- mouse_gtf[mouse_gtf$type == "transcript"]

# 2. Compute TSS positions from the transcript annotations
mouse_TSS <- promoters(mouse_gtf_transcripts, upstream = 0, downstream = 1)

# 3. Calculate the distance from each cCRE in human_ccre_gr to its nearest TSS
nearest_hits <- distanceToNearest(mouse_ccre_gr, mouse_TSS)

# 4. Convert your cCRE GRanges to a data frame and add the distance information
ccre_df <- as.data.frame(mouse_ccre_gr)
ccre_df$nearest_tss_distance <- NA
ccre_df$nearest_tss_distance[queryHits(nearest_hits)] <- mcols(nearest_hits)$distance

# 5. Group by cCRE category (stored in column 'V5') and compute the mean distance, ordering by mean distance
mean_distance <- ccre_df %>%
  group_by(V5) %>%
  summarise(mean_distance = mean(nearest_tss_distance, na.rm = TRUE)) %>%
  arrange(mean_distance)

# 6. Plot the mean distances by cCRE category with mean values shown in the middle of each bar
ggplot(mean_distance, aes(x = reorder(V5, mean_distance), y = mean_distance, fill = V5)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = round(mean_distance, 0)), vjust = 1.5, color = "black", size = 5) +
  scale_fill_brewer(palette = "Set3") +
  labs(x = "cCRE Category", y = "Mean distance to nearest TSS (bp)",
       title = "Mean Distance to Nearest TSS by cCRE Category") +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14, face = "bold"),
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5))

# how many PLS cCREs are within 200bp of a mouse TSS
# Filter for PLS or PLS,CTCF-bound cCREs
target_ccres <- ccre_df %>% 
  filter(V5 %in% c("PLS", "PLS,CTCF-bound"))

# Count how many of these have a nearest TSS distance of 200bp or less
num_within_200bp <- sum(target_ccres$nearest_tss_distance <= 200, na.rm = TRUE)
total_target_ccres <- nrow(target_ccres)

cat("Number of PLS or PLS,CTCF-bound cCREs within 200bp of a TSS:", num_within_200bp, "\n")
cat("Total number of PLS or PLS,CTCF-bound cCREs:", total_target_ccres, "\n")

# Dog
# 1. Subset the GTF to protein-coding transcripts
dog_gtf_transcripts <- dog_gtf[dog_gtf$type == "transcript"]
## modify seqnames of dog genome
current_seqlevels <- seqlevels(dog_gtf_transcripts)
new_seqlevels <- ifelse(grepl("^chr", current_seqlevels),
                        current_seqlevels,
                        paste0("chr", current_seqlevels))
seqlevels(dog_gtf_transcripts) <- new_seqlevels

# 2. Compute TSS positions from the transcript annotations
dog_TSS <- promoters(dog_gtf_transcripts, upstream = 0, downstream = 1)

# 3. Calculate the distance from each cCRE in human_ccre_gr to its nearest TSS
nearest_hits <- distanceToNearest(dog_ccre_gr, dog_TSS)

# 4. Convert your cCRE GRanges to a data frame and add the distance information
ccre_df <- as.data.frame(dog_ccre_gr)
ccre_df$nearest_tss_distance <- NA
ccre_df$nearest_tss_distance[queryHits(nearest_hits)] <- mcols(nearest_hits)$distance

# 5. Group by cCRE category (stored in column 'V5') and compute the mean distance, ordering by mean distance
mean_distance <- ccre_df %>%
  group_by(V5) %>%
  summarise(mean_distance = mean(nearest_tss_distance, na.rm = TRUE)) %>%
  arrange(mean_distance)

# 6. Plot the mean distances by cCRE category with mean values shown in the middle of each bar
ggplot(mean_distance, aes(x = reorder(V5, mean_distance), y = mean_distance, fill = V5)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = round(mean_distance, 0)), vjust = 1.5, color = "black", size = 5) +
  scale_fill_brewer(palette = "Set3") +
  labs(x = "cCRE Category", y = "Mean distance to nearest TSS (bp)",
       title = "Mean Distance to Nearest TSS by cCRE Category") +
  theme_minimal() +
  theme(legend.position = "none",
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 14, face = "bold"),
        plot.title = element_text(size = 16, face = "bold", hjust = 0.5))

# how many PLS cCREs are within 200bp of a dog TSS
# Filter for PLS or PLS,CTCF-bound cCREs
target_ccres <- ccre_df %>% 
  filter(V5 %in% c("PLS", "PLS,CTCF-bound"))

# Count how many of these have a nearest TSS distance of 200bp or less
num_within_200bp <- sum(target_ccres$nearest_tss_distance <= 200, na.rm = TRUE)
total_target_ccres <- nrow(target_ccres)

cat("Number of PLS or PLS,CTCF-bound cCREs within 200bp of a TSS:", num_within_200bp, "\n")
cat("Total number of PLS or PLS,CTCF-bound cCREs:", total_target_ccres, "\n")

####

# Look at the nearest genes for cCREs in each species

# Human
human_nearest_hits <- distanceToNearest(human_ccre_gr, human_TSS)
human_ccre_df <- as.data.frame(human_ccre_gr)
human_ccre_df$nearest_tss_distance <- NA
human_ccre_df$nearest_tss_distance[queryHits(human_nearest_hits)] <- mcols(human_TSS)$gene_name[subjectHits(human_nearest_hits)]

names(human_ccre_df)[names(human_ccre_df) == "nearest_tss_distance"] <- "nearest_gene_human"

# Mouse
mouse_nearest_hits <- distanceToNearest(mouse_ccre_gr, mouse_TSS)
mouse_ccre_df <- as.data.frame(mouse_ccre_gr)
mouse_ccre_df$nearest_tss_distance <- NA
mouse_ccre_df$nearest_tss_distance[queryHits(mouse_nearest_hits)] <- mcols(mouse_TSS)$gene_name[subjectHits(mouse_nearest_hits)]

names(mouse_ccre_df)[names(mouse_ccre_df) == "nearest_tss_distance"] <- "nearest_gene_mouse"

# Dog
dog_nearest_hits <- distanceToNearest(dog_ccre_gr, dog_TSS)
dog_ccre_df <- as.data.frame(dog_ccre_gr)
dog_ccre_df$nearest_tss_distance <- NA
dog_ccre_df$nearest_tss_distance[queryHits(dog_nearest_hits)] <- mcols(dog_TSS)$gene_name[subjectHits(dog_nearest_hits)]

names(dog_ccre_df)[names(dog_ccre_df) == "nearest_tss_distance"] <- "nearest_gene_dog"

# 
merged_df <- merge(human_ccre_df[, c("V4", "nearest_gene_human")],
                   mouse_ccre_df[, c("V4", "nearest_gene_mouse")],
                   by = "V4")
merged_df <- merge(merged_df,
                   dog_ccre_df[, c("V4", "nearest_gene_dog")],
                   by = "V4")
head(merged_df)



merged_df <- merge(human_ccre_df[, c("V4", "nearest_gene_human")],
                   dog_ccre_df[, c("V4", "nearest_gene_dog")],
                   by = "V4")
head(merged_df)



plot_cCRE_region("ColoncCRE101760", flank = 10000, human_gtf_transcripts, human_ccre_gr, species = "Human")
plot_cCRE_region("ColoncCRE101760", flank = 10000, mouse_gtf_transcripts, mouse_ccre_gr, species = "Mouse")
plot_cCRE_region("ColoncCRE101760", flank = 10000, dog_gtf_transcripts, dog_ccre_gr, species = "Dog")




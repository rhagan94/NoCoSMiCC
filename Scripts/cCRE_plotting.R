################################################################
##### Visualising colon cCREs
################################################################

### Load libraries
library(ggbio)
library(GenomicRanges)
library(ggplot2)
library(patchwork)
library(dplyr)
library(RColorBrewer)
library(rtracklayer)
library(GenomeInfoDb)
library(gggenes)

# Load the cCREs and annotation files for each species
## Human
human_ccre_bed = read.table("/project/home/p201120/ryan/cCRE_pipeline/outputs/LiftOver/colon-epithelial-lifted-cCREs.bed")
human_ccre_gr <- makeGRangesFromDataFrame(human_ccre_bed, keep.extra.columns=TRUE,
                                          start.field="V2",
                                          end.field = "V3",
                                          seqnames.field="V1")
human_ccre_gr

human_gtf <- import("/mnt/tier2/project/p201120/ryan/refs/gtf/Homo_sapiens.GRCh38.113.gtf.gz")
human_gtf_genes <- human_gtf[human_gtf$type == "gene" & human_gtf$gene_biotype == "protein_coding"]
human_genes_gr <- as(human_gtf_genes, "GRanges")
human_genes_gr

## Mouse
mouse_ccre_bed = read.table("/project/home/p201120/ryan/cCRE_pipeline/outputs/LiftOver/sc_cCRE_mouse_lifted_0.5.bed")
mouse_ccre_gr <- makeGRangesFromDataFrame(mouse_ccre_bed, keep.extra.columns=TRUE,
                                          start.field="V2",
                                          end.field = "V3",
                                          seqnames.field="V1")
mouse_ccre_gr

mouse_gtf <- import("/mnt/tier2/project/p201120/ryan/refs/gtf/Mus_musculus.GRCm39.113.gtf.gz")
mouse_gtf_genes <- mouse_gtf[mouse_gtf$type == "gene" & mouse_gtf$gene_biotype == "protein_coding"]
mouse_genes_gr <- as(mouse_gtf_genes, "GRanges")
mouse_genes_gr

## Dog
dog_ccre_bed = read.table("/project/home/p201120/ryan/cCRE_pipeline/outputs/LiftOver/sc_cCRE_dog_lifted_0.5.bed")
dog_ccre_gr <- makeGRangesFromDataFrame(dog_ccre_bed, keep.extra.columns=TRUE,
                                        start.field="V2",
                                        end.field = "V3",
                                        seqnames.field="V1")
dog_ccre_gr

dog_gtf <- import("/mnt/tier2/project/p201120/ryan/refs/gtf/Canis_lupus_familiarisboxer.Dog10K_Boxer_Tasha.113.gtf.gz")
dog_gtf_genes <- dog_gtf[dog_gtf$type == "gene" & dog_gtf$gene_biotype == "protein_coding"]
dog_genes_gr <- as(dog_gtf_genes, "GRanges")
dog_genes_gr

## modify seqnames
current_seqlevels_dog <- seqlevels(dog_gtf_genes)
new_seqlevels <- ifelse(grepl("^chr", current_seqlevels_dog),
                        current_seqlevels,
                        paste0("chr", current_seqlevels_dog))

# Set the new seqlevels to the GRanges object
seqlevels(dog_gtf_genes) <- new_seqlevels

current_seqlevels_human <- seqlevels(human_gtf_genes)
new_seqlevels <- ifelse(grepl("^chr", current_seqlevels_human),
                        current_seqlevels,
                        paste0("chr", current_seqlevels_human))

# Set the new seqlevels to the GRanges object
seqlevels(human_gtf_genes) <- new_seqlevels

current_seqlevels_mouse <- seqlevels(mouse_gtf_genes)
new_seqlevels <- ifelse(grepl("^chr", current_seqlevels_mouse),
                        current_seqlevels,
                        paste0("chr", current_seqlevels_mouse))

# Set the new seqlevels to the GRanges object
seqlevels(mouse_gtf_genes) <- new_seqlevels

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

human_cCRE <- plot_cCRE_region("ColoncCRE71", flank = 6800, human_gtf_genes, human_ccre_gr, species = "Human")
human_cCRE
mouse_cCRE <- plot_cCRE_region("ColoncCRE71", flank = 6800, mouse_gtf_genes, mouse_ccre_gr, species = "Mouse")
mouse_cCRE
dog_cCRE <- plot_cCRE_region("ColoncCRE71", flank = 6800, dog_gtf_genes, dog_ccre_gr, species = "Dog")
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

####################################################################################
# Plot proportions of cCREs across species with distance/ category conserved
####################################################################################

# Human
# 
human_gtf_transcripts <- human_gtf[human_gtf$type == "transcript"]

# 
human_TSS <- promoters(human_gtf_transcripts, upstream = 0, downstream = 1)

# 
nearest_hits <- distanceToNearest(human_ccre_gr, human_TSS)

# 
ccre_df <- as.data.frame(human_ccre_gr)
ccre_df$nearest_tss_distance <- NA
ccre_df$nearest_tss_distance[queryHits(nearest_hits)] <- mcols(nearest_hits)$distance

# 
mean_distance <- ccre_df %>%
  group_by(V5) %>%
  summarise(mean_distance = mean(nearest_tss_distance, na.rm = TRUE)) %>%
  arrange(mean_distance)

# 
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
# 
mouse_gtf_transcripts <- mouse_gtf[mouse_gtf$type == "transcript"]

#
mouse_TSS <- promoters(mouse_gtf_transcripts, upstream = 0, downstream = 1)

#
nearest_hits <- distanceToNearest(mouse_ccre_gr, mouse_TSS)

#
ccre_df <- as.data.frame(mouse_ccre_gr)
ccre_df$nearest_tss_distance <- NA
ccre_df$nearest_tss_distance[queryHits(nearest_hits)] <- mcols(nearest_hits)$distance

#
mean_distance <- ccre_df %>%
  group_by(V5) %>%
  summarise(mean_distance = mean(nearest_tss_distance, na.rm = TRUE)) %>%
  arrange(mean_distance)

#
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
#
dog_gtf_transcripts <- dog_gtf[dog_gtf$type == "transcript"]
## modify seqnames of dog genome
current_seqlevels <- seqlevels(dog_gtf_transcripts)
new_seqlevels <- ifelse(grepl("^chr", current_seqlevels),
                        current_seqlevels,
                        paste0("chr", current_seqlevels))
seqlevels(dog_gtf_transcripts) <- new_seqlevels

#
dog_TSS <- promoters(dog_gtf_transcripts, upstream = 0, downstream = 1)

#
nearest_hits <- distanceToNearest(dog_ccre_gr, dog_TSS)

#
ccre_df <- as.data.frame(dog_ccre_gr)
ccre_df$nearest_tss_distance <- NA
ccre_df$nearest_tss_distance[queryHits(nearest_hits)] <- mcols(nearest_hits)$distance

#
mean_distance <- ccre_df %>%
  group_by(V5) %>%
  summarise(mean_distance = mean(nearest_tss_distance, na.rm = TRUE)) %>%
  arrange(mean_distance)

#
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


###### 2026 Update #######
# =============================================================================
# cCRE synteny plot (horizontal 3-species panel)
# Final reproducible version used for ColoncCRE189745
# =============================================================================

suppressPackageStartupMessages({
  library(GenomicRanges)
  library(ggplot2)
  library(gggenes)
  library(patchwork)
  library(grid)
  library(scales)
})

plot_cCRE_synteny_horizontal <- function(ccre_id,
                                         flank_left    = 10000,
                                         flank_right   = 120000,
                                         human_gtf     = human_gtf_genes,
                                         mouse_gtf     = mouse_gtf_genes,
                                         dog_gtf       = dog_gtf_genes,
                                         human_ccre    = human_ccre_gr,
                                         mouse_ccre    = mouse_ccre_gr,
                                         dog_ccre      = dog_ccre_gr) {

  color_mapping  <- c("selected" = "#4CBB17", "other" = "grey70")
  species_colors <- c("Human" = "#7F3D17", "Mouse" = "#1A5276", "Dog" = "#1D6A39")

  make_species_panel <- function(ccre_gr, gtf_genes, species_label) {

    ccre_selected <- ccre_gr[mcols(ccre_gr)$name == ccre_id]
    if (length(ccre_selected) == 0) {
      message(species_label, ": cCRE not found — skipping.")
      return(NULL)
    }

    p_start <- max(1L, as.integer(start(ccre_selected)) - flank_left)
    p_end   <- as.integer(end(ccre_selected)) + flank_right
    p_chr   <- as.character(seqnames(ccre_selected))

    plot_range <- GRanges(
      seqnames = p_chr,
      ranges   = IRanges(start = p_start, end = p_end)
    )

    genes_sub <- subsetByOverlaps(gtf_genes, plot_range)
    ccre_sub  <- subsetByOverlaps(ccre_gr,   plot_range)

    genes_df <- data.frame(
      start     = pmax(as.integer(start(genes_sub)), p_start),
      end       = pmin(as.integer(end(genes_sub)),   p_end),
      strand    = as.character(strand(genes_sub)),
      gene_name = as.character(mcols(genes_sub)$gene_name),
      stringsAsFactors = FALSE
    )
    genes_df <- genes_df[!is.na(genes_df$gene_name), ]
    genes_df <- genes_df[!duplicated(genes_df$gene_name), ]

    ccre_df <- data.frame(
      start = as.integer(start(ccre_sub)),
      end   = as.integer(end(ccre_sub)),
      name  = as.character(mcols(ccre_sub)$name),
      type  = as.character(mcols(ccre_sub)$type),
      stringsAsFactors = FALSE
    )
    ccre_df$highlight <- ifelse(ccre_df$name == ccre_id, "selected", "other")

    sel <- ccre_df[ccre_df$name == ccre_id, ]
    region_label <- paste0(
      p_chr, ": ",
      format(p_start, big.mark = ","),
      " - ",
      format(p_end, big.mark = ",")
    )
    sp_col <- species_colors[[species_label]]

    # -------------------------------------------------------------------------
    # cCRE track
    # -------------------------------------------------------------------------
    ccre_track <- ggplot(ccre_df) +
      geom_rect(
        aes(xmin = start, xmax = end, ymin = 0.45, ymax = 0.55, fill = highlight),
        colour = "black",
        linewidth = 0.3
      ) +
      geom_segment(
        data = sel,
        aes(
          x    = (start + end) / 2,
          xend = (start + end) / 2,
          y    = 0.78,
          yend = 0.63
        ),
        arrow     = arrow(type = "closed", length = unit(0.08, "inches")),
        colour    = "red",
        linewidth = 0.8
      ) +
      scale_fill_manual(values = color_mapping) +
      scale_x_continuous(
        limits = c(p_start, p_end),
        labels = scales::comma
      ) +
      scale_y_continuous(limits = c(0.2, 0.85)) +
      labs(title = species_label) +
      theme_classic() +
      theme(
        plot.title      = element_text(
          size   = 15,
          face   = "bold",
          colour = sp_col,
          hjust  = 0.5
        ),
        axis.text.x     = element_blank(),
        axis.ticks.x    = element_blank(),
        axis.title.x    = element_blank(),
        axis.line.x     = element_blank(),
        axis.title.y    = element_blank(),
        axis.text.y     = element_blank(),
        axis.ticks.y    = element_blank(),
        axis.line.y     = element_blank(),
        legend.position = "none",
        plot.margin     = margin(4, 6, 0, 6)
      ) +
      coord_cartesian(clip = "off")

    # -------------------------------------------------------------------------
    # Gene track
    # -------------------------------------------------------------------------
    gene_track <- ggplot(
      genes_df,
      aes(xmin = start, xmax = end, y = gene_name)
    ) +
      geom_gene_arrow(
        data      = subset(genes_df, strand == "+"),
        fill      = "#800020",
        colour    = "grey30",
        linewidth = 0.2
      ) +
      geom_gene_arrow(
        data      = subset(genes_df, strand == "-"),
        fill      = "#2171b5",
        colour    = "grey30",
        linewidth = 0.2,
        forward   = FALSE
      ) +
      geom_text(
        aes(x = (start + end) / 2, label = gene_name),
        position = position_nudge(y = 0.35),
        size     = 3.5,
        colour   = "grey20"
      ) +
      theme_genes() +
      scale_x_continuous(
        limits = c(p_start, p_end),
        labels = scales::comma,
        guide  = guide_axis(check.overlap = TRUE)
      ) +
      labs(x = region_label) +
      theme(
        axis.title.x = element_text(
          size   = 10,
          face   = "bold",
          colour = "grey40",
          vjust  = -0.3
        ),
        axis.text.x  = element_text(
          size   = 9,
          colour = "grey40",
          margin = margin(t = 6)
        ),
        axis.ticks.x = element_line(colour = "grey70"),
        axis.title.y = element_blank(),
        axis.text.y  = element_blank(),
        axis.ticks.y = element_blank(),
        plot.margin  = margin(0, 6, 18, 6)
      ) +
      coord_cartesian(clip = "off")

    # -------------------------------------------------------------------------
    # Stack cCRE + gene track, wrap in species-coloured box
    # -------------------------------------------------------------------------
    inner <- ccre_track / gene_track +
      plot_layout(heights = c(0.4, 1))

    inner + plot_annotation(
      theme = theme(
        plot.background = element_rect(
          colour    = sp_col,
          fill      = "white",
          linewidth = 1.2
        ),
        plot.margin = margin(10, 10, 14, 10)
      )
    )
  }

  # ---------------------------------------------------------------------------
  # Build species panels
  # ---------------------------------------------------------------------------
  human_panel <- make_species_panel(human_ccre, human_gtf, "Human")
  mouse_panel <- make_species_panel(mouse_ccre, mouse_gtf, "Mouse")
  dog_panel   <- make_species_panel(dog_ccre,   dog_gtf,   "Dog")

  cat_label <- as.character(
    mcols(human_ccre[mcols(human_ccre)$name == ccre_id])$type
  )

  # ---------------------------------------------------------------------------
  # Assemble final horizontal figure
  # ---------------------------------------------------------------------------
  combined <- (
    wrap_elements(human_panel) |
    wrap_elements(mouse_panel) |
    wrap_elements(dog_panel)
  ) +
    plot_annotation(
      title = paste0(ccre_id, "  |  Category: ", cat_label),
      theme = theme(
        plot.title = element_text(
          hjust  = 0.5,
          size   = 18,
          face   = "bold",
          colour = "grey15"
        ),
        plot.margin = margin(15, 15, 15, 15)
      )
    )

  print(combined)
  invisible(combined)
}

# =============================================================================
# Example run
# =============================================================================

p <- plot_cCRE_synteny_horizontal(
  "ColoncCRE189745",
  flank_left  = 10000,
  flank_right = 120000
)

ggsave(
  file.path(wd, "results/cCRE189745_synteny_horizontal.png"),
  plot   = p,
  width  = 12,
  height = 6,
  dpi    = 300
)

ggsave(
  file.path(wd, "results/cCRE189745_synteny_horizontal.pdf"),
  plot   = p,
  width  = 12,
  height = 6,
  device = cairo_pdf
)

message("Saved!")


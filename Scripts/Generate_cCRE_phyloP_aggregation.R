
# -- Kaili Fan
# fankaili.bio@gmail.com
# Weng Lab, UMass Chan Med
# Jun, 2022

# This script is for making aggregation plot for given cCRE phyloP matrix.

library(ggplot2)

#args = commandArgs(trailingOnly=TRUE)
#phyloP100_matrix = args[1]
phyloP241_matrix = "/Users/ryanhagan/NoCoSMiCC/Zoonomia/Outputs/phyloP241_matrix.mtx"
#cCRE_file = "/Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/hg38-rCARs/MaxZ/hg38-cCREs-v3.bed"
#cCRE_file = "/Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/LiftOver/cCRE_lifted"
#cCRE_file = "/Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/LiftOver/cCRE_dog_lifted"
#cCRE_file = "/Users/ryanhagan/NoCoSMiCC/Synteny/syntenic_cCREs.bed"
#cCRE_file = "/Users/ryanhagan/NoCoSMiCC/Synteny/nonsyntenic_cCREs.bed"

#cCRE_file = "/Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/LiftOver/cCRE_unlifted"
#cCRE_file = "/Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/LiftOver/cCRE_dog_unlifted"

#cCRE_file = "/Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/hg38-rCARs/MaxZ/liftover-syntenic-cCREs.bed"

cCRE_file = "/Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/LiftOver/cCRE_mouse_lifted_strict"

cCRE_file = "/Users/ryanhagan/NoCoSMiCC/colon_TAD_overlaps.bed"

cCRE_file = "/Users/ryanhagan/NoCoSMiCC/TAD_Analysis/Outputs/cCRE_domains_overlap.bed"
cCRE_file = "/Users/ryanhagan/NoCoSMiCC/TAD_Analysis/Outputs/cCRE_domains_no_overlap.bed"

cCRE_file = "/Users/ryanhagan/NoCoSMiCC/TAD_Analysis/Outputs/G1_cCREs_in_TADs.bed"
cCRE_file = "/Users/ryanhagan/NoCoSMiCC/TAD_Analysis/Outputs/G1_cCREs_NOT_in_TADs.bed"

cCRE_file = "/Users/ryanhagan/NoCoSMiCC/Zoonomia/Conserved_colon_G1_cCREs.bed"
cCRE_file = "/Users/ryanhagan/NoCoSMiCC/Zoonomia/Evolving_colon_G2_cCREs.bed"
cCRE_file = "/Users/ryanhagan/NoCoSMiCC/Zoonomia/Primate_colon_G3_cCREs.bed"

cCRE_file = "/Users/ryanhagan/Downloads/exon_cCRE_overlap.bed"
cCRE_file = "/Users/ryanhagan/Downloads/exon_cCRE_complete_overlap.bed"
cCRE_file = "/Users/ryanhagan/Downloads/exon_cCRE_no_overlap.bed"

cCRE_file = "/Users/ryanhagan/NoCoSMiCC/HiC_Analysis/compartment_A_cCREs.bed"
cCRE_file = "/Users/ryanhagan/NoCoSMiCC/HiC_Analysis/compartment_B_cCREs.bed"
cCRE_file = "/Users/ryanhagan/NoCoSMiCC/HiC_Analysis/compartment_I_cCREs.bed"

cCRE_file = "/Users/ryanhagan/NoCoSMiCC/Zoonomia/hg38_random_10k_250bp.bed"

## new phase of the analysis - separate single cell and bulk
cCRE_file = "/Users/ryanhagan/NoCoSMiCC/HiC_Analysis/single_cell_compartment_A_cCREs.bed"
cCRE_file = "/Users/ryanhagan/NoCoSMiCC/HiC_Analysis/single_cell_compartment_B_cCREs.bed"
cCRE_file = "/Users/ryanhagan/NoCoSMiCC/HiC_Analysis/single_cell_compartment_I_cCREs.bed"
cCRE_file = "/Users/ryanhagan/NoCoSMiCC/Synteny/sc_compA_lift_synt.bed"
cCRE_file = "/Users/ryanhagan/NoCoSMiCC/Zoonomia/sc_liftover_compA_synt_G1_cCREs.bed"
cCRE_file = "/Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/LiftOver/sc_compA_lifted.bed"
cCRE_file = "/Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/sc-only-rCARs/sc_bulk_nonOverlap_cCREs.bed"



cCRE_file = "/Users/ryanhagan/NoCoSMiCC/HiC_Analysis/bulk_compartment_A_cCREs.bed"
cCRE_file = "/Users/ryanhagan/NoCoSMiCC/HiC_Analysis/bulk_compartment_B_cCREs.bed"
cCRE_file = "/Users/ryanhagan/NoCoSMiCC/HiC_Analysis/bulk_compartment_I_cCREs.bed"
cCRE_file = "/Users/ryanhagan/NoCoSMiCC/Synteny/bulk_compA_lift_synt.bed"

output_fig = "/Users/ryanhagan/NoCoSMiCC/Zoonomia/Outputs/PhyloP_Aggregation_plot"

cCRE_file = "/Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/sc-only-rCARs/MaxZ/hg38-cCREs-sc-only-lifted.bed"
cCRE_file = "/Users/ryanhagan/NoCoSMiCC/ENCODE_outputs/sc-only-rCARs/MaxZ/hg38-cCREs-sc-only.bed"

df_1 <- read.table(cCRE_file)
head(df_1)
df_2 <- read.table(cCRE_file)
head(df_2)

excluded_ids <- df_1$V4

# Subset df_2 to keep only rows where V4 is NOT in df_1$V4
df_2_subset <- df_2[!df_2$V4 %in% excluded_ids, ]

df <- df_2_subset

cCRE_file = "/Users/ryanhagan/NoCoSMiCC/Zoonomia/sc-G1-lifted-cCREs.bed"
cCRE_file = "/Users/ryanhagan/NoCoSMiCC/Zoonomia/sc-G1-cCREs.bed"
cCRE_file = "/Users/ryanhagan/NoCoSMiCC/Zoonomia/sc-G2-cCREs.bed"
cCRE_file = "/Users/ryanhagan/NoCoSMiCC/Zoonomia/sc-G3-cCREs.bed"

####################
#dat1 = read.table(phyloP100_matrix, row.names = 1)
df = read.table(cCRE_file)
dat2 = read.table(phyloP241_matrix, row.names =1)

#dat1_pls = dat1[df[df$V5=="PLS",]$V4,]
#d1_pls = data.frame(apply(dat1_pls,2,mean))
#d1_pls$ccre = "PLS"
#d1_pls$loci = 1:500
#d1_pls$phyloP = "phyloP100"
#
#dat1_pels = dat1[df[df$V5=="pELS",]$V4,]
#d1_pels = data.frame(apply(dat1_pels,2,mean))
#d1_pels$ccre = "pELS"
#d1_pels$loci = 1:500
#d1_pels$phyloP = "phyloP100"
#
#dat1_dels = dat1[df[df$V5=="dELS",]$V4,]
#d1_dels = data.frame(apply(dat1_dels,2,mean))
#d1_dels$ccre = "dELS"
#d1_dels$loci = 1:500
#d1_dels$phyloP = "phyloP100"
#
#dat1_dnase = dat1[df[df$V5=="DNase-H3K4me3",]$V4,]
#d1_dnase = data.frame(apply(dat1_dnase,2,mean))
#d1_dnase$ccre = "DNase-H3K4me3"
#d1_dnase$loci = 1:500
#d1_dnase$phyloP = "phyloP100"
#
#dat1_ctcf = dat1[df[df$V5=="CTCF-only",]$V4,]
#d1_ctcf = data.frame(apply(dat1_ctcf,2,mean))
#d1_ctcf$ori = "CTCF-only"
#d1_ctcf$ccre = 1:500
#d1_ctcf$phyloP = "phyloP100"
#
dat2_pls = dat2[df[df$V5=="PLS" | df$V5=="PLS,CTCF-bound",]$V4,]
d2_pls = data.frame(apply(dat2_pls,2,mean))
d2_pls$cCRE = "PLS"
d2_pls$loci = 1:500
d2_pls$phyloP = "phyloP241"
#
dat2_pels = dat2[df[df$V5=="pELS" | df$V5=="pELS,CTCF-bound",]$V4,]
d2_pels = data.frame(apply(dat2_pels,2,mean))
d2_pels$cCRE = "pELS"
d2_pels$loci = 1:500
d2_pels$phyloP = "phyloP241"
#
dat2_dels = dat2[df[df$V5=="dELS" | df$V5=="dELS,CTCF-bound",]$V4,]
d2_dels = data.frame(apply(dat2_dels,2,mean))
d2_dels$cCRE = "dELS"
d2_dels$loci = 1:500
d2_dels$phyloP = "phyloP241"
#
dat2_dnase = dat2[df[df$V5=="DNase-H3K4me3" | df$V5=="DNase-H3K4me3,CTCF-bound",]$V4,]
d2_dnase = data.frame(apply(dat2_dnase,2,mean))
d2_dnase$cCRE = "DNase-H3K4me3"
d2_dnase$loci = 1:500
d2_dnase$phyloP = "phyloP241"
#
dat2_ctcf = dat2[df[df$V5=="CTCF-only,CTCF-bound",]$V4,]
d2_ctcf = data.frame(apply(dat2_ctcf,2,mean))
d2_ctcf$ori = "CTCF-only"
d2_ctcf$cCRE = 1:500
d2_ctcf$phyloP = "phyloP241"

dat2_ca= dat2[df[df$V5=="CA-only",]$V4,]
dat2_ca = data.frame(apply(dat2_ca,2,mean))
dat2_ca$ori = "CA-only"
dat2_ca$cCRE = 1:500
dat2_ca$phyloP = "phyloP241"

#colnames(d1_pls) = colnames(d1_pels) = colnames(d1_dels) = colnames(d1_dnase) = colnames(d1_ctcf) = c("score", "ccre","loci","phyloP")
colnames(d2_pls) = colnames(d2_pels) = colnames(d2_dels) = colnames(d2_dnase) = colnames(d2_ctcf) = colnames(dat2_ca) = c("score", "cCRE","loci","phyloP")
d = data.frame(rbind(
  #d1_pls,d1_pels,d1_dels,d1_dnase,d1_ctcf,
                     d2_pls,d2_pels,d2_dels,d2_dnase,d2_ctcf, dat2_ca))
d$cCRE = factor(d$cCRE, levels=c("PLS","pELS","dELS","DNase-H3K4me3","CTCF-only", "CA-only"))
d$phyloP = factor(d$phyloP, levels=c("phyloP241"))

ggplot(d, aes(x=loci, y=score,col=cCRE)) +
  geom_smooth(method="loess",span=0.1, se=FALSE) +
  theme_bw() +
  xlab("Distance from cCRE centre") + ylab("phyloP score") +
  scale_x_continuous(breaks=c(0,250,500),
                     labels=c("-250bp","center","250bp")) +
  geom_vline(xintercept = 250, linetype="dashed", col="gray") +
  scale_color_manual(values=c("#FF0000","#FFA700","#FFCD00","#FFAAAA","#00b0f0", "#bd7ebe")) +
  geom_hline(yintercept = 0.186608, col="#a1a1a1", linetype="solid") +
  geom_hline(yintercept = 0.0817535, col="#a1a1a1", linetype="dashed") +
  geom_hline(yintercept = mean(d$score), col="black", linetype="dashed") +
  ggtitle("G3-overlapping colon cCREs") + 
 ylim(-0.1, 1.6)+
  theme(plot.title = element_text(hjust = 0.5, face="bold"),
        axis.text.x = element_text(face="bold"),
        axis.text.y = element_text(face="bold"),
        axis.title.x = element_text(face="bold"),
        axis.title.y = element_text(face="bold"))

#gsave(output_fig, width=5)


## Pie charts

# Read the cCRE file
cCRE_file <- "/Users/ryanhagan/NoCoSMiCC/Zoonomia/sc-G2-cCREs.bed"
df <- read.table(cCRE_file, header = FALSE, stringsAsFactors = FALSE)

# Define column names
colnames(df) <- c("Chromosome", "Start", "End", "ID", "Category")

# Rename and categorize as in your script
df$Category <- ifelse(df$Category %in% c("PLS", "PLS,CTCF-bound"), "PLS", 
                      ifelse(df$Category %in% c("pELS", "pELS,CTCF-bound"), "pELS", 
                             ifelse(df$Category %in% c("dELS", "dELS,CTCF-bound"), "dELS", 
                                    ifelse(df$Category %in% c("DNase-H3K4me3", "DNase-H3K4me3,CTCF-bound"), "DNase-H3K4me3", 
                                           ifelse(df$Category %in% c("CTCF-only", "CTCF-only,CTCF-bound"), "CTCF-only", 
                                                  ifelse(df$Category == "CA-only", "CA-only", df$Category))))))

cCRE_counts <- table(df$Category)
cCRE_df <- as.data.frame(cCRE_counts)
colnames(cCRE_df) <- c("Category", "Count")

# Define color scheme
cCRE_colors <- c(
  "PLS" = "#FF0000",         # Red
  "pELS" = "#FFA500",        # Orange
  "dELS" = "#FFCD00",        # Royal Blue
  "DNase-H3K4me3" = "#FFAAAA", # Green
  "CTCF-only" = "#00b0f0",   # Purple
  "CA-only" = "#bd7ebe"      # Brown
)

# Generate the pie chart
cCRE_df$Percent <- cCRE_df$Count / sum(cCRE_df$Count) * 100
library(ggrepel)
ggplot(cCRE_df, aes(x = "", y = Count, fill = Category)) +
  geom_bar(width = 1, stat = "identity", color = "white") +
  coord_polar(theta = "y") +
  scale_fill_manual(values = cCRE_colors) +
  theme_void() +
  geom_text_repel(aes(label = paste0(Category, "\n", Count, " (", round(Percent, 1), "%)")), 
                  position = position_stack(vjust = 0.5), 
                  color = "black", 
                  box.padding = 0.5, 
                  point.padding = 0.5) +
  theme(legend.position = "none") 


# Assuming cCRE_df has 'Percent' and 'Count' already calculated
cCRE_df$Percent <- cCRE_df$Count / sum(cCRE_df$Count) * 100

# Create the pie chart with labels inside each wedge
ggplot(cCRE_df, aes(x = "", y = Count, fill = Category)) +
  geom_bar(width = 1, stat = "identity", color = "white") + 
  coord_polar(theta = "y") +  # Keeps the pie chart
  scale_fill_manual(values = cCRE_colors) + 
  theme_void() + 
  # Add the labels showing the Category and percentage of the total
  geom_text(aes(label = paste0(Category, "\n", round(Percent, 1), "%")),
            position = position_stack(vjust = 0.5), 
            color = "black") +  # Label text color
  theme(legend.position = "none")  # Hide the legend

# This script is for making aggregation plot for given cCRE phyloP matrix - adapted from:

# -- Kaili Fan
# fankaili.bio@gmail.com
# Weng Lab, UMass Chan Med
# Jun, 2022

# Load ggplot
library(ggplot2)

phyloP241_matrix = "/project/home/p201120/ryan/cCRE_pipeline/outputs/phyloP241_matrix.mtx"
cCRE_file="/project/home/p201120/ryan/cCRE_pipeline/ENCODE_outputs/colon-epithelial-cCREs-Z1.28.bed"
lifted_cCRE_file = "/project/home/p201120/ryan/cCRE_pipeline/outputs/LiftOver/colon-epithelial-lifted-cCREs.bed"
non_lifted_cCRE_file="/project/home/p201120/ryan/cCRE_pipeline/outputs/LiftOver/colon-epithelial-nonlifted-cCREs.bed"

# Process the data
df = read.table(cCRE_file)
dat2 = read.table(phyloP241_matrix, row.names =1)
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
dat2_dnase = dat2[df[df$V5=="CA-H3K4me3",]$V4,]
d2_dnase = data.frame(apply(dat2_dnase,2,mean))
d2_dnase$cCRE = "CA-H3K4me3"
d2_dnase$loci = 1:500
d2_dnase$phyloP = "phyloP241"
# 
dat2_ctcf  <- dat2[df[df$V5 == "CA-CTCF", ]$V4, ]
d2_ctcf <- data.frame(apply(dat2_ctcf, 2, mean))
d2_ctcf$cCRE <- "CA-CTCF"
d2_ctcf$loci <- 1:500
d2_ctcf$phyloP <- "phyloP241"
# 
dat2_ca_sub <- dat2[df[df$V5 == "CA", ]$V4, ]
d2_ca <- data.frame(apply(dat2_ca_sub, 2, mean))
d2_ca$cCRE <- "CA-only"
d2_ca$loci <- 1:500
d2_ca$phyloP <- "phyloP241"

#colnames(d1_pls) = colnames(d1_pels) = colnames(d1_dels) = colnames(d1_dnase) = colnames(d1_ctcf) = c("score", "ccre","loci","phyloP")
colnames(d2_pls) = colnames(d2_pels) = colnames(d2_dels) = colnames(d2_dnase) = colnames(d2_ctcf) = colnames(d2_ca) = c("score", "cCRE","loci","phyloP")
d = data.frame(rbind(
  #d1_pls,d1_pels,d1_dels,d1_dnase,d1_ctcf,
                     d2_pls,d2_pels,d2_dels,d2_dnase,d2_ctcf, d2_ca))
d$cCRE = factor(d$cCRE, levels=c("PLS","pELS","dELS","CA-H3K4me3","CA-CTCF", "CA-only"))
d$phyloP = factor(d$phyloP, levels=c("phyloP241"))

ggplot(d, aes(x=loci, y=score,col=cCRE)) +
  geom_smooth(method="loess",span=0.1, se=FALSE) +
  theme_bw() +
  xlab("Distance from cCRE centre") + ylab("phyloP score") +
  scale_x_continuous(breaks=c(0,250,500),
                     labels=c("-250bp","center","250bp")) +
  geom_vline(xintercept = 250, linetype="dashed", col="gray") +
  # scale_color_manual(values=c("#FF0000","#FFA700","#FFCD00","#FFAAAA","#00b0f0", "#bd7ebe")) +
  scale_color_manual(values = c(
    "PLS"        = "#E7298A",
    "pELS"       = "#7570B3",
    "dELS"       = "#D95F02",
    "CA-H3K4me3" = "#E6AB02",
    "CA-CTCF"    = "#66A61E",
    "CA-only"    = "#1B9E77"
  )) +
  geom_hline(yintercept = 0.186608, col="#a1a1a1", linetype="solid") +
  geom_hline(yintercept = 0.0817535, col="#a1a1a1", linetype="dashed") +
  geom_hline(yintercept = mean(d$score), col="black", linetype="dashed") +
  ggtitle("Lifted") + 
 ylim(-0.1, 1.6)+
  theme(plot.title = element_text(hjust = 0.5, face="bold"),
        axis.text.x = element_text(face="bold"),
        axis.text.y = element_text(face="bold"),
        axis.title.x = element_text(face="bold"),
        axis.title.y = element_text(face="bold"))
ggsave("phyloP_cCREs.png")
#gsave(output_fig, width=5)


## Pie charts

# Define column names
colnames(df) <- c("Chromosome", "Start", "End", "ID", "Category")

# Rename and categorize as in your script
df$Category <- ifelse(df$Category %in% c("PLS", "PLS,CTCF-bound"), "PLS", 
                      ifelse(df$Category %in% c("pELS", "pELS,CTCF-bound"), "pELS", 
                             ifelse(df$Category %in% c("dELS", "dELS,CTCF-bound"), "dELS", 
                                    ifelse(df$Category %in% c("CA-H3K4me3", "CA-H3K4me3,CTCF-bound"), "CA-H3K4me3", 
                                           ifelse(df$Category %in% c("CA-CTCF", "CA-CTCF,CTCF-bound"), "CA-CTCF", 
                                                  ifelse(df$Category == "CA", "CA-only", df$Category))))))

cCRE_counts <- table(df$Category)
cCRE_df <- as.data.frame(cCRE_counts)
colnames(cCRE_df) <- c("Category", "Count")

# Define color scheme
cCRE_colors <-  c(
    "PLS"        = "#E7298A",
    "pELS"       = "#7570B3",
    "dELS"       = "#D95F02",
    "CA-H3K4me3" = "#E6AB02",
    "CA-CTCF"    = "#66A61E",
    "CA-only"    = "#1B9E77"
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
ggsave("Pie_cCREs.png")

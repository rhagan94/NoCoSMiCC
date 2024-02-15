################################################
################################################

## Script name: _01_Preprocess_cCREs.R

## Purpose: To facilitate the construction of a colonic tissue specific cis regulatory element map.

## Author: Ryan Hagan

## Institute: Royal College of Surgeons in Ireland (RCSI)

## Date created: Feb 2024

## Email: ryanhagan@rcsi.com

################################################
################################################

## Load libraries
library(dplyr)
library(ggplot2)
library(RColorBrewer)

## Create directory and function for saving figures
figsdir <- "./figures/"
if (!dir.exists(figsdir)) {
  dir.create(figsdir)
}

saveFigs <- function(p, width, height) {
  figfile <- file.path(figsdir, sprintf(
    "%s.png",
    gsub("\\.", "_", deparse(substitute(p)))
  ))
  ggsave(figfile, plot = p, width = width, height = height, units = "in", device = "png")
}

## Obtain the filenames of each mouse biosample
mouse_filenames <- list.files(path = "./raw_data/encode_cres/mouse/", pattern = ".bed.gz", full.names = TRUE)
mouse_filenames

## Pre-process the list of cCREs

## Read, subset and filter
lapply(mouse_filenames, function(x) {
  # Read in the raw files
  colon_cCREs <- subset(read.table(x))
  # Subset columns
  colon_cCREs_filt <- subset(colon_cCREs, select = c(V1,V2,V3,V4,V10)) %>%
    filter(V10 != "Unclassified" & V10 != "Low-DNase")
}) -> list_cCREs
list_cCREs

## Set new column names
colnames <- c("Chr", "Start", "End", "cCRE_name", "Classification")
list_cCREs <- lapply(list_cCREs, setNames, colnames)
list_cCREs

## Subset to include only cCRE names
subset_list_cCREs <- lapply(list_cCREs, function(x) x %>% 
                              select(cCRE_name))

## Find unique cCREs
unique_mouse_cCREs <- unique(unlist(subset_list_cCREs))
unique_mouse_cCREs

#unique_mouse_cCREs <- unique(bind_rows(list_cCREs)$cCRE_name)

## Prepare for plotting
merged_mouse_cCREs <- bind_rows(list_cCREs, .id = "Biosample") %>%
  as_tibble() %>%
  group_by(Biosample, Classification) %>%
  count(Classification, sort = TRUE)
merged_mouse_cCREs

options(scipen=10000)
# Define the number of colors you want
nb.cols <- 12
mycolors <- colorRampPalette(brewer.pal(8, "Set1"))(nb.cols)
mouse_cCREs <- ggplot(merged_mouse_cCREs, aes(fill=Classification, y=n, x=Biosample)) + 
  geom_bar(position="stack", stat="identity")+
  scale_fill_manual(values = mycolors) +
  labs(fill = "cCRE Group", y = "Number of cCREs", x = "Biosample")+
  ggtitle("Mouse (Total = 153,949)")+
  theme_bw()+
  theme(legend.position = "right",
        legend.title = element_text(face = "bold"),
        plot.title = element_text(face= "bold", colour= "black" ),
        axis.title.x = element_text(face="bold", colour = "black"),    
        axis.title.y = element_text(face="bold", colour = "black"))
mouse_cCREs
saveFigs(mouse_cCREs, 9, 6)

## Repeat for human biosamples

## Obtain the filenames of each human biosample
human_filenames <- list.files(path = "./raw_data/encode_cres/human/", pattern = ".bed.gz", full.names = TRUE)
human_filenames

## Pre-process the list of cCREs

## Read, subset and filter
lapply(human_filenames, function(x) {
  # Read in the raw files
  colon_cCREs <- subset(read.table(x))
  # Subset columns
  colon_cCREs_filt <- subset(colon_cCREs, select = c(V1,V2,V3,V4,V10)) %>%
    filter(V10 != "Unclassified" & V10 != "Low-DNase")
}) -> list_cCREs
list_cCREs

## Set new column names
colnames <- c("Chr", "Start", "End", "cCRE_name", "Classification")
list_cCREs <- lapply(list_cCREs, setNames, colnames)
list_cCREs

## Identify unique cCREs
unique_human_cCREs <- unique(bind_rows(list_cCREs)$cCRE_name)
length(unique_human_cCREs)

## Prepare for plotting
merged_human_cCREs <- bind_rows(list_cCREs, .id = "Biosample") %>%
  as_tibble() %>%
  group_by(Biosample, Classification) %>%
  count(Classification, sort = TRUE)
merged_human_cCREs

## Set colour params
# Define the number of colors you want
nb.cols <- 17
mycolors <- colorRampPalette(brewer.pal(8, "Set1"))(nb.cols)

merged_human_cCREs$Biosample <- as.factor(merged_human_cCREs$Biosample)
human_cCREs <- ggplot(merged_human_cCREs, aes(fill=Classification, y=n, x=Biosample)) + 
  geom_bar(position="stack", stat="identity")+
  scale_fill_manual(values = mycolors) +
  labs(fill = "cCRE Group", y = "Number of cCREs", x = "Biosample")+
  ggtitle("Human (Total = 436,780)")+
  theme_bw()+
  theme(legend.position = "right",
        legend.title = element_text(face = "bold"),
        plot.title = element_text(face= "bold", colour= "black" ),
        axis.title.x = element_text(face="bold", colour = "black"),    
        axis.title.y = element_text(face="bold", colour = "black"))
human_cCREs
saveFigs(human_cCREs, 12, 6)






table(duplicated(bind_rows(list_cCREs)$cCRE_name))




lapply(human_filenames, function(x) {
  colon_cCREs <- subset(read.table(x))
  return(colon_cCREs)
}) -> list_cCREs

length(list_cCREs[[4]]$V10)





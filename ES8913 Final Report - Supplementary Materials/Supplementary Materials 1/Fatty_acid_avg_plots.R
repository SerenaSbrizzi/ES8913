
## code to create phylogeny figure for fatty acid data

# clear workspace
rm(list = ls())

# set working directory to location on your machine **CHANGE TO YOUR FOLDER**
setwd("~/Desktop/R/Final_phylo_docs")


# install and call libraries that you will need to use
library(ape)
library(tidyverse)
library(geiger)
library(cowplot)
library(RColorBrewer)
library(ggtree)
library(gdata)
library(viridis)



# read in phylogeny in newick format, plot it for a first check, look at tips
phy <- read.tree("final_phy.txt")
plot(phy)
phy$tip.label


# read in tip data (as csv, "taxon" column needs to match the tip names in the phylogeny ***exactly***)
# this column will then be assigned as the rownames of the data to align it with the phylogeny

data <- read.csv("combined_avg_data.csv", header = TRUE, na.strings = "nm")
head(data)
str(data)
rownames(data) <- data$taxon


# check to see if taxon labels are in tips, if FALSE then need to check!
data$taxon %in% phy$tip.label
phy$tip.label %in% data$taxon

# check what's not matching, by making new data frame with ordered columns

check_data <- data.frame(taxon = sort(data$taxon), tips = sort(phy$tip.label))
check_data
# then look at object for typos


# create base plot of the tree in vertical format, can change the tip label size etc.
# will need to adjust the ofset to leave room for the mapping
plot_base <- ggtree(phy) + geom_tiplab(size = 1.3, offset = 0.5)
plot_base


plot_list<- list()

for(i in 2:ncol(data)) {
  curr_df <- data.frame(data[,i])
  rownames(curr_df) <- data$taxon
  
  
  plot_list[[i-1]]<-plot((gheatmap(
    plot_base,
    curr_df,
    offset = -0.6,
    width = 0.03,
    colnames = FALSE,
    font.size = 4,
    family = "test",
    hjust = 0.5) +
    scale_fill_viridis(option = "plasma", begin = 1, end = 0, na.value="white", aes(colour = "Average FAME")) +
    ggtitle(colnames(data)[i])))

      
  Sys.sleep(1)
}

for (i in 2:ncol(data)) {
  file_name = paste ("FA_avg_plot_",(colnames(data)[i]), ".pdf", sep="")
  pdf(file_name, width = 13, height = 10)
  print(plot_list[[i-1]])
  dev.off()
}




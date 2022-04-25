# clear workspace
rm(list = ls())

# set working directory 
setwd("~/Desktop/R/Final_phylo_docs")

# install and call libraries 
library(ape)
library(Rphylopars)
library(data.tree)
library(dplyr)

#open phylogeny 
phy <- read.tree("final_phy.txt")
is.binary(phy)
phy1 <- multi2di(phy, random = TRUE)
is.binary(phy1)
phy1 <- collapse.singles(phy1, root.edge = FALSE)
is.binary(phy1)
phy1 <- multi2di(phy)
phy1$branch.length = NULL
phy1 <- ape::compute.brlen(phy1, method = "Grafen")
summary(phy1$edge.length)
plot(phy1)



#open dataset 
FA_data <- read.csv("combined_avg_data_signal.csv", header = TRUE, na.strings = "nm")

head(FA_data)
str(FA_data)
rownames(FA_data) <- FA_data$species


FA_data$species %in% phy1$tip.label
phy1$tip.label %in% FA_data$species


chi_list<- list()
df_list<-list()
p_list<-list()
final_list<-list()
star_list <- list()
lambda_list<- list()

for(i in 2:ncol(FA_data)) {
  curr_df <- data.frame(FA_data[,c(1,i)])
  rownames(curr_df) <- FA_data$species
  set.seed(i)
  
  p_lambda <- phylopars(trait_data = curr_df, tree = phy1, model = "lambda") ##calculating pagel's lambda for observed data, measuring the degree to which taxa drives trait distribution 

  
  p_star <- phylopars(trait_data = curr_df ,tree = phy1, model = "star") #fit data to star model 
  
  chi_square <- as.double(2*(logLik(p_lambda)-logLik(p_star))) ##comparing log-likelihood: comparing thefir between the lambda and star models to the bserveddata 
  
  degrees_freedom <- (p_lambda$npars - p_star$npars)
  
  p_value <- pchisq(q = chi_square,df = degrees_freedom, lower.tail = FALSE)

  chi_list[[i-1]]<-chi_square
  df_list[[i-1]]<- degrees_freedom
  p_list[[i-1]]<- p_value
  star_list[[i-1]] <- p_star
lambda_list[[i-1]] <-p_lambda

 final_list[[i-1]]<- list(chi=chi_square, df=degrees_freedom, p=p_value, star = p_star, lambda= p_lambda)
 print(final_list[[i-1]])
 
}

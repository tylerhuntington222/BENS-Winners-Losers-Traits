
#########################
# Set directory         
#########################

# set working directory
setwd("/Users/gcclab2/Dropbox/Interface/Dados/Non-Geospatial/Winners_Losers_Traits/Analysis/Data")

# Specify directory for results and figures

results = "~/Dropbox/Interface/Dados/Non-Geospatial/Winners_Losers_Traits/Analysis/Data/results" 
figures = "~/Dropbox/Interface/Dados/Non-Geospatial/Winners_Losers_Trait/Analysis/Data/figures" 

#########################
# Load libraries   
#########################

library(xlsx)

#################
# Set data
#################

# read biodiversity data
beetle_bd.df = read.xlsx("Biodiversity_27_09_16.xlsx", sheetIndex = 2, header = TRUE)

# read bird biodiversity data
#bird_bd.df = read.xlsx("Bird_species_abundance.xlsx", sheetIndex = 1, header = TRUE)

# read beetle trait data
#beetle_trait.df = read.xlsx("Species_traits_23_09_16.xlsx", sheetIndex = 1, header = TRUE)

# read bird trait data
#bird_trait.df = read.xlsx("Bird_species_traits.xlsx", sheetIndex = 1, header = TRUE)

head(beetle_bd.df)

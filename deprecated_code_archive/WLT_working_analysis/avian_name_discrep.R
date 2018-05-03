
# set working directory
setwd("/Users/nicholslab/Dropbox/Interface/Dados/Non-Geospatial/Winners_Losers_Traits/Analysis/Data")

# load libraries

library(xlsx)
library(reshape)
library(plyr)
library(reshape2)
library(dplyr)
library (rgdal)
library (rgeos)
library (spdep)
library (plyr)
library(ggplot2)


# load avian species names in original bird species traits xls
o.tr.df = read.xlsx ( "Bird_species_traits_21_11_16.xlsx" , sheetIndex = 1 , header = TRUE )

# subset column of species names and save as o.df
o.df <- data.frame ( o.tr.df [,2] ) 
colnames(o.df) <- c('species')
head(o.df)

# load avian species names in most recent (corrected) bird species traits xls
c.tr.df = read.xlsx ( "Bird_species_traits_orig.xlsx" , sheetIndex = 1 , header = TRUE )

# create df with species names and save as c.df
c.df <- data.frame (c.tr.df [,2])
colnames(c.df) <- c('species')
head (c.df)
# anti join o.df and c.df to determine naming discrepencies

discrep <- anti_join(c.df , o.df , by = 'species')


discrep

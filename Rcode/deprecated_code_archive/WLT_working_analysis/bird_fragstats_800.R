#########################
# Set directory         
#########################

setwd("/Users/nicholslab/Dropbox/Interface/Dados/Non-Geospatial/Winners_Losers_Traits/Analysis/Data/a.pts.b800.LC.shapes")

#########################
#Load libraries   
#########################

library (rgdal)
library (rgeos)
library (spdep)
library (plyr)
library(ggplot2)


#########################
# Notes 
#########################
# area here will be in square meters

TAh = TAm / 10000

#################
# Set data
#################

# read CSV with bird sample point codes
a.pts1.v <- read.csv("/Users/nicholslab/Dropbox/Interface/Dados/Non-Geospatial/Winners_Losers_Traits/Analysis/Data/Bird_Sampling_Points/a.pts.csv", header = F)

a.pts.v <- as.factor(a.pts1.v[2:157,2])
a.pts.v

# set up empty df to write results into
df = NULL


for (point in a.pts.v){
  
  # set up file locations and names to read from
  
  #path = paste("/Users/nicholslab/Dropbox/Interface/Dados/Non-Geospatial/Winners_Losers_Traits/Analysis/Data/a.pts.b800.dissolves/", point,".b.800.dissolve/" , point, ".b.800.dissolve.shp", point, sep = "")

  #file = paste(point, ".b.800.dissolve", sep = "")
  
  # read in file
  
  #ptLS <- subset ( readOGR(dsn = path.expand(path), layer = file, [2,] ))
  
  file <- paste("/Users/nicholslab/Dropbox/Interface/Dados/Non-Geospatial/Winners_Losers_Traits/Analysis/Data/a.pts.b800.dissolves/", point,".b800.dissolve/" , point, ".b800.dissolve.shp", sep = "")
  
  layer <- paste (point, ".b800.dissolve", sep = "")
  
  buff <- readOGR (file , layer)
  
  #Calculate total edge for each LC type

  edge <- gBoundary(buff, byid = T) # calculates length of boundary edge
  perims <- data.frame(gLength(edge, byid = T)) # returns length of geometry in current projection
  
  # subset TE for forest
  
  TE.forest <- perims[2,]
  
  #Calculate percent area of each LC type 
  
  perArea <- data.frame(buff$Area/2010537)
  
  perArea.FC <- perArea [2,]
  
  # bind together columns for output CSV
  
  buff.stats <- as.data.frame(cbind (point, TE.forest, perArea.FC))
  
  df = rbind(buff.stats, df)
  
  }

# define ds as fragstats.df

fragstats.df <- data.frame(df)


# Write fragstats.df to CSV

write.csv (fragstats.df, file = "a.pts.b800.fragstats.csv")

buff$Area
  
sum(122368.68, 325645.03, 320297.96,  60244.39, 286938.65, 895041.91)
perArea
sum(perArea)


#########################

# Goal: recalculate Fragstats for 12 Cantareira landscapes

#########################

#########################
# Set directory         
#########################

setwd("~/Dropbox/Interface/Dados/Geospatial/Private_LEPAC/Final products/Pitfalls/Buffers_forest_pitfalls/400m") #kelley
# setwd("~/Dropbox/Grants/Active/Interface/Dados/Geospatial/Private_LEPAC/Final products/Pitfalls/Buffers_forest_pitfalls/400m") #liz

#######################
# Define functions
####################

# Multiple plot function
#
# ggplot objects can be passed in ..., or to plotlist (as a list of ggplot objects)
# - cols:   Number of columns in layout
# - layout: A matrix specifying the layout. If present, 'cols' is ignored.
#
# If the layout is something like matrix(c(1,2,3,3), nrow=2, byrow=TRUE),
# then plot 1 will go in the upper left, 2 will go in the upper right, and
# 3 will go all the way across the bottom.
#
multiplot <- function(..., plotlist=NULL, file, cols=1, layout=NULL) {
  library(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # If layout is NULL, then use 'cols' to determine layout
  if (is.null(layout)) {
    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
                     ncol = cols, nrow = ceiling(numPlots/cols))
  }
  
  if (numPlots==1) {
    print(plots[[1]])
    
  } else {
    # Set up the page
    grid.newpage()
    pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))
    
    # Make each plot, in the correct location
    for (i in 1:numPlots) {
      # Get the i,j matrix positions of the regions that contain this subplot
      matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))
      
      print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }
}

#########################
#Load libraries   
#########################

rm(list = ls())
library (rgdal)
library (rgeos)
library (spdep)
library (plyr)
library(ggplot2)

#################
# Notes
#################
# area here will be in square meters

TAh = TAm / 10000

#################
# Set data
#################

# set up empty df to write results into
df = NULL

landscapes <- c("148", "215", "263", "266", "282", "291", "317", "329", "333", "335", "359", "399")

for (LS in landscapes){
  
  # set up file locations and names to read from
  
  path = paste("~/Dropbox/Interface/Dados/Geospatial/Private_LEPAC/Final products/Pitfalls/Buffers_forest_pitfalls/400m/b_400_",LS,"/b_400_",LS,"_final", sep ="")
  # path = paste("~/Dropbox/Grants/Active/Interface/Dados/Geospatial/Private_LEPAC/Final products/Pitfalls/Buffers_forest_pitfalls/400m/b_400_", LS,"/b_400_",LS,"_final", sep ="") #liz
  file = paste("b_400_", LS, "_dissolve", sep="")
  
  # read in file
  
  pLS <- subset(readOGR(dsn = path.expand(path), layer = file), Classe == "forest")

#############
# total edge
##############
  
  edge <- gBoundary(pLS, byid = T) # calculates length of boundary edge
  perims <- gLength(edge, byid = T) # returns length of geometry in current projection
  pLS$perims <- perims

######################
# mean patch isolation
####################

  # Calculated as mean nearest neighbor distance per trap
  # Calculate distance between a patch and all other patches within a buffer

  # Calculate smallest edge distance between patches in pLS
  Fdist <- list()
  for(i in 1:dim(pLS)[1]) { 
    pDist <- vector()
    for(j in 1:dim(pLS)[1]) { 
      if(pLS$Code[i] == pLS$Code[j]) { # do only for patches in same trap
        if(pLS$perims[i] != pLS$perims[j]){ # don't compare a patch to itself
        pDist <- append(pDist, gDistance(pLS[i,],pLS[j,])) 
      }}}
    Fdist[[i]] <- pDist
  } 

  # reads out file with a list of neighbor distances for each patch
  # zeros are from polygons sitting up against main polygon
  # polygons from buffers with no other polygons have logical(0)

  # Extract nearest neighbor distance
  PolyDist <- unlist(lapply(Fdist, FUN=function(x) min(x)[1])) 
  # returned infinity for logical 0s, change this to 0
  PolyDist[is.infinite(PolyDist)] <- 0 

  pLS$PolyDist <- PolyDist
  
####################
# Number of patches
###################

  # create a column with a "1" for each patch
  # when we ddply for point, sum up these 1s, get # of patches

  pLS$num_p <- 1

##################
# ddply per trap
##################
  
  # coerce into df
  pLS.df <- as.data.frame(pLS)
  # summarize fragstats per trap
  pLS.trap <- ddply(pLS.df, .(Landscape, Point, Trap_num, Code_min, Code), summarize,
                   TE = sum(perims),
                   MPI = mean(PolyDist),
                   NP = sum(num_p))

  df = rbind(df, pLS.trap)
}

################
# save results as .csv
################

data = "~/Dropbox/Interface/Dados/Non-Geospatial/Biodiversity/Analysis/Data" #Kelley
filename <- paste(data, "Fragstats.csv", sep = "/")
write.csv(df, file = filename, row.names = FALSE)

############################
# make PDF with histograms of each metric per LS
##########################

figures = "~/Dropbox/Interface/Dados/Non-Geospatial/Biodiversity/Analysis/Data/figures" #Kelley

filename <- paste(figures, "Fragstats.pdf", sep = "/")
pdf(file = filename, width = 20, height = 40)

par(mfrow=c(3,1))
# b, l, t, r
par(oma=c(3,3,3,3))
par(mar=c(17,20,5,20))

p1 <- ggplot(df, aes(x = TE)) + 
      geom_histogram( binwidth = 150) +
      facet_grid(Landscape ~ .) +
      xlab("Total Edge") +
      theme(text = element_text(size = 30))

p2 <- ggplot(df, aes(x = MPI)) + 
      geom_histogram( binwidth = 10) +
      facet_grid(Landscape ~ .) +
      xlab("Mean Patch Isolation") +
      theme(text = element_text(size = 30))

p3 <- ggplot(df, aes(x = NP)) + 
      geom_histogram( binwidth = 1) +
      facet_grid(Landscape ~ .) +
      xlab("Number of Patches") +
      scale_x_continuous( breaks = c(0,1,2,3,4,5,6,7,8,9)) +
      theme(text = element_text(size = 30))

multiplot(p1, p2, p3, cols = 3)

dev.off()


############################################## CREATE MERGED LS BASEMAP ##############################################

# read in centroids for landscapes

ls.centroids1 = readOGR("/Users/nicholslab/Dropbox/Interface/Dados/Non-Geospatial/Winners_Losers_Traits/Analysis/Data/geospatial/LS_centroids.shp", "LS_centroids")
ls.centroids <- spTransform (ls.centroids1, CRS("+init=epsg:32723"))

# create empty 5km buffer polygons around LS centroids
ls.centroids_buffer5k <- gBuffer(ls.centroids, width=5000, byid=TRUE, quadsegs=100)

# write shapefile of dissolved 5km empty buffer polygon
writeOGR(obj=ls.centroids_buffer5k , dsn="ls.centroids_buffer5k f", layer="ls.centroids_buffer5k", driver="ESRI Shapefile") 

plot(ls.centroids_buffer5k)


coordinates((ls148.buff.3k))
# read 5km buffer land cover polygon 




# dissolve 

empty.3k.buff <- unionSpatialPolygons(ls148.buff.3k, ls148.buff.3k$Landscape, avoidGEOS = T)

empty.3k.buff <- unionSpatialPolygons(ls148.buff.3k, c("Landscape", "area", "ID", "Classe"))

help(unionSpatialPolygons)

assign (paste("a.", i , ".b800.LC" , sep = "") , subset (a.all.b800.LC.insct, a.all.b800.LC.insct$ponto == i))

print(paste("a.", i ,".b800.LC" , sep = ""))

}

# write shapefile for individual buffer LC intersect product

# for (i in a.pt.LS.codes.v) {
writeOGR(paste("a.", i ,".b800.LC" , sep = ""), dsn = paste("a.", i ,".b.800.LC", sep = ""), layer = paste("a." , i ,".b_800.LC", sep = ""), driver = "ESRI Shapefile")

}

writeOGR( obj = a.P359P02P.OK.b800.LC, dsn = "a.P359P02P.OK.b800.LC", layer = "a.P359P02P.OK.b800.LC", driver = "ESRI Shapefile")

summary(a.P359P02P.OK.b800.LC)

print(a.pt.LS.codes.v)
head(a.all.b800.LC.insct)

############################################## LOOP FOR SUBSETTING AVIAN POINT BUFFER INTERSECT PRODUCTS (V. 2) ##############################################

for (ponto in a.all.b800.LC.insct) {
  
  # subset individual point buffer from a.buff
  
  paste("a.", ponto ,".b800.LC") <- subset (a.all.b800.LC.insct, a.all.b800.LC.insct$ponto == ponto)
  
  # write shapefile for individual buffer LC intersect product
  
  writeOGR(obj = paste("a.", ponto ,".b800.LC", dsn = paste("~/Dropbox/Interface/Dados/Geospatial/Non-Geospatial/Winners_Losers_Traits/Analysis/Data/Bird_Sampling_Points/final_buffer_insct_products/", "a.", ponto ,"_b_800_LC", sep = ""), 
                       layer = paste(ponto ,"_b_800_LC", sep =""), 
                       driver = "ESRI Shapefile"))}

head(a.pt.LS.codes.df)
head(a.all.b800.LC.insct)


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
  
  
  library(PBSmapping)
  
  p148.1 = readOGR(dsn = "/Users/nicholslab/Dropbox/Interface/Dados/Non-Geospatial/Winners_Losers_Traits/Analysis/Data/a.pts.b800.LC.shapes/a.ENT.P148P10M.b.800.LC/a.ENT.P148P10M.b_800.LC.shp", layer = "a.ENT.P148P10M.b_800.LC")
  summary (p148.1)
  Class <- as.character(p148.1$Class)
  Class
  
  
  p148.1.PS = SpatialPolygons2PolySet(p148.1)
  
  p148.1.dslv <- unionSpatialPolygons(p148.1 , Class)
  
  gArea(p148.1.dslv, T)
  
  SpatialPolygonsDataFrame(p148.1.dslv)
  # make sure row names match
  row.names(p148.1.dsv.df) <- as.character(1:length(p148.1.dslv))
  
  # Extract the data you want (the larger geography)
  p148.1.dsv.df1 <- unique(p148.1$Class)
  p148.1.dsv.df <- as.data.frame(p148.1.dsv.df1)
  
  colnames(p148.1.dsv.df) <- "Class"  # your data will probably have more than 1 row!
  
  # And add the data back in
  p148.1.dsv.SPDF <- SpatialPolygonsDataFrame(p148.1.dslv , p148.1.dsv.df)
  
  
  LC.areas.df <- data.frame(gArea(p148.1.dsv.SPDF, TRUE))
  
  
  rownames(LC.areas.df) <- c("eucalyptus" , "forest" , "human_sett"  ,"pasture" , "water")
  
  gArea(p148.1, T)
  
  summary (p148.1.dslv, TRUE)
  
  summary (p148.1, T)
  
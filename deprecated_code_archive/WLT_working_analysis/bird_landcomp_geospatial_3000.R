# Geospatial Analyses for Bird Samping Points

# Goals

# 1) Load shapefile with all bird sampling points

# 2) Create empty buffers at 3000m radius around each point

# 3) Intersect buffers with landscape rasters

# 4) Dissolve polygons by like land cover type


# set working directory

setwd("/Users/Tyler/Dropbox/Interface/Dados/Non-Geospatial/Winners_Losers_Traits/Analysis/Data/Bird_Sampling_Points")

# set results directory
results = "~/Dropbox/Interface/Dados/Non-Geospatial/Winners_Losers_Traits/Analysis/Data/results" 

# load packages

library (PBSmapping)
library (sp) 
library (rgdal) 
library (raster)
library (rgeos) 
library (maptools)
library (xlsx)
library (reshape)
library (plyr)
library (reshape2)
library (dplyr)
library (raster)
library (GISTools)
library (gpclib)
library (rgdal)
library (spdep)
library (ggplot2)



############################################## LOAD DATA ##############################################

# load bird sampling points shapefile
a.pts.geo1 <- readOGR("/Users/Tyler/Dropbox/Interface/Dados/Non-Geospatial/Winners_Losers_Traits/Analysis/Data/Bird_Sampling_Points/Point_counts_Cantareira.shp", "Point_counts_Cantareira")

# load basemap of collated landscapes
all.ls.merge1 <- readOGR("/Users/Tyler/Dropbox/Interface/Dados/Geospatial/Private_LEPAC/Final products/Mapeamento/merged_land_use/all_landscapes_merge.shp" , "all_landscapes_merge")

# clean up topology errors
all.ls.base <- gBuffer(all.ls.merge1, byid=TRUE, width=0.0)

# check to make sure no topology errors 
gIsValid(all.ls.merge)



############################################## DATA MANIP ##############################################

# Retransform data to WGS 84/UTM 23S 
a.pts<- spTransform (a.pts.geo1, CRS("+init=epsg:32723"))

# create buffer polygons of 3000m radius around each point
a.buff <- gBuffer(a.pts, width=3000, byid=TRUE, quadsegs=100)

# write shapefile a.buff
writeOGR(obj=a.buff, dsn="a.buff", layer="a.buff", driver="ESRI Shapefile") 

# double check areas of buffer polygons 
round(sqrt(gArea(a.buff, byid=TRUE)/pi))

plot(a.buff)

# clean up topology errors of collated landscapes basemap

all.ls.base <- gBuffer(all.ls.merge1, byid=TRUE, width=0.0)

# check to make sure no topology errors 
gIsValid(all.ls.base)

# Retransform data to WGS 84/UTM 23S 
ls.base <- spTransform (all.ls.base, CRS("+init=epsg:32723"))



# subset individual empty point buffer from a.buff

#a.148.buff.emty = a.buff.emty[a.buff.emty$ponto=="P148P17P OK",]

#a.148.buff.emty


# intersect individual empty buffer polygons with base map

#a.148.buff_intersect <- raster::intersect(a.148.buff , ls.base)

#plot (a.148.buff_intersect)





# intersect buffer polygons with base map
a.allpts_insct<- raster::intersect(a.buff , ls.base)

writeOGR(obj=a.allpts_insct, dsn="a.allpts_insct", layer="a.allpts_insct", driver="ESRI Shapefile")


# subset individual point buffer from a.buff

a.148.buff_insct = a.allpts_insct[a.allpts_insct$ponto=="P148P17P OK",]



writeOGR(obj=a.148.buff_insct , dsn="a.148.buff_insct", layer="a.148.buff_insct", driver="ESRI Shapefile")

plot(a.148.buff_insct)

############################################## LOOP FOR SUBSETTING AVIAN POINT BUFFER INTERSECT PRODUCTS ##############################################

# load merged avian sampling point 3000m buffers intersected universal LC basemap

a.all.b3000.LC.insct <- readOGR("/Users/Tyler/Dropbox/Interface/Dados/Non-Geospatial/Winners_Losers_Traits/Analysis/Data/Bird_Sampling_Points/a.allpts_insct" , "a.allpts_insct")

a.all.b3000.LC.insct$ponto <- gsub (" ", ".",a.all.b3000.LC.insct$ponto, fixed = TRUE) 

a.all.b3000.LC.insct$ponto

# Create vector containing point LS codes

a.pt.LS.codes <- (a.all.b3000.LC.insct$ponto)

a.pt.LS.codes.v1 <- as.factor(a.pt.LS.codes)
a.pt.LS.codes.v <- unique(a.pt.LS.codes.v1)



# Enter loop

for (i in a.pt.LS.codes.v) {
  
  # subset individual point buffer from a.buff
  
  writeOGR(subset (a.all.b3000.LC.insct, a.all.b3000.LC.insct$ponto == i), dsn = paste("a.", i ,".b.3000.LC", sep = ""), layer = paste("a." , i ,".b.3000.LC", sep = ""), driver = "ESRI Shapefile")
  
}

############################################## LOOP FOR DISSOVLING IND AVIAN POINT BUFFER INTERSECT PRODUCTS ##############################################

for (pt in a.pt.LS.codes.v) {
  
  point = readOGR(dsn = paste("/Users/Tyler/Dropbox/Interface/Dados/Non-Geospatial/Winners_Losers_Traits/Analysis/Data/a.pts.b3000.LC.shapes/a.", pt, ".b.3000.LC/a." , pt ,".b.3000.LC.shp" , sep = ""), layer = paste ("a." , pt , ".b.3000.LC" , sep = ""))
  
  Class <- as.character(point$Class)
 
  ps <- SpatialPolygons2PolySet(point)
  
  dsv <- unionSpatialPolygons(point , Class)
  
  LC.areas.df <- data.frame(gArea(dsv, TRUE))
 
  colnames( LC.areas.df ) <- c("Area")
  
  unq_Class <-   as.character(unique(Class))
  
  rownames (LC.areas.df) <- unq_Class
  
  spdf <- SpatialPolygonsDataFrame(dsv,  LC.areas.df)
  
  writeOGR(spdf, dsn = paste(pt, ".b3000.dissolve", sep = ""), layer = paste(pt, ".b3000.dissolve", sep = ""), driver = "ESRI Shapefile")
}



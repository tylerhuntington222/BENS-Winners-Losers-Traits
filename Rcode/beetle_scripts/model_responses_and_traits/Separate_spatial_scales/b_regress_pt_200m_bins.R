#####################################################################################################################################
################################## REGRESSION ANALYSIS OF BEETLES USING POINT LEVEL ABUNDANCES BINNED ############################### 
#####################################################################################################################################
# USES FC FROM 200m BUFFERS AROUND POINTS

# set working directory
setwd("/Users/Tyler/Dropbox/Interface/Dados/Non-Geospatial/Winners_Losers_Traits/Analysis/Data")

# load beetle dataframes
b.df <- readRDS ( "b.df.rds" )

head(b.df)

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
library (caper)
library (ape)
library (MuMIn)


################################## MANIPULATE DATA FOR THIS PARTICULAR APPROACH ################################## 

# create df with unique entry for each beetle species and associated attributes
b.sp.df <- ddply ( b.df, .(Species , Diel , Nest , Biogeo , BM , Hab , Diet, N.tot), summarize, mean(N.tot))
b.sp.df <- b.sp.df[,1:8]


### change perFC column in b.df to contain 200m buffer FC levels

# import 200m buffer FC data
fc.pt.200m.df <- read.csv("Landscape_composition_pitfalls.csv", header = T)

# subset for 200m buffer
fc.pt.200m.df <- subset (fc.pt.200m.df, fc.pt.200m.df$Buffer_m == 200)

# subset for forest cover
fc.pt.200m.df <- subset (fc.pt.200m.df, fc.pt.200m.df$Class == 'forest')

# collapse traps to point level, averaging trap buffer FC values
fc.pt.200m.df <- ddply (fc.pt.200m.df, .(Landscape, Point, Code_min), 
                        summarize, perFC = mean(Per_area))

colnames(fc.pt.200m.df) <- c( "Landscape", "Point" , "TrapID" , "pt.200m.perFC" )

# subset for only TrapID and 200m FC cols
fc.pt.200m.df <- subset (fc.pt.200m.df, select = c( "TrapID", "pt.200m.perFC" ))

# merge 200m pt buffer FC data with b.df
b.df <- merge(b.df, fc.pt.200m.df, by = "TrapID")

# replace perFC column contents with pt.200m.perFC data
b.df$perFC <- b.df$pt.200m.perFC

# remove extranneous pt.200m.perFC column
b.df <- b.df[,1:15]

### Binning abundances 
head(b.df)

length(unique(b.df$TrapID))
length(unique(b.df$perFC))

perFC.all.levels <- unique(b.df$perFC)

perFC.breaks <- c(0,10,20,30,40,50,60,70,80,90,100)

perFC.binCodes <- .bincode(perFC.all.levels, perFC.breaks, right = F)

# bind perFC values to bin assignments in dataframe
perFC.bin.key <- data.frame(perFC.all.levels, perFC.binCodes) 
colnames(perFC.bin.key) <- c("perFC.bin" , "FC.binCode")

# merge FC bin key with b.df
b.df <- merge (b.df, perFC.bin.key)

# aggregate abbundances at FC bin level
b.df <- ddply (b.df , .(Species, Diel, Nest, Biogeo, Diet, BM, N.tot, FC.binCode), summarise,
               perFC.bin = mean(perFC), N.bin = sum(N))



print( b.df$N.bin)

################################################ Plot N ~ FC for all each beetle species

# make plot matrix of individual beetle species N ~ FC
par(mfrow=c(2,3))
par (pty = "s")
par(mar=c(2, 2,  2,  2))

for (i in b.sp.df$Species) {
  
  # create df with entries as all captures of a particular species
  b.par.df <- subset (b.df , b.df$Species == i & b.df$Presence == 1) # omitting/including zeros doesn't seem to affect model
  
  # plot N ~ FC for particular species
  plot(b.par.df$perFC.bin, b.par.df$N.bin, main = i)
  
  print (i)
  
}

################################################ CALCULATE N ~ FC RESPONSES FOR ALL BEETLE SPECIES

# create blank beetle FC response df
b.respFC.df <- data.frame(
  Species = factor(), 
  Intercept = double(),
  Inter.error = double(),
  Inter.t.val = double(),
  Inter.p.val = double(), 
  Slope = double (),
  Slope.error = double (),
  Slope.t.val = double(),
  Slope.p.val = double()
) 
# set levels of species var
levels(b.respFC.df$Species) <- unique (b.sp.df$Species)

# initialize row 1 as first slot of b.respFC.df to be filled
row = 1

# iterate response calculation over all beetle species and store in b.respFC.df
for (i in b.sp.df$Species) {
  
  # create df with entries as all captures of a particular species
  b.par.df <- subset (b.df , b.df$Species == i ) # & b.df$Presence == 1
  
  # plot N ~ FC for particular species
  b.lm.par <- lm (b.par.df$N.bin ~ b.par.df$perFC.bin)
  
  print (summary (b.lm.par)$coefficients)
  
  b.respFC.df[row, 1] <- (i)
  b.respFC.df[row,2:5] <- summary (b.lm.par)$coefficients[1,]
  b.respFC.df[row,6:9] <- summary (b.lm.par)$coefficients[2,]
  
  row <- row + 1
}


print (b.respFC.df)
summary (b.respFC.df)

hist (b.respFC.df$Slope, breaks = 20)

# species with significant slope coefficients
b.sigResp <- subset (b.respFC.df, b.respFC.df$Slope.p.val < 0.05)

hist (b.sigResp$Slope, breaks = 15)
hist (b.sigResp$Slope.p.val, breaks = 15)


#################################################################### 
# OPTION B: CALCULATING ABUNDANCES AT TRANSECT LEVEL
#################################################################### 

# load beetle and avian dataframes
b.df <- readRDS ( "b.df.rds" )

# calculate abundances at transect level
b.df <- ddply (b.df , .(Species, Landscape, Diel, Nest, Biogeo, Diet, BM, Hab, perFC), 
               summarise, N.ls = sum(N), LS.tot = mean(LS.total))

# create a column for presence absence in a landcape
b.df$LS.presence <- ifelse ( paste (b.df$N.ls) > 0, 1, 0 )

head(b.df)

sort(b.df$N.ls)

################################################ Plot N ~ FC for all each beetle species

# make plot matrix of individual beetle species N ~ FC
par(mfrow=c(2,3))
par (pty = "s")
par(mar=c(2, 2,  2,  2))

for (i in b.sp.df$Species) {
  
  # create df with entries as all captures of a particular species
  b.par.df <- subset (b.df , b.df$Species == i & b.df$LS.presence == 1)
  
  # plot N ~ FC for particular species
  plot(b.par.df$perFC, (b.par.df$N.bin), main = i)
  
  print (i)
  
}

################################################ CALCULATE N ~ FC RESPONSES 

# create blank beetle FC response df
b.respFC.df <- data.frame(
  Species = factor(), 
  Intercept = double(),
  Inter.error = double(),
  Inter.t.val = double(),
  Inter.p.val = double(), 
  Slope = double (),
  Slope.error = double (),
  Slope.t.val = double(),
  Slope.p.val = double()
) 
# set levels of species var
levels(b.respFC.df$Species) <- unique (b.sp.df$Species)

# initialize row 1 as first slot of b.respFC.df to be filled
row = 1

# iterate response calculation over all beetle species and store in b.respFC.df
for (i in b.sp.df$Species) {
  
  # create df with entries as all captures of a particular species
  b.par.df <- subset (b.df , b.df$Species == i )
  
  # plot N ~ FC for particular species
  b.lm.par <- lm (b.par.df$N.bin ~ b.par.df$perFC)
  
  print (summary (b.lm.par)$coefficients)
  
  b.respFC.df[row, 1] <- (i)
  b.respFC.df[row,2:5] <- summary (b.lm.par)$coefficients[1,]
  b.respFC.df[row,6:9] <- summary (b.lm.par)$coefficients[2,]
  
  row <- row + 1
}

print(b.respFC.df)

################################################ MERGE N ~ FC RESPONSES WITH b.sp.df

b.sp.df <- merge(b.respFC.df, b.sp.df)

head(b.sp.df)

# make row names of b.sp.df be species names
rownames(b.sp.df) <- b.sp.df$Species


################################################ MODEL RESPONSES ~ TRAITS | PHYLO

#### Import phylo tree object from Kelley here
b.tre <- read.tree("b.sp.tre")
b.tre$node.label <- NULL

b.comp <- comparative.data(b.tre, b.sp.df, Species, vcv=TRUE)

# build model
mod <- pgls(Slope ~ log(BM) + Diet + Diel + Nest , b.comp)

# + Hab  --- generates error
# + Biogeo ---generates error

print(mod)



# generate all possible models and rank by AIC

all.poss.mod <- dredge(mod, beta = c("none", "sd", "partial.sd"), evaluate = TRUE,
                       rank = "AICc", fixed = NULL)

print (all.poss.mod)


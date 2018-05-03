
#########################
# Set directory         
#########################

# set working directory
setwd("/Users/tyler/Dropbox/Interface/Dados/Non-Geospatial/Winners_Losers_Traits/Analysis/Data")

list.files()
# Specify directory for results and figures

results = "~/Dropbox/Interface/Dados/Non-Geospatial/Winners_Losers_Traits/Analysis/results" 
figures = "~/Dropbox/Interface/Dados/Non-Geospatial/Winners_Losers_Trait/Analysis/figures" 

#########################
# Load libraries   
#########################
library(readxl)
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
library (gdata)

#################
# Read in data
#################

########################## BIODIVERSITY DATA ##########################

# read beetle biodiversity data
b.bdv.df2 = read_excel("beetle_biodiversity_27_09_16.xlsx", sheet = 2, col_names = TRUE)

b.bdv.df1 <- reshape::rename ( b.bdv.df2 , replace = c ("Code"="TrapID") )

########################## TRAIT DATA #################################

# read beetle trait data
b.tr.df1 = read_excel("beetle_species_traits_23_09_16.xlsx", sheet = 1, col_names = TRUE)

# rename Code to code
b.tr.df1 <- reshape::rename (b.tr.df1, replace = c ( "Code" = "code" ) )

########################## LANDSCAPE DATA #############################

# key for TrapID_min codes
TrapID.key.df<- read_excel("Trap_Code_Key.xlsx", sheet = 1, col_names = T)

# # read beetle FC data for all landscapes 
# b.fc.df1 <- read.xls("LS_summary_stats.xls", sheet = 3, header = T)
# 
# # subset landscape and perFC columns
# b.fc.df1 <- subset (b.fc.df1, select = c( "Landscape" , "X3km_FC." )) 
# 
# # rename columns
# b.fc.df <- reshape::rename(b.fc.df1, c(X3km_FC. = "LS_perFC"))
# 
# # read beetle landscape TE data (200m buffer) and rename column headers appropriately
# # b.te.df1 <- read.csv("Fragstats_200.csv", header = T)
# # b.te.df<- reshape::rename(b.te.df1, replace = c("Code"="TrapID", "Code_min"="TrapID_min" ))


# read in forest cover data
a.landcomp.df <- read.csv("Landscape_composition_pitfalls.csv", header = T)

# subset for 200m buffer
fc.pt.200m.df <- subset (a.landcomp.df, a.landcomp.df$Buffer_m == 200)

# subset for forest fragments
fc.pt.200m.df <- subset (fc.pt.200m.df, fc.pt.200m.df$Class == 'forest')

# collapse traps to point level, averaging trap buffer FC values
fc.pt.200m.df <- ddply (fc.pt.200m.df, .(Landscape, Point, Code_min), 
                        summarize, perFC_200 = mean(Per_area))




# subset for 800m buffer
fc.pt.400m.df <- subset (a.landcomp.df, a.landcomp.df$Buffer_m == 400)

# subset for forest fragments
fc.pt.400m.df <- subset (fc.pt.400m.df, fc.pt.400m.df$Class == 'forest')

# collapse traps to point level, averaging trap buffer FC values
fc.pt.400m.df <- ddply (fc.pt.400m.df, .(Landscape, Point, Code_min), 
                        summarize, perFC_400 = mean(Per_area))


# subset for 800m buffer
fc.pt.800m.df <- subset (a.landcomp.df, a.landcomp.df$Buffer_m == 800)

# subset for forest fragments
fc.pt.800m.df <- subset (fc.pt.800m.df, fc.pt.800m.df$Class == 'forest')

# collapse traps to point level, averaging trap buffer FC values
fc.pt.800m.df <- ddply (fc.pt.800m.df, .(Landscape, Point, Code_min), 
                        summarize, perFC_800 = mean(Per_area))



# subset for 1000m buffer
fc.pt.1000m.df <- subset (a.landcomp.df, a.landcomp.df$Buffer_m == 1000)

# subset for forest fragments
fc.pt.1000m.df <- subset (fc.pt.1000m.df, fc.pt.1000m.df$Class == 'forest')

# collapse traps to point level, averaging trap buffer FC values
fc.pt.1000m.df <- ddply (fc.pt.1000m.df, .(Landscape, Point, Code_min), 
                         summarize, perFC_1000 = mean(Per_area))


# subset for 3000m buffer
fc.pt.3000m.df <- subset (a.landcomp.df, a.landcomp.df$Buffer_m == 3000)

# subset for forest fragments
fc.pt.3000m.df <- subset (fc.pt.3000m.df, fc.pt.3000m.df$Class == 'forest')

# collapse traps to point level, averaging trap buffer FC values
fc.pt.3000m.df <- ddply (fc.pt.3000m.df, .(Landscape, Point, Code_min), 
                         summarize, perFC_3000 = mean(Per_area))


# merge 200, 400, 800, 1000, and 3000 perFC dataframes
b.lc.df_merge.df <- merge(fc.pt.200m.df, fc.pt.400m.df)
b.lc.df_merge.df <- merge(b.lc.df_merge.df, fc.pt.800m.df)
b.lc.df_merge.df <- merge(b.lc.df_merge.df, fc.pt.1000m.df)
b.fc.df <- merge(b.lc.df_merge.df, fc.pt.3000m.df)

# rename Code_min col name to TrapID_min
names(b.fc.df)[names(b.fc.df) == "Code_min"] <- "TrapID_min"


b.fc.df <- subset(b.fc.df, select = c("Landscape", 
                                      "Point",
                                      "TrapID_min", 
                                      "perFC_200", 
                                      "perFC_400", 
                                      "perFC_800", 
                                      "perFC_1000",
                                      "perFC_3000"))

####################################### MANIP BEETLE ABUNDANCE DATA 

# subset forest captures from beetle biodiversity data
b.bdv.fcap.df1 <- subset ( b.bdv.df1 , Habitat == "M" )

# merge TrapID_min codes with TrapIDs in beetle_bd1.df
b.bdv.df <- merge ( b.bdv.fcap.df1 , TrapID.key.df, by = "TrapID")

# convert to long format, summing captures at trap level in "value"
b.bdv.df_melt = melt(b.bdv.df, id.vars=c(1:16,79), measure.vars=17:78, variable_name="code")

# rename "value" to "N.trap" (captures at trap level)
names(b.bdv.df_melt)[match("value", names(b.bdv.df_melt))]="N.trap" 

# get rid of observations that are 0 bc of missing trap 
b.bdv.df_melt1 <- subset(b.bdv.df_melt, b.bdv.df_melt$missing==0)
b.bdv.df_melt1$Lat <- as.numeric(b.bdv.df_melt1$Lat)
b.bdv.df_melt1$Long <- as.numeric(b.bdv.df_melt1$Long)

### summarize abundance per point

# for each species within each point, sum up N over all traps, and stores in 'N.trap' column
b.sum.df = ddply(b.bdv.df_melt1, .(Landscape, Point, code, TrapID_min), summarize, N = sum(N.trap), Lat = mean(Lat), Long = mean(Long))

# determine presence/absence for a species at a point and assigns 1 (present) or 0 (absent) to 'presence' column
b.sum.df$presence = ifelse ( b.sum.df$N >= 1 , 1 , 0 )

# create df of species codes of all beetle species captured
b.spec.capd <- data.frame (unique (b.sum.df$code))
colnames( b.spec.capd ) <- c( "code" )

####################################### CLEAN UP BEETLE TRAIT DATA

# Subset beetle traits of interest
b.tr.df <- b.tr.df1[,c(6:9,12,14,20)]

# subset captured species from trait dataframe
b.tr.df <- merge (b.spec.capd, b.tr.df , by = "code")

# In Diel, consolidate NAs and assign to true categorical <NA>
is.na (b.tr.df$Diel) [b.tr.df$Diel =="Na"] <- TRUE
is.na (b.tr.df$Diel) [b.tr.df$Diel =="NA"] <- TRUE

# In Diel, Replace "D/N" with "N" 
b.tr.df$Diel[b.tr.df$Diel =="D/N"] <- "N"

# set Diel as a factor
b.tr.df$Diel <- factor(b.tr.df$Diel)


# In Nest, assign NA to true categorical <NA>
is.na (b.tr.df$Nest) [b.tr.df$Nest == "NA"] <- TRUE

# set Nest as a factor
b.tr.df$Nest <- factor(b.tr.df$Nest)

# clean up levels of Biogeo_assoc variable
is.na (b.tr.df$Biogeo_assoc)[b.tr.df$Biogeo_assoc == "?" ] <- TRUE
is.na (b.tr.df$Biogeo_assoc)[b.tr.df$Biogeo_assoc == "NA" ] <- TRUE
b.tr.df$Biogeo_assoc[b.tr.df$Biogeo_assoc == "Cerrado, Pantanal and disturbed areas in the MA and Amazon"] <- "G"
b.tr.df$Biogeo_assoc[b.tr.df$Biogeo_assoc == "MA and transition to Cerrado"] <- "T"
b.tr.df$Biogeo_assoc[b.tr.df$Biogeo_assoc == "MA, Cerrado, Chaco, Pampa"] <- "G"
b.tr.df$Biogeo_assoc[b.tr.df$Biogeo_assoc == "Widely distributed "] <- "G"
b.tr.df$Biogeo_assoc[b.tr.df$Biogeo_assoc == "C"] <- "CER"
b.tr.df$Biogeo_assoc[b.tr.df$Biogeo_assoc == "Cerrado, tending to cerradão"] <- "CER"
b.tr.df$Biogeo_assoc[b.tr.df$Biogeo_assoc == "MA and transition to Cerrado and disturbed areas in the MA"] <- "T"
b.tr.df$Biogeo_assoc[b.tr.df$Biogeo_assoc == "Widely distributed in open areas"] <- "G"
b.tr.df$Biogeo_assoc[b.tr.df$Biogeo_assoc == "Ca"] <- "CER"
b.tr.df$Biogeo_assoc[b.tr.df$Biogeo_assoc == "MA, Cerrado"] <- "T"

# group CER associated beetles into T bc only one species
b.tr.df$Biogeo_assoc[b.tr.df$Biogeo_assoc == "CER"] <- "T"

# set biogeo to factor
b.tr.df$Biogeo_assoc <- factor(b.tr.df$Biogeo_assoc)

# make BM_mean_mg variable numeric
b.tr.df$BM_mean_mg <- as.double ( paste ( b.tr.df$BM_mean_mg ) )
b.tr.df$BM_mean_mg <- round (b.tr.df$BM_mean_mg, digits = 3)

# clean up Hab_pref factor levels
b.tr.df$Hab_pref[b.tr.df$Hab_pref == "Forest specialist" ] <- "C"
b.tr.df$Hab_pref[b.tr.df$Hab_pref == " Forest specialist" ] <- "C"
b.tr.df$Hab_pref[b.tr.df$Hab_pref == "Forest specialist" ] <- "C"
b.tr.df$Hab_pref[b.tr.df$Hab_pref == "Forest specialist " ] <- "C"
b.tr.df$Hab_pref[b.tr.df$Hab_pref == "Campo cerrado" ] <- "O"
b.tr.df$Hab_pref[b.tr.df$Hab_pref == "Campo limpo" ] <- "O"
b.tr.df$Hab_pref[b.tr.df$Hab_pref == "Campo rupestre" ] <- "O"
b.tr.df$Hab_pref[b.tr.df$Hab_pref == "Forest" ] <- "C"
b.tr.df$Hab_pref[b.tr.df$Hab_pref == "Generalist" ] <- "G"
b.tr.df$Hab_pref[b.tr.df$Hab_pref == "Generalista" ] <- "G"
b.tr.df$Hab_pref[b.tr.df$Hab_pref == "Open environ. specialist" ] <- "O"
b.tr.df$Hab_pref[b.tr.df$Hab_pref == "Prefers forest" ] <- "C"
b.tr.df$Hab_pref[b.tr.df$Hab_pref == "Prefers open areas" ] <- "O"
is.na (b.tr.df$Hab_pref)[b.tr.df$Hab_pref == "NA" ] <- TRUE

# set Hab_pref to factor
b.tr.df$Hab_pref <- factor(b.tr.df$Hab_pref)

# clean up Diet var
is.na (b.tr.df$Diet)[b.tr.df$Diet == "?" ] <- TRUE
is.na (b.tr.df$Diet)[b.tr.df$Diet == "NA" ] <- TRUE

# set Diet var as factor
b.tr.df$Diet <- factor(b.tr.df$Diet)
b.tr.df$Diet <- droplevels ( b.tr.df$Diet )

############################################## MANIP BEETLE LANDSCAPE COMPOSITION DATA ############################################## 


### Create single df for beetle land comp data

# average TE per point per landscape
#mean.b.te.df = ddply(b.te.df, .(Landscape, Point, TrapID_min ), summarize, TE = mean(TE))

# merge TE data with FC data at point level by TrapID_min
#b.lcomp.df2 = merge(mean.b.te.df, mean.b.fc.df,   by= "TrapID_min")

# subset relevant columns
#b.lcomp.df1 <- b.lcomp.df2[,c(1:4,7,8)]

# rename columns appropriately
#b.lcomp.df <- reshape::rename(b.lcomp.df1, replace = c("Landscape.x" = "Landscape", "Point.x"= "Point"))

# subset percent forest cover per point from b.landcomp.df
#b.perFC.df = subset(b.lcomp.df , select = c("TrapID_min" , "Per_area"))


############################################## MERGE BEETLE TRAIT & ABUNDANCE DATA AT POINT LEVEL ############################################## 

# merge b.sum.df and b.tr.df by species code
b.tr.bdv.df = merge ( b.sum.df , b.tr.df ,   by = "code" )

# make "Point" variable a factor
b.tr.bdv.df$Point <- factor ( paste (b.tr.bdv.df$Point) )

summary (b.tr.bdv.df)


############################################## MERGE BEETLE  TRAIT & ABUNDANCE WITH BEETLE LANDCOMP DATA AT POINT LEVEL ######################### 

# merge a.tr.bdf.df with percent forest cover per point data 
b.df3 <- merge ( b.fc.df, b.tr.bdv.df)

names(b.fc.df)

# make 'Lanscape' in b.df3 a factor
b.df3$Landscape <- as.factor ( b.df3$Landscape )

# assign b.df3 to b.df 
b.df <- b.df3


############################################## FINAL CLEANUP OF b.df ############################################## 
head (b.df)

# rename columns appropriately 
b.df <- reshape::rename ( b.df, replace = c ("TrapID_min" = "TrapID" , 
                                             "code" = "Species" , 
                                             "presence" = "Presence" , 
                                             "Biogeo_assoc" = "Biogeo" , 
                                             "BM_mean_mg" = "BM" , 
                                             "Hab_pref" = "Hab" ))
                                        
# set Lat and Long as numerical variables
b.df$Lat <- as.numeric ( paste ( b.df$Lat ) )
b.df$Long <- as.numeric ( paste ( b.df$Long ) )

# determine which species were captured fewer than three times across all landscapes
b.spec.sumCaps <- ddply(b.df, .(Species), summarize, N.tot = sum(N))

# merge N.tot values for each species with entries in b.df
b.df <- merge (b.df , b.spec.sumCaps , by = "Species")

# drop unused levels of Species factor in b.df
b.df$Species <- droplevels(b.df$Species)

### calculate how many landscapes each species was found in

# create b.df variant with no species/point pairings where no capture occurred
b.df.caps <- subset ( b.df , b.df$Presence == 1)

ls.count.df <-ddply ( b.df.caps, .(Species , Landscape) , 
                      summarise, N.total = mean(N.tot))

ls.count.df <- data.frame (table (ls.count.df$Species , ls.count.df$Landscape))
colnames(ls.count.df) <- c("Species" , "Landscape", "Presence")

# create df with num of landscapes each species was cap'd in
b.sp.ls.count.df <-ddply ( ls.count.df, .(Species) , 
                           summarise, LS.total = sum(Presence))

# merge b.sp.ls.count.df with b.df
b.df <- merge (b.df, b.sp.ls.count.df)

### calculate how many points each species was found at
b.df.caps <- subset ( b.df , b.df$Presence == 1)

pt.count.df <-ddply ( b.df.caps, .(Species , TrapID) , 
                      summarise, N.total = mean(N.tot))

pt.count.df <- data.frame (table (pt.count.df$Species , pt.count.df$TrapID))
colnames(pt.count.df) <- c("Species" , "TrapID", "Presence")

# create df with column containing number of points each species was sighted in
b.sp.pt.count.df <-ddply ( pt.count.df, .(Species) , 
                           summarise, Pts.total = sum(Presence))

# merge point count column with b.df
b.df <- merge (b.df, b.sp.pt.count.df)


### export final b.df object
# save b.df as RDS file
saveRDS ( b.df , file = "b.df.rds" ) 


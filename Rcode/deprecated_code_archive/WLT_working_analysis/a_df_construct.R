#########################
# AVIAN DATA ORGANIZATION
#########################

# PURPOSE: Organize avian abundance, trait, biogeographic and landscape
# composition/structure data into single data frame (a.df) for 
# subsequent analyses.

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

# load packages

library(readxl)
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
library(xlsx)



########################## SET AVIAN DATA ##########################

# read avian biodiversity data
a.bdv.df = data.frame(read_excel("Bird_species_abundance_18_10_16.xlsx", sheet = 1, col_names = TRUE ))

# read avian trait data
a.tr.df1 = data.frame(read_excel( "Bird_species_traits_21_11_16.xlsx" , sheet = 1 , col_names = TRUE ))

# read avian biogeographic association data from bird encyclo
a.biogeo.df1 <- data.frame(read_excel ("avian_encyclopedia_19_1_17.xls", sheet = 1, col_names = TRUE))


########################## LANDSCAPE DATA ##########################

# key for TrapID_min codes
TrapID.key.df<- data.frame(read_excel("Trap_Code_Key.xlsx", sheet = 1, col_names = T))

# read avian FC data
a.landcomp.df <- read.csv("Landscape_composition_pitfalls.csv", header = T)


############################################## MANIP AVIAN ABUNDANCE DATA ############################################## 
  
# convert to long form, summing counts of each species at each trap location 
a.bdv.df_melt <- melt ( a.bdv.df , id.vars = c ( 1:5 ) , measure.vars = 9:232 , variable_name = "code" , header = TRUE )

# rename "value" column to "N.count"
names ( a.bdv.df_melt ) [ match ("value", names(a.bdv.df_melt ) ) ] = "N.count" 

# Replace underscores and spaces in "code" column (species names) with periods
a.bdv.df_melt$code <- gsub ( "_" , "." , a.bdv.df_melt$code )
a.bdv.df_melt$code <- gsub ( " " , "." , a.bdv.df_melt$code )

### summarize abundance per point per landscape, per habitat use subset

# for each species within each point, sum up N over all traps and write to new column "N"
a.sum.df = ddply (a.bdv.df_melt , .( Landscape , Habitat , Point , code , TrapID_min ) , summarize , N = sum ( N.count ) )

#  creates a presence/absence column per species per point and determines presence/absence for all species/point [s] and 
a.sum.df$presence = ifelse ( a.sum.df$N >= 1 , 1 , 0 )

# subset forest captures
a.sum.fcaps.df <- subset ( a.sum.df , Habitat == "M" )

# TEMP: determine how many avian species were sampled
b.cap.summary <- ddply ( a.sum.fcaps.df, .( code ) , summarize , N.total = sum ( N ) )

length((which(b.cap.summary$N.total !=0)))

############################################## MANIP AVIAN TRAIT DATA ############################################## 

# In avian trait dataframe, rename species "Code" column to "code"
names(a.tr.df1)[1] <- "code"

# replace underscores with periods in species code column
a.tr.df1$code <- gsub ( "_" , "." , a.tr.df1$code )
a.tr.df1$code <- gsub ( " " , "." , a.tr.df1$code )

# Subset avian traits of interest
a.tr.df <- a.tr.df1 [ , c("code", "Species", "English.name", 
                          "Body.size.g", "Eggs_minmax_mean", 
                          "Nesting" , "diet_plasticity",  
                          "diet_primary", "hab_assoc") ]



# create granivore/herbivore level of diet_primary factor (GH)
levels(a.tr.df$diet_primary) <- c(levels(a.tr.df$diet_primary), "GH")

# assign granivores to GH level
a.tr.df$diet_primary[a.tr.df$diet_primary == "GRA"] <-"G_H"

# assign herbivores to GH level
a.tr.df$diet_primary[a.tr.df$diet_primary == "HER"] <- "G_H"

# Assign carnivores to NA (because so few of them)
a.tr.df$diet_primary[a.tr.df$diet_primary == "CAR"] <- NA

# make diet_primary a factor
a.tr.df$diet_primary <- factor(a.tr.df$diet_primary)

# NEST: Group open and semi-open species
a.tr.df$Nesting[a.tr.df$Nesting == "SO"] <- "OP"

# make Nesting a factor
a.tr.df$Nesting <- factor(a.tr.df$Nesting)

# Change levels of "Nesting" factor
is.na ( a.tr.df$Nesting ) [ a.tr.df$Nesting == "" ]  <- TRUE

# Drop extranneous levels of Nesting factor
a.tr.df$Nesting <- droplevels ( a.tr.df$Nesting )


# ############################################## MERGE AVIAN ABUNDANCE & TRAIT DATA ON POINT LEVEL ################################# 
# # Identify which avian species are in a.tr.df but not a.sum.fcaps.df (i.e. non-overlap between dfs)
# a.tr.bdv.antiJoin <- anti_join(a.sum.fcaps.df, a.tr.df, by = "code")
# unique (a.tr.bdv.antiJoin$code)


# Merge avian abundance and trait data
a.tr.bdv.df = merge ( a.sum.fcaps.df , a.tr.df ,  by = "code" )


# summarize by species codes
#ddply (a.tr.bdv.df, .(code), summarize, bm = mean (Body.size.g) )


############################################## MANIP AVIAN LANDSCAPE COMPOSITION DATA ############################################## 

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
a.lc.df_merge.df <- merge(fc.pt.200m.df, fc.pt.400m.df)
a.lc.df_merge.df <- merge(a.lc.df_merge.df, fc.pt.800m.df)
a.lc.df_merge.df <- merge(a.lc.df_merge.df, fc.pt.1000m.df)
a.lc.df_merge.df <- merge(a.lc.df_merge.df, fc.pt.3000m.df)

# rename Code_min col name
#names(a.lc.df_merge.df)[3] <- "TrapID_min"



############################################## MERGE AVIAN ABUNDANCE & TRAIT DATA WITH LANDCOMP DATA AT POINT LEVEL ############################################## 

# merge dfs
a.df2 <- merge ( a.tr.bdv.df , a.lc.df_merge.df)


# subset relevant columns
a.df <- a.df2 [ , c("TrapID_min", "code", "Landscape", "Habitat", 
                     "Point" , "N", "presence", "Species",
                     "Body.size.g", 
                     "Eggs_minmax_mean", "Nesting", 
                     "diet_plasticity", 
                     "diet_primary", "hab_assoc", 
                     "Point", 
                    "perFC_200",
                    "perFC_400",
                    "perFC_800",
                    "perFC_1000",
                    "perFC_3000")]

# set variables in a.df as factors
a.df$Species <- as.factor ( a.df$Species )
a.df$English.name <- as.factor ( a.df$English.name )
a.df$code <- as.factor ( a.df$code )

# create new "species.genus" column to correspond with a.biogeo.df
a.df$species.genus <- gsub (" " ,  "." , a.df$Species)


############################################## MERGE a.df WITH AVIAN BIOGEO ASSOC DATA  ##############################################

# rename "code" column of a.biogeo.df to "species.genus"
a.biogeo.df <- reshape::rename (a.biogeo.df1, replace = c("code" = "species.genus"))

#subset relevant columns of a.biogeo.df
a.biogeo.df <- subset ( a.biogeo.df, select = c("species.genus" , "NHAB" , "NZOO" , "biogeo_assoc.v1"))

# TEMP: identify non-overlap 
# antiJoin <- anti_join( a.df, a.biogeo.df , by = "species.genus")
# head(antiJoin)
# ddply (antiJoin, .(species.genus), summarize, N = sum (N))
# 
# ####

# Merge a.df with biogeo assoc data
a.df <- merge ( a.df , a.biogeo.df)

# Biogeographic association: Group transition species with generalists
a.df$biogeo_assoc.v1[a.df$biogeo_assoc.v1 == "T"] <- "GEN"

# TEMP: view a.df with each row as single species 
#all.spec.df <- ddply (a.df, .(a.df$species.genus, a.df$biogeo_assoc), summarize, bm = mean (Point) )

# rename NHAB and NZOO columns in a.df
a.df <- reshape::rename (a.df, replace = c("NHAB" = "hab_plasticity" , "NZOO"= "biogeo_plasticity", 
                                           "biogeo_assoc.v1" = "biogeo_assoc" ))


############################################## MERGE LAT LONG DATA WITH a.df ############################################## 

coord.df <- readRDS("trap_coords.df.rds")
names(coord.df)[1] <- "TrapID_min"

# Merge coord.df with a.df on TrapID_min
a.df <- merge ( coord.df , a.df )

############################################## FINAL CLEANUP OF a.df ############################################## 

# Subset columns of interest
a.df <- subset ( a.df, select = 
c ( "TrapID_min", "code", "species.genus", "Landscape" , "Point" , "N" , 
    "presence" , "Nesting" , "biogeo_assoc" , 
    "biogeo_plasticity" , "Body.size.g" , 
  "Eggs_minmax_mean" , "diet_plasticity" , 
  "diet_primary", "hab_assoc" , "hab_plasticity" , 
  "Lat" , "Long" , "perFC_200", "perFC_400", 
  "perFC_800", "perFC_1000", "perFC_3000" ) ) 

# Rename columns appropriately
a.df <- reshape::rename ( a.df, replace = c 
                          ("TrapID_min" = "TrapID" , "code" = "Species" , "presence" = "Presence" , "Nesting" = "Nest" , 
                          "biogeo_assoc" = "Biogeo" , "biogeo_plasticity" = "BiogeoPlas" , "Body.size.g" = "BM" , 
                          "hab_assoc" = "Habitat" , "Eggs_minmax_mean" = "Fecund" , "diet_plasticity" = "DietPlas" , 
                          "hab_plasticity" = "HabPlas" , "diet_primary" = "Diet"))

# set factor variables to factors
a.df$TrapID <- factor(a.df$TrapID)
a.df$Landscape <- factor(a.df$Landscape)
a.df$Point <- factor(a.df$Point)
a.df$Biogeo <- factor(a.df$Biogeo)
a.df$Diet <- factor(a.df$Diet)
a.df$Habitat <- factor(a.df$Habitat)


# Redefine levels of "Biogeo" factor in a.df
levels ( a.df$Biogeo ) <- c ( levels (a.df$Biogeo) , "G" )

# Assign entries with Biogeo level GEN to G
a.df$Biogeo [a.df$Biogeo == "GEN"] <- "G"
 
# Drop unused levels of Biogeo factor
a.df$Biogeo <- droplevels ( a.df$Biogeo )

# Redefine levels of "Habitat" factor in a.df
levels ( a.df$Habitat ) <- c ( levels (a.df$Habitat) , "C" )

# Assign entries with Habitat level GEN to G
a.df$Habitat [a.df$Habitat == "F"] <- "C"

# Drop unused levels of Biogeo factor
a.df$Habitat <- droplevels ( a.df$Habitat )


# remove duplicate row entries in a.df
a.df <- unique ( a.df )

# determine which species were captured fewer than three times across all landscapes
a.spec.sumCaps <- ddply(a.df, .(Species), summarize, N.tot = sum(N))

# merge N.tot values for each species with entries in a.df
a.df <- merge (a.df , a.spec.sumCaps)


### calculate how many landscapes each species was found in

# create a.df variant with no species/point pairings where no capture occurred
a.df.caps <- subset ( a.df , a.df$Presence == 1)

ls.count.df <-ddply ( a.df.caps, .(Species , Landscape) , 
                      summarise, N.total = mean(N.tot))

ls.count.df <- data.frame (table (ls.count.df$Species , ls.count.df$Landscape))
colnames(ls.count.df) <- c("Species" , "Landscape", "Presence")

# create df with num of landscapes each species was cap'd in
a.sp.ls.count.df <-ddply ( ls.count.df, .(Species) , 
                           summarise, LS.total = sum(Presence))

# merge b.sp.ls.count.df with b.df
a.df <- merge (a.df, a.sp.ls.count.df)


### calculate how many points each species was found at
a.df.caps <- subset ( a.df , a.df$Presence == 1)

pt.count.df <-ddply ( a.df.caps, .(Species , TrapID) , 
                      summarise, N.total = mean(N.tot))

pt.count.df <- data.frame (table (pt.count.df$Species , pt.count.df$TrapID))
colnames(pt.count.df) <- c("Species" , "TrapID", "Presence")

# create df with column containing number of points each species was sighted in
a.sp.pt.count.df <-ddply ( pt.count.df, .(Species) , 
                           summarise, Pts.total = sum(Presence))

# merge point count column with a.df
a.df <- merge (a.df, a.sp.pt.count.df)

# change Species column to contain full species.genus strings
a.species.codes <- a.df$Species 
a.df$Species <- a.df$species.genus

# put species codes in a.species.codes columns of df
names(a.df)[names(a.df)=="species.genus"] <- "a.species.codes"
a.df$a.species.codes <- a.species.codes


# save a.df as RDS object
saveRDS ( a.df , file = "a.df.rds" )

a.df$Diet




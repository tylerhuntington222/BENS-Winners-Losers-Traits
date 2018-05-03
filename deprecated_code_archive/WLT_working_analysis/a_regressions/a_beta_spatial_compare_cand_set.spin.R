#####################################################################################################################################
################################## AVIAN N~FC REGRESSION MODELS COMPARISON #########################################################
#####################################################################################################################################

# PURPOSE: compare AVIA responses (N.point ~ FC) at different buffer sizes around forest points (800m, 3000m )

# AUTHORS: Tyler Huntington, Elizabeth Nichols, Larissa Boesing, Jean Paul Metzger


# set working directory
setwd("/Users/tyler/Dropbox/Interface/Dados/Non-Geospatial/Winners_Losers_Traits/Analysis/Data")

# load base beetle dataframe
a.df <- readRDS("a.df.rds")

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
library (lme4)



################################### READ IN & ORGANIZE DATA 

# subset a.df for species that were captured 3 or more times
a.df <- subset ( a.df , a.df$N.tot > 3 )

# subset a.df for species that were sighted at three distinct points
a.df <- subset ( a.df , a.df$Pts.total > 3 )


# subset out species with unclear data trends
a.df <- subset (a.df, a.df$Species!= "Ramphas.toco")
a.df <- subset (a.df, a.df$Species!= "Phae.eury")
a.df <- subset (a.df, a.df$Species!= "Phae.pret")
a.df <- subset (a.df, a.df$Species!= "Myiornis.aur")
a.df <- subset (a.df, a.df$Species!= "Malac.stri")


a.df$Species <- factor(a.df$Species)

a.df$Species <- droplevels ( a.df$Species  )

# create df with unique entry for each beetle species and associated attributes
a.sp.df <- ddply ( a.df, .(Species , Nest , Biogeo , BiogeoPlas , BM , 
                           Fecund, DietPlas, Diet, Habitat, HabPlas, LS.total, Pts.total), 
                   summarize, N.tot = mean(N.tot))

### Run N ~ FC multivariate regressions for each species using FC vals at each buffer size as IVs

# initialize blank df for top spatial models per species
a_spatial_mod.df <- data.frame("Species" = factor(), 
                               "Top_Spatial_Mod" = factor())

# initialize model counter vars
b1000_count <- 0
b800_count <- 0
b3000_count <- 0
b400_count <- 0
b200_count <- 0

for (i in a.sp.df$Species) {
  
  # initialize species row to be inserted in top model df
  a_species_row <- NULL
  
  # create df with entries as all captures of a particular species
  a.par.df <- subset (a.df , a.df$Species == i) 
  
  # regress N ~ FC for particular species at multiple spatial scales
  a.lm.par <- glmer (a.par.df$N ~ a.par.df$perFC_200
                     + a.par.df$perFC_400
                     + a.par.df$perFC_800
                     + a.par.df$perFC_1000
                     + a.par.df$perFC_3000
                     + (1|a.par.df$Landscape), family = poisson)
  
  options(na.action = "na.fail")
  a.par.spatial.mod.comp <- dredge(a.lm.par, beta = c("none", "sd", "partial.sd"), evaluate = TRUE,
                                   m.lim = c(1, 1), rank = "AICc", fixed = NULL)
  
  a.par.top.models <- get.models(a.par.spatial.mod.comp, subset = delta < 2)
  
  cand_set_codes <- rownames (summary(a.par.top.models))
  
  # KEY FOR MODEL NUMBERS IN DREDGE TABLE
  # Model 2 = 1000m buffer
  # Model 17 = 800m buffer
  # Model 5 = 3000m buffer
  # Model 9 = 400m buffer
  # Model 3 = 200m buffer
  
  
  # use model code to assign top model for particular species
  if ( "2" %in% cand_set_codes) {
    b1000_count <- b1000_count + 1
  }  
  if ("17" %in% cand_set_codes) {
    b800_count <- b800_count + 1
  } 
  if ("5" %in% cand_set_codes) {
    b3000_count <- b3000_count + 1
  } 
  if ("9" %in% cand_set_codes) {
    b400_count <- b400_count + 1
  } 
  if ("3" %in% cand_set_codes) {
    b200_count <- b200_count + 1
  }
  
}

################################################ ANALYZE RESULTS OF SPATIAL SCALE COMPS

# plot distribution of top spatial models across beetle species
par(mar = c(5, 5, 5, 5) + 0.1)
counts = c(b200_count, b400_count, b800_count, b1000_count, b3000_count)
labels = c("200", "400", "800", "1000", "3000")
barplot <- barplot(counts, space = 0.25, names.arg = labels, ylim = c(0 ,50),
                   main = "Birds: Candidate Set Approach\n for Determining Spatial Scale",
                   xlab = "Buffer Radius (m)",
                   ylab = "# Species w/ Model in Candidate Set")
text(barplot, counts, labels = counts, pos = 3)
barplot

# print counts as a table
print(data.frame(b200_count, b400_count, b800_count, b1000_count, b3000_count))


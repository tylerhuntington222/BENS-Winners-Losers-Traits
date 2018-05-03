#####################################################################################################################################
################################## REGRESSION ANALYSIS: TRAIT CORRELATES OF RESPONSES TO FOREST COVER IN DUNG BEETLES ###############
#####################################################################################################################################

# OVERVIEW:

# Step 1:
# GLMM: N ~ FC | Landscape
# (BEETLE ABUNDANCES SUMMED AT THE POINT LEVEL)
# (FC CALCULATED AT SPECIES-SPECIFIC SPATIAL SCALE DETERMINED BY INITIAL MODELING STEP)

# Step 2:
# PGLS: Response ~ Traits | Phylo

#####################################################################################################################################
options(max.print = 1000)
dev.off()

# set working directory
setwd("/Users/nicholslab/Dropbox/Interface/Dados/Non-Geospatial/Winners_Losers_Traits/Analysis/Data")

list.files()


# load libraries
library(knitr)
library(phytools)
library(car)
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
library(Hmisc)



############### BEETLES ############### 

# Read in beetle abundance and trait data
b.df <- readRDS ( "b.df.rds" )

names(b.df)
########## PREP BEETLE DATA FOR MODELING ########## 

# create df with unique entry for each beetle species and associated attributes
b.sp.df <- ddply ( b.df, .(Species , Diel , Nest , Biogeo , BM , Hab , Diet, Pts.total, LS.total, N.tot), 
                   summarise, mean(N.tot))

# subset for species that were captured at at three distinct points
b.sp.df <- subset ( b.sp.df , b.sp.df$Pts.total >= 3 )

# subset for species that were captured in at least two different landscapes
b.sp.df <- subset ( b.sp.df , b.sp.df$LS.total >= 2 )

# re-factor species column of df
b.sp.df$Species <- factor(b.sp.df$Species)


########## PLOT N ~ FC FOR EACH BEETLE SPECIES

# # make plot matrix of individual beetle species N ~ FC
# 
# for (i in b.sp.df$Species) {
# 
#   # create df with entries as all captures of a particular species
#   b.par.df <- subset (b.df , b.df$Species == i) # omitting/including zeros doesn't seem to affect model
# 
#   # plot N ~ FC for particular species
#   plot(b.par.df$perFC_200, b.par.df$N, main = i)
# 
#   print (i)
# 
# }

########## REMOVE BEETLE SPECIES THAT DO NOT MEET MODELING CRITERIA ########## 


# remove species with insufficient values to run
b.sp.df = b.sp.df[!b.sp.df$Species == "Canthon.aff.semiopacus", ]
# b.sp.df = b.sp.df[!b.sp.df$Species == "Deltochilum.dentipes", ]
# b.sp.df = b.sp.df[!b.sp.df$Species == "Uroxys.sp1", ]
# b.sp.df = b.sp.df[!b.sp.df$Species == "Dichotomius.depressicollis", ]
# b.sp.df = b.sp.df[!b.sp.df$Species == "Scybalocanthon.nigriceps", ]

# # remove species with fewer than 2 captures at a point
b.sp.df = b.sp.df[!b.sp.df$Species == "Dichotomius.sp5", ]
b.sp.df = b.sp.df[!b.sp.df$Species == "Uroxys.sp4", ]
b.sp.df = b.sp.df[!b.sp.df$Species == "Coprophanaeus.cerberus", ]

# # remove species with very few captures
# b.sp.df = b.sp.df[!b.sp.df$Species == "Dichotomius.sp2", ]
# b.sp.df = b.sp.df[!b.sp.df$Species == "Canthon.aff.luctuosus", ]

# check that it worked
b.sp.df$Species = factor(b.sp.df$Species) # weirdly important step. 



########## DETERMINE BEST SPATIAL SCALE PER BEETLE SPECIES ########## 

# initialize blank df for top spatial models per species
b_sp_spatial_scales.df <- data.frame("Species" = factor(), 
                               "TopSpatialScale" = factor())


# initialize model counter vars
b1000_count_b_tm <- 0
b800_count_b_tm <- 0
b3000_count_b_tm <- 0
b400_count_b_tm <- 0
b200_count_b_tm <- 0


for (i in levels(b.sp.df$Species)) {
  
  # initialize species row to be inserted in top model df
  b_species_row <- NULL
  
  # create df with entries as all captures of a particular species
  b.par.df <- subset (b.df , b.df$Species == i) 
  
  # model N ~ FC for particular species at all spatial scales
  b.par.spatial.glmm <- glmer (b.par.df$N ~ b.par.df$perFC_200
                     + b.par.df$perFC_400
                     + b.par.df$perFC_800
                     + b.par.df$perFC_1000
                     + b.par.df$perFC_3000
                     + (1|b.par.df$Landscape), family = poisson)
  
  options(na.action = "na.fail")
  b.par.spatial.mod.comp <- dredge(b.par.spatial.glmm, beta = c("none", "sd", "partial.sd"), evaluate = TRUE,
                                   m.lim = c(1,1), rank = "AICc", fixed = NULL)
  
  b.par.top.model <- get.models(b.par.spatial.mod.comp, subset = delta == 0)
  
  top_model_code <- rownames (summary(b.par.top.model))
  
  # KEY FOR MODEL NUMBERS IN DREDGE TABLE
  # Model 2 = 1000m buffer
  # Model 17 = 800m buffer
  # Model 5 = 3000m buffer
  # Model 9 = 400m buffer
  # Model 3 = 200m buffer
  
  
  # use model code to assign top model for particular species
  if (top_model_code == "2") {
    top_mod <- 1000
    b1000_count_b_tm = b1000_count_b_tm + 1 
    
  } else if (top_model_code == "17") {
    top_mod <- 800
    b800_count_b_tm = b800_count_b_tm + 1 
    
  } else if (top_model_code == "5") {
    top_mod <- 3000
    b3000_count_b_tm = b3000_count_b_tm + 1
    
  } else if (top_model_code == "9") {
    top_mod <- 400
    b400_count_b_tm = b400_count_b_tm + 1
    
  } else if (top_model_code == "3") {
    top_mod <- 200
    b200_count_b_tm = b200_count_b_tm + 1
  }
  
  b_species_row <- data.frame(factor(i), factor(top_mod))
  b_sp_spatial_scales.df <- rbind (b_sp_spatial_scales.df, b_species_row)
  #print(b_sp_spatial_scales.df)
}

# rectify column names of b_spatial_mod.df
colnames(b_sp_spatial_scales.df) <- c("Species", "TopSpatialScale")

# merge top spatial model df with b.sp.df
b.sp.df <- merge(b.sp.df, b_sp_spatial_scales.df)

# merge top spatial model df with b.df
b.df <- merge(b.df, b_sp_spatial_scales.df) 


# plot distribution of top spatial models across beetle species
par(mar = c(4, 4, 4, 4) + 0.8)
b.counts = c(b200_count_b_tm, b400_count_b_tm, b800_count_b_tm, b1000_count_b_tm, b3000_count_b_tm)
labels = c("200", "400", "800", "1000", "3000")
title = ("Distribution of Optimal Spatial Scales for Modeling\
Dung Beetle Species Responses to Habitat Loss")
b.spatial_barplot <- barplot(b.counts, space = 0.15, names.arg = labels, ylim = c(0 , 15),
                   main = title,
                   cex.main = 1.3,
                   xlab = "Buffer Radius (m)",
                   ylab = "Number of Species",
                   cex.names = 1,
                   cex.axis = 1,
                   cex.lab = 1.2,
                   col = "turquoise3",
                   border = F)

text(b.spatial_barplot, b.counts, labels = b.counts, cex=0.9, pos = 3)

# # set up dataframe for plotting
# b.sp_scale_counts.df <-  data.frame(b200_count_b_tm, 
#                                       b400_count_b_tm, 
#                                       b800_count_b_tm, 
#                                       b1000_count_b_tm, 
#                                       b3000_count_b_tm)
# 
# b.sp_scale_counts.df <- data.frame(t(b.sp_scale_counts.df))
# rownames(b.sp_scale_counts.df) <- c('200', '400', '800', '1000', '3000')
# b.sp_scale_counts.df$b.spatial.scale <- factor(rownames(b.sp_scale_counts.df))
# 
# ggplot (b.sp_scale_counts.df, aes(x = reorder(b.spatial.scale, c(1,2,3,4,5)), y = b.counts)) +
#   geom_bar(stat = 'identity', fill = 'springgreen3') +
#   # title
#   ggtitle("Distribution of Optimal Spatial Scales for Modeling
# Dung Beetle Species Responses to Habitat Loss") +
#   theme(plot.title = element_text(hjust = 0.5, family = "sans", face = "bold")) +
#   theme(plot.title = element_text(size=20)) + 
#   # axes 
#   xlab("Buffer Radius Used to Calculate Forest Cover (m)") +
#   theme(axis.title.x = element_text(size=15, family = "sans")) +
#   ylab("Number of Species") +
#   theme(axis.title.y = element_text(size=15, family = "sans")) 
#   
# 
#   
#   
# 
# 
# str(b.sp_scale_counts.df)
# 
# b.sp_scale_counts.df$b.spatial.scale


########## CALCULATE SPECIES SPECIFIC N ~ FC RESPONSES FOR BEETLES

# create blank beetle FC response df
b.responses.df <- data.frame(
  Species = character(), 
  Intercept = double(),
  Inter.error = double(),
  Inter.t.val = double(),
  Inter.p.val = double(), 
  b.beta = double(),
  b.beta.error = double (),
  b.beta.t.val = double(),
  b.beta.p.val = double(),
  stringsAsFactors = F
)

# initialize row 1 as first slot of responses.df to be filled
row = 1

# iterate response calculation over all beetle species and store in b.responses.df
for (i in b.sp.df$Species) {
  
  # create df with entries as all captures of a particular (par) species
  b.par.df <- subset (b.df , b.df$Species == i) 
  
  # scale longitude to the mean
  b.par.df$Long.mscale = scale(b.par.df$Long, center = TRUE, scale = FALSE)
  
  # model N ~ FC at species specific spatial scale
  if (b.par.df$TopSpatialScale == "200"){
    
    # scale FC at 200m scale to mean
    b.par.df$perFC_200.mscale = scale(b.par.df$perFC_200, center = TRUE, scale = FALSE)
    
    # model N ~ FC for particular species
    b.glmm.par <- glmer (N ~ perFC_200.mscale  + Long.mscale + (1|Landscape), data = b.par.df, family = poisson)
    
  } else if (b.par.df$TopSpatialScale == "400"){
    
    # scale FC at 400m scale to mean
    b.par.df$perFC_400.mscale = scale(b.par.df$perFC_400, center = TRUE, scale = FALSE)
    
    # model N ~ FC for particular species
    b.glmm.par <- glmer (N ~ perFC_400.mscale  + Long.mscale + (1|Landscape), data = b.par.df, family = poisson)
    
  } else if (b.par.df$TopSpatialScale == "800") {
    
    # scale FC at 800m scale to mean
    b.par.df$perFC_800.mscale = scale(b.par.df$perFC_800, center = TRUE, scale = FALSE)
    
    # model N ~ FC for particular species
    b.glmm.par <- glmer (N ~ perFC_800.mscale  + Long.mscale + (1|Landscape), data = b.par.df, family = poisson)
    
  } else if (b.par.df$TopSpatialScale == "1000") {
    
    # scale FC at 1000m scale to mean
    b.par.df$perFC_1000.mscale = scale(b.par.df$perFC_1000, center = TRUE, scale = FALSE)
    
    # model N ~ FC for particular species
    b.glmm.par <- glmer (N ~ perFC_1000.mscale  + Long.mscale + (1|Landscape), data = b.par.df, family = poisson)
    
  } else if (b.par.df$TopSpatialScale == "3000") 
    {
    
    # scale FC at 3000m scale to mean
    b.par.df$perFC_3000.mscale = scale(b.par.df$perFC_3000, center = TRUE, scale = FALSE)
    
    # model N ~ FC for particular species
    b.glmm.par <- glmer (N ~ perFC_3000.mscale  + Long.mscale + (1|Landscape), data = b.par.df, family = poisson)}

  
  # add row to response df with model output for particular species
  b.responses.df[row, 1] <- as.character(i)
  b.responses.df[row,2:5] <- summary (b.glmm.par)$coefficients[1,]
  b.responses.df[row,6:9] <- summary (b.glmm.par)$coefficients[2,]
  
  # increase row counter by one 
  row <- row + 1
}

# round p-values in response df
b.responses.df$b.beta.p.val <- round(b.responses.df$b.beta.p.val, 3)

# create histogram of responses
b.winner_loser.hist <- hist(b.responses.df$b.beta, 
                        main = "Distribution of Beetle Responses to Habitat Loss",
                        xlab = "Response Coefficient",
                        ylab = "Frequency",
                        ylim = c(0,10),
                        col = "gray",
                        breaks = 30) 
abline(v = 0, col = "red")


text(0,9,  paste("Winners            Losers"), cex=1.1)



################################################ MERGE N ~ FC RESPONSES WITH b.sp.df

b.sp.df <- merge(b.responses.df, b.sp.df)

# make row names of b.sp.df be species names
rownames(b.sp.df) <- b.sp.df$Species

################################################ MODEL RESPONSES ~ TRAITS | PHYLO

#### Import phylo tree object 
b.tre <- read.tree("b.sp.tre")
b.tre$node.label <- NULL

# generate comparative phylogenetic autocorrelation object
b.comp <- comparative.data(b.tre, b.sp.df, Species, vcv=TRUE)

# build pgls model
b_glmm.pgls <- pgls(b.beta ~ BM + Diel + Nest + Diet, b.comp )

# print model summary
summary(b_glmm.pgls)
plot(b_glmm.pgls)

# + Hab  --- generates error
# + Biogeo -- generates error


# generate all possible models and rank by AIC

b_pgls_mods.table <- dredge(b_glmm.pgls, beta = c("none", "sd", "partial.sd"), evaluate = TRUE,
                                rank = "AICc", fixed = NULL)

print(b_pgls_mods.table)


plot(log(b.sp.df$BM) ~ b.sp.df$b.beta)

losers <- subset(b.sp.df, b.sp.df$b.beta > 0)
meanBM_losers <- mean(losers$BM)

winners <- subset(b.sp.df, b.sp.df$b.beta < 0)
meanBM_winners <- mean(winners$BM)

# create boxplot of winners v. losers log(BM)
plot.new()
b.bm_boxplot <- boxplot( log(winners$BM), log(losers$BM),
                         main = "Boxplot of Winner vs. Loser Dung Beetle Body Masses",
                         ylab = "log(BM)",
                         names = c("Winners", "Losers"),
                         col = c("skyblue2", "brown2"))



############### BIRDS ############### 

# Read in bird abundance and trait data
a.df <- readRDS ( "a.df.rds" )

names(a.df)

########## MANIPULATE BIRD DATA FOR MODELING ##########

# create df with unique entry for each bird species and associated attributes
a.sp.df <- ddply ( a.df, .(Species, Nest, Biogeo, BiogeoPlas, BM, Fecund, DietPlas, Diet,
                           Habitat , HabPlas, Pts.total, LS.total), summarise, N.tot = mean(N.tot))

# subset for species that were captured at at three distinct points
a.sp.df <- subset ( a.sp.df , a.sp.df$Pts.total >= 3 )

# subset for species that were captured in at least two different landscapes
a.sp.df <- subset ( a.sp.df , a.sp.df$LS.total >= 2 )

# remove bird species with fewer than two sightings at a point
a.sp.df = a.sp.df[!a.sp.df$Species == "Tyrannus.melancholicus", ]
a.sp.df = a.sp.df[!a.sp.df$Species == "Trogon.surrucura", ]
a.sp.df = a.sp.df[!a.sp.df$Species == "Thlypopsis.sordida", ]
a.sp.df = a.sp.df[!a.sp.df$Species == "Syndactyla.rufosuperciliata", ]
a.sp.df = a.sp.df[!a.sp.df$Species == "Stephanophorus.diadematus", ]
a.sp.df = a.sp.df[!a.sp.df$Species == "Sclerurus.scansor", ]
a.sp.df = a.sp.df[!a.sp.df$Species == "Pachyramphus.viridis", ]
a.sp.df = a.sp.df[!a.sp.df$Species == "Myiozetetes.similis", ]
a.sp.df = a.sp.df[!a.sp.df$Species == "Hylophilus.poicilotis", ]
a.sp.df = a.sp.df[!a.sp.df$Species == "Herpsilochmus.rufimarginatus", ]
a.sp.df = a.sp.df[!a.sp.df$Species == "Herpetotheres.cachinnans", ]
a.sp.df = a.sp.df[!a.sp.df$Species == "Elaenia.flavogaster", ]
a.sp.df = a.sp.df[!a.sp.df$Species == "Dryocopus.lineatus", ]
a.sp.df = a.sp.df[!a.sp.df$Species == "Cranioleuca.pallida", ]
a.sp.df = a.sp.df[!a.sp.df$Species == "Columbina.talpacoti", ]
a.sp.df = a.sp.df[!a.sp.df$Species == "Colonia.colonus", ]
a.sp.df = a.sp.df[!a.sp.df$Species == "Colaptes.melanochloros", ]
a.sp.df = a.sp.df[!a.sp.df$Species == "Cnemotriccus.fuscatus", ]
a.sp.df = a.sp.df[!a.sp.df$Species == "Aphantochroa.cirrochloris", ]

# remove species with insufficient values to run
a.sp.df = a.sp.df[!a.sp.df$Species == "Calli.ameth", ]
a.sp.df = a.sp.df[!a.sp.df$Species == "Gral.varia", ]
a.sp.df = a.sp.df[!a.sp.df$Species == "Pseroco.dec", ]
a.sp.df = a.sp.df[!a.sp.df$Species == "Trogon.rufus", ]

# re-factor Species var
a.sp.df$Species = factor(a.sp.df$Species) 




# re-factor species column of df
a.sp.df$Species <- factor(a.sp.df$Species)


# ########## PLOT N ~ FC FOR EACH BIRD SPECIES ########## 
# 
# # make plot matrix of individual bird species N ~ FC
# 
# for (i in a.sp.df$Species) {
# 
#   # create df with entries as all captures of a particular species
#   a.par.df <- subset (a.df , a.df$Species == i) # omitting/including zeros doesn't seem to affect model
# 
#   # plot N ~ FC for particular species
#   plot(a.par.df$perFC_200, a.par.df$N, main = i)
# 
#   print (i)
# 
# }
# 
# ########## DETERMINE BEST SPATIAL SCALE PER BIRD SPECIES

# initialize blank df for top spatial models per species
a_sp_spatial_scales.df <- data.frame("Species" = factor(), 
                                     "TopSpatialScale" = factor())


# initialize model counter vars
b1000_count_a_tm <- 0
b800_count_a_tm <- 0
b3000_count_a_tm <- 0
b400_count_a_tm <- 0
b200_count_a_tm <- 0


for (i in levels(a.sp.df$Species)) {
  
  # initialize species row to be inserted in top model df
  a_species_row <- NULL
  
  # create df with entries as all captures of a particular species
  a.par.df <- subset (a.df , a.df$Species == i) 
  
  # model N ~ FC for particular species at all spatial scales
  a.par.spatial.glmm <- glmer (a.par.df$N ~ a.par.df$perFC_200
                               + a.par.df$perFC_400
                               + a.par.df$perFC_800
                               + a.par.df$perFC_1000
                               + a.par.df$perFC_3000
                               + (1|a.par.df$Landscape), family = poisson)
  
  options(na.action = "na.fail")
  a.par.spatial.mod.comp <- dredge(a.par.spatial.glmm, beta = c("none", "sd", "partial.sd"), evaluate = TRUE,
                                   m.lim = c(1,1), rank = "AICc", fixed = NULL)
  
  a.par.top.model <- get.models(a.par.spatial.mod.comp, subset = delta == 0)
  
  top_model_code <- rownames (summary(a.par.top.model))
  
  # KEY FOR MODEL NUMBERS IN DREDGE TABLE
  # Model 2 = 1000m buffer
  # Model 17 = 800m buffer
  # Model 5 = 3000m buffer
  # Model 9 = 400m buffer
  # Model 3 = 200m buffer
  
  
  # use model code to assign top model for particular species
  if (top_model_code == "2") {
    top_mod <- 1000
    b1000_count_a_tm = b1000_count_a_tm + 1 
    
  } else if (top_model_code == "17") {
    top_mod <- 800
    b800_count_a_tm = b800_count_a_tm + 1 
    
  } else if (top_model_code == "5") {
    top_mod <- 3000
    b3000_count_a_tm = b3000_count_a_tm + 1
    
  } else if (top_model_code == "9") {
    top_mod <- 400
    b400_count_a_tm = b400_count_a_tm + 1
    
  } else if (top_model_code == "3") {
    top_mod <- 200
    b200_count_a_tm = b200_count_a_tm + 1
  }
  
  a_species_row <- data.frame(factor(i), factor(top_mod))
  a_sp_spatial_scales.df <- rbind (a_sp_spatial_scales.df, a_species_row)
  #print(a_sp_spatial_scales.df)
}

# rectify column names of a_spatial_mod.df
colnames(a_sp_spatial_scales.df) <- c("Species", "TopSpatialScale")

# merge top spatial model df with a.sp.df
a.sp.df <- merge(a.sp.df, a_sp_spatial_scales.df)

# merge top spatial model df with a.df
a.df <- merge(a.df, a_sp_spatial_scales.df) 


# plot distribution of top spatial models across bird species
par(mar = c(4, 4, 4, 4) + 0.8)
a.counts = c(b200_count_a_tm, b400_count_a_tm, b800_count_a_tm, b1000_count_a_tm, b3000_count_a_tm)
labels = c("200", "400", "800", "1000", "3000")
title = ("Distribution of Optimal Spatial Scales for Modeling\
Bird Species Responses to Habitat Loss")
a.sp_scale_dist.plot <- barplot(a.counts, space = 0.15, names.arg = labels, ylim = c(0 , 30),
                   main = title,
                   cex.main = 1.3,
                   xlab = "Buffer Radius (m)",
                   ylab = "Number of Species",
                   cex.names = 1,
                   cex.axis = 1,
                   cex.lab = 1.2,
                   col = "springgreen3",
                   border = F)

text(a.sp_scale_dist.plot , a.counts, labels = a.counts, cex=0.9, pos = 3)


########## CALCULATE SPECIES SPECIFIC N ~ FC RESPONSES FOR birds

# create blank bird FC response df
a.responses.df <- data.frame(
  Species = character(), 
  Intercept = double(),
  Inter.error = double(),
  Inter.t.val = double(),
  Inter.p.val = double(), 
  a.beta = double(),
  a.beta.error = double (),
  a.beta.t.val = double(),
  a.beta.p.val = double(),
  stringsAsFactors = F
)


# initialize row 1 as first slot of responses.df to be filled
row = 1

# iterate response calculation over all bird species and store in a.responses.df
for (i in a.sp.df$Species) {
  
  # create df with entries as all captures of a particular (par) species
  a.par.df <- subset (a.df , a.df$Species == i) 
  
  # scale longitude to the mean
  a.par.df$Long.mscale = scale(a.par.df$Long, center = TRUE, scale = FALSE)
  
  # model N ~ FC at species specific spatial scale
  if (a.par.df$TopSpatialScale == "200"){
    
    # scale FC at 200m scale to mean
    a.par.df$perFC_200.mscale = scale(a.par.df$perFC_200, center = TRUE, scale = FALSE)
    
    # model N ~ FC for particular species
    a.glmm.par <- glmer (N ~ perFC_200.mscale  + Long.mscale + (1|Landscape), data = a.par.df, family = poisson)
    
  } else if (a.par.df$TopSpatialScale == "400"){
    
    # scale FC at 400m scale to mean
    a.par.df$perFC_400.mscale = scale(a.par.df$perFC_400, center = TRUE, scale = FALSE)
    
    # model N ~ FC for particular species
    a.glmm.par <- glmer (N ~ perFC_400.mscale  + Long.mscale + (1|Landscape), data = a.par.df, family = poisson)
    
  } else if (a.par.df$TopSpatialScale == "800") {
    
    # scale FC at 800m scale to mean
    a.par.df$perFC_800.mscale = scale(a.par.df$perFC_800, center = TRUE, scale = FALSE)
    
    # model N ~ FC for particular species
    a.glmm.par <- glmer (N ~ perFC_800.mscale  + Long.mscale + (1|Landscape), data = a.par.df, family = poisson)
    
  } else if (a.par.df$TopSpatialScale == "1000") {
    
    # scale FC at 1000m scale to mean
    a.par.df$perFC_1000.mscale = scale(a.par.df$perFC_1000, center = TRUE, scale = FALSE)
    
    # model N ~ FC for particular species
    a.glmm.par <- glmer (N ~ perFC_1000.mscale  + Long.mscale + (1|Landscape), data = a.par.df, family = poisson)
    
  } else if (a.par.df$TopSpatialScale == "3000") 
  {
    
    # scale FC at 3000m scale to mean
    a.par.df$perFC_3000.mscale = scale(a.par.df$perFC_3000, center = TRUE, scale = FALSE)
    
    # model N ~ FC for particular species
    a.glmm.par <- glmer (N ~ perFC_3000.mscale  + Long.mscale + (1|Landscape), data = a.par.df, family = poisson)}
  
  
  # add row to response df with model output for particular species
  a.responses.df[row, 1] <- as.character(i)
  a.responses.df[row,2:5] <- summary (a.glmm.par)$coefficients[1,]
  a.responses.df[row,6:9] <- summary (a.glmm.par)$coefficients[2,]
  
  # increase row counter by one 
  row <- row + 1
}

# round p-values in response df
a.responses.df$a.beta.p.val <- round(a.responses.df$a.beta.p.val, 3)

# merge avian responses with a.sp.df
a.sp.df <- merge(a.responses.df, a.sp.df)

# add column to a.sp.df to ID winners and losers
a.sp.df$Winner_Loser = NULL
a.sp.df$Winner_Loser[a.sp.df$a.beta > 0] <- "Loser"
a.sp.df$Winner_Loser[a.sp.df$a.beta < 0] <- "Winner"
a.sp.df$Winner_Loser <- factor(a.sp.df$Winner_Loser)

# make row names of a.sp.df be species names
rownames(a.sp.df) <- a.sp.df$Species

# create histogram of responses
plot.new()
a.winner_loser.hist <- hist (a.sp.df$a.beta, main = "Distribution of Avian Species Responses to Habitat Loss",
      xlab = "Response Coefficient",
      ylab = "Frequency",
      breaks = 30,
      ylim = c(0,20),
      xlim = c(-0.6,0.2),
      col = "grey")

abline(v = 0, col = "red")

text(0,18,  paste("Winners            Losers"), cex=1.1)



########## MODEL RESPONSES ~ TRAITS | PHYLO ##########

#### Import phylo tree object 
a.trees <- read.nexus("avian_complete_species.tre")

# calculate consensus tree
a.tre <- ls.consensus(a.trees, start=NULL, tol=1e-12, quiet=FALSE)

# alternative function for calculating consensus tree -- VERY long runtime
#a.tre <- averageTree(a.trees, start=NULL, tol=1e-12, quiet=FALSE)

# arrange a.sp.df for generating comp phylo autocorrelation object
a.sp.df$Species <- gsub(".", "_", a.sp.df$Species, fixed = T)

# set row names to species names in format that matches tree object
rownames(a.sp.df) <- paste(a.sp.df$Species)

# generate comparative phylogenetic autocorrelation object
a.comp <- comparative.data(a.tre, a.sp.df, Species, vcv=TRUE, force.root = T)

# build pgls model
a_glmm.pgls <- pgls(a.beta ~ Nest + BiogeoPlas + BM + Diet + Habitat, a.comp)
                    
# print model summary
summary(a_glmm.pgls)
plot(a_glmm.pgls)

# generate all possible models and rank by AIC

a_pgls_mods.table <- dredge(a_glmm.pgls, 
                            beta = c("none", "sd", "partial.sd"), 
                            evaluate = TRUE,
                            rank = "AICc", 
                            fixed = NULL)

print(a_pgls_mods.table, digits = 2)


# exploratory plots
plot(a.sp.df$Biogeo, a.sp.df$a.beta)

# group winners and losers based on beta
losers <- subset(a.sp.df, a.sp.df$a.beta > 0)
meanBM_losers <- mean(losers$BM)

winners <- subset(a.sp.df, a.sp.df$a.beta < 0)
meanBM_winners <- mean(winners$BM)

# create boxplot of winners v. losers log(BM)
boxplot(log(winners$BM), log(losers$BM))






  
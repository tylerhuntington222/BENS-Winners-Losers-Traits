##############################################################################
# a+b_trait_response_models_at_spp_specific_spatial_scales.R

# Generates global model of FC responses ~ traits for avian and beetle 
# species using FC responses calculated at species-specific 
# spatial scales (i.e. the
# spatial scale that best modeled N ~ FC for a given species).

# High-level overview of script workflow:

# for (model_system in [avian_species, beetle_species]):

# Step 1:
# GLMM: N ~ FC | Landscape
# (ABUNDANCES SUMMED AT THE POINT LEVEL)
# (FC CALCULATED AT SPECIES-SPECIFIC SPATIAL SCALE)

# Step 2:
# PGLS: Response ~ Traits | Phylo
# b_df_construct.R

# INPUT: Cleaned avian and beetle data.

# OUTPUT: Numerous plots visualizing model results.  

# AUTHOR: Tyler Huntington
###############################################################################

# load libraries
library (reshape)
library (plyr)
library (reshape2)
library (dplyr)
library (rgdal)
library (rgeos)
library (spdep)
library (plyr)
library (ggplot2)
library (caper)
library (ape)
library (MuMIn)
library (lme4)
library (Hmisc)
library (phytools)

# Set WLT Data folder as working directory.
# Filepaths are relative, so running the below code
# without modification
# on any machine should set the wd properly.
tryCatch (
  {
    this.dir <- dirname(parent.frame(2)$ofile)
  },
  error = function(x) {
    this.dir <<- dirname(rstudioapi::getActiveDocumentContext()$path)
  }
)
setwd(this.dir)
setwd("../..")

################################# BEETLES ################################# 

# Read in beetle abundance and trait data
b.df <- readRDS ( "Intermediates/b.df.rds" )

###### MANIPULATE DATA FOR MODELING ######

# create df with unique entry for each beetle species and associated attributes
b.sp.df <- ddply (b.df, 
                  .(Species , Diel , Nest , Biogeo , BM , Hab , 
                    Diet, Pts.total, LS.total, N.tot), 
                  summarise, mean.N.tot = mean(N.tot))


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

########## DETERMINE BEST SPATIAL SCALE PER BEETLE SPECIES

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
  tryCatch({
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
    b.par.spatial.mod.comp <- dredge(b.par.spatial.glmm, 
                                     beta = c("none", "sd", "partial.sd"), 
                                     evaluate = TRUE,
                                     m.lim = c(1,1), rank = "AICc", 
                                     fixed = NULL)
    
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
    print(b_sp_spatial_scales.df)
  }, 
  error = function(e) {
    print(i)
    print(e)
  }
  )
}

# rectify column names of b_spatial_mod.df
colnames(b_sp_spatial_scales.df) <- c("Species", "TopSpatialScale")

# merge top spatial model df with b.sp.df
b.sp.df <- merge(b.sp.df, b_sp_spatial_scales.df)

# merge top spatial model df with b.df
b.df <- merge(b.df, b_sp_spatial_scales.df) 


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


# factorize species col
b.sp.df$Species = factor(b.sp.df$Species) 

# initialize row 1 as first slot of responses.df to be filled
row = 1

# iterate response calculation over all beetle species and store in b.responses.df
for (i in b.sp.df$Species) {
  tryCatch({
    print(sprintf("Working on species %s", i))
    
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
    
    # increase row counter
    row <- row + 1
  }, error = function(e) {
    print(i)
    print(e)
  })
  
}

# round p-values in response df
b.responses.df$b.beta.p.val <- round(b.responses.df$b.beta.p.val, 3)

# create histogram of responses
b.winner_loser.hist <- hist(b.responses.df$b.beta, 
                            main = "Distribution of Beetle Responses to Habitat Loss",
                            xlab = "Disturbance Response",
                            ylab = "Frequency",
                            ylim = c(0.5,10),
                            col = "turquoise3",
                            breaks = 35,
                            border = "gray") 

abline(v = 0, col = "red")
minor.tick(nx=10, ny=5, tick.ratio=0.3)

text(0,9,  paste("Winners               Losers   "), cex=1.1)


############# MERGE N ~ FC RESPONSES WITH b.sp.df

b.sp.df <- merge(b.responses.df, b.sp.df)

# make row names of b.sp.df be species names
rownames(b.sp.df) <- b.sp.df$Species

############ MODEL RESPONSES ~ TRAITS | PHYLO

#### Import phylo tree object 
b.tre <- read.tree("Data/b.sp.tre")
b.tre$node.label <- NULL

# subset species to model 

# subset for species that were captured at at three distinct points
b.sp.df.mod <- subset ( b.sp.df , b.sp.df$Pts.total >= 3 )

# subset for species that were captured in at least two different landscapes
b.sp.df.mod <- subset ( b.sp.df.mod , b.sp.df$LS.total >= 2)

## subset for species that had > 1 cap for at least one point
# init vector to store results
max.caps <- c(NULL)
# iterate over species to find max caps at a point per species
for (s in b.sp.df.mod$Species) {
  s.df <- b.df[b.df$Species == s, ]
  max.caps.at.pt <- max(s.df$N)
  max.caps <- c(max.caps, max.caps.at.pt)
}
# bind to b.sp.df
b.sp.df.mod$max.pt.caps <- max.caps

b.sp.df.mod <- subset(b.sp.df.mod, b.sp.df.mod$max.pt.caps > 2)

# re-factor species column of df
b.sp.df.mod$Species <- factor(b.sp.df.mod$Species)

# set up vars for modeling
b.sp.df.mod$Biogeo <- as.character(b.sp.df.mod$Biogeo)

# refactor
b.sp.df.mod$Hab <- factor(b.sp.df.mod$Hab)

# generate comparative phylogenetic autocorrelation object
b.comp <- comparative.data(b.tre, b.sp.df.mod, Species, vcv=TRUE)

# BEETLE MODEL
b_glmm.pgls <- pgls(b.beta ~ BM + Nest + Hab + Diet + Biogeo, b.comp )
summary(b_glmm.pgls)


# df <- data.frame(b_glmm.pgls$data)
# 
# write.table(b_glmm.pgls, file = "sumstats.txt", sep = ",", 
#             quote = FALSE, row.names = F)


# generate all possible models and rank by AIC

b_pgls_mods.table <- dredge(b_glmm.pgls, beta = c("none", "sd", "partial.sd"), 
                            evaluate = TRUE,
                            rank = "AICc", fixed = NULL)

# explore corrs between ind vars
# variation in beta?
head(b.sp.df.mod)
plot(b.sp.df.mod$TopSpatialScale , b.sp.df.mod$BM)
plot(b.sp.df.mod$TopSpatialScale, b.sp.df.mod$b.beta)
plot(b.sp.df.mod$Nest, b.sp.df.mod$b.beta)
plot(b.sp.df.mod$BM, b.sp.df.mod$b.beta)

hist(b.sp.df.mod$BM, breaks = 100)
# order this in terms of spatial scale
# but in general, that's interesting?


#-----------------------------------------------------------------------------#
###### BIRDS ######

# OVERVIEW:

# Step 1:
# GLMM: N ~ FC | Landscape
# (BEETLE ABUNDANCES SUMMED AT THE POINT LEVEL)
# (FC CALCULATED AT SPECIES-SPECIFIC SPATIAL SCALE DETERMINED BY INITIAL MODELING STEP)

# Step 2:
# PGLS: Response ~ Traits | Phylo
############

########## BIRDS

# Read in beetle abundance and trait data
a.df <- readRDS ( "Intermediates/a.df.rds" )

# replace abbreviated species name field with full species genus string
a.df$Species <- NA
a.df$Species <- a.df$species.genus

################################## MANIPULATE DATA FOR MODELING ##############

# create df with unique entry for each beetle species and associated attributes
a.sp.df <- ddply (a.df, 
                  .(Species , Nest , Biogeo, 
                    BiogeoPlas, BM, Habitat , N.tot,
                    Diet, Pts.total, LS.total, Fecund), 
                  summarise, mean.N.tot = mean(N.tot))

names(a.sp.df)[names(a.sp.df) == "species.genus"] <- "Species"
a.sp.df$Species <- as.character(a.sp.df$Species)

# subset for species that were captured at at three distinct points
a.sp.df <- subset ( a.sp.df , a.sp.df$Pts.total >= 3 )

# subset for species that were captured in at least two different landscapes
a.sp.df <- subset ( a.sp.df , a.sp.df$LS.total >= 2)

## subset for species that had > 1 cap for at least one point
# init vector to store results
max.caps <- c(NULL)
# iterate over species to find max caps at a pt per species
for (s in a.sp.df$Species) {
  s.df <- a.df[a.df$Species == s, ]
  max.caps.at.pt <- max(s.df$N)
  max.caps <- c(max.caps, max.caps.at.pt)
}
# bind to a.sp.df
a.sp.df$max.pt.caps <- max.caps

a.sp.df <- subset(a.sp.df, a.sp.df$max.pt.caps > 2)

# re-factor species column of df
a.sp.df$Species <- factor(a.sp.df$Species)


########## PLOT N ~ FC FOR EACH AVIAN SPECIES

# # make plot matrix of individual beetle species N ~ FC
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

########## DETERMINE BEST SPATIAL SCALE PER BEETLE SPECIES

# initialize blank df for top spatial models per species
b_sp_spatial_scales.df <- data.frame("Species" = factor(), 
                                     "TopSpatialScale" = factor())


# initialize model counter vars
b1000_count_a_tm <- 0
b800_count_a_tm <- 0
b3000_count_a_tm <- 0
b400_count_a_tm <- 0
b200_count_a_tm <- 0


for (i in levels(a.sp.df$Species)) {
  
  # initialize species row to be inserted in top model df
  b_species_row <- NULL
  
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
  a.par.spatial.mod.comp <- dredge(a.par.spatial.glmm, 
                                   beta = c("none", "sd", "partial.sd"), 
                                   evaluate = TRUE,
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
  
  b_species_row <- data.frame(factor(i), factor(top_mod))
  b_sp_spatial_scales.df <- rbind (b_sp_spatial_scales.df, b_species_row)
  print(b_sp_spatial_scales.df)
}

# rectify column names of b_spatial_mod.df
colnames(b_sp_spatial_scales.df) <- c("Species", "TopSpatialScale")

# merge top spatial model df with a.sp.df
a.sp.df <- merge(a.sp.df, b_sp_spatial_scales.df)

# merge top spatial model df with a.df
a.df <- merge(a.df, b_sp_spatial_scales.df) 



########## CALCULATE SPECIES SPECIFIC N ~ FC RESPONSES FOR BEETLES

# create blank beetle FC response df
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

# # remove species with insufficient values to run
# a.sp.df = a.sp.df[!a.sp.df$Species == "Canthon.aff.semiopacus", ]
# a.sp.df = a.sp.df[!a.sp.df$Species == "Deltochilum.dentipes", ]
# a.sp.df = a.sp.df[!a.sp.df$Species == "Uroxys.sp1", ]
# a.sp.df = a.sp.df[!a.sp.df$Species == "Dichotomius.depressicollis", ]
# a.sp.df = a.sp.df[!a.sp.df$Species == "Scybalocanthon.nigriceps", ]

# # remove species with fewer than 2 captures at a point
# a.sp.df = a.sp.df[!a.sp.df$Species == "Dichotomius.sp5", ]
# a.sp.df = a.sp.df[!a.sp.df$Species == "Uroxys.sp4", ]
# a.sp.df = a.sp.df[!a.sp.df$Species == "Coprophanaeus.cerberus", ]


# check that it worked
a.sp.df$Species = factor(a.sp.df$Species) #weirdly important step. 

# initialize row 1 as first slot of responses.df to be filled
row = 1

# iterate response calculation over all beetle species and store in a.responses.df
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
    a.par.df$perFC_3000.mscale = scale(a.par.df$perFC_3000, 
                                       center = TRUE, scale = FALSE)
    
    # model N ~ FC for particular species
    a.glmm.par <- glmer (N ~ perFC_3000.mscale  + Long.mscale + 
                           (1|Landscape), data = a.par.df, family = poisson)}
  
  
  # add row to response df with model output for particular species
  
  a.responses.df[row, 1] <- as.character(i)
  a.responses.df[row,2:5] <- summary (a.glmm.par)$coefficients[1,]
  a.responses.df[row,6:9] <- summary (a.glmm.par)$coefficients[2,]
  
  # increase row counter
  row <- row + 1
}

# round p-values in response df
a.responses.df$a.beta.p.val <- round(a.responses.df$a.beta.p.val, 3)


############### MERGE N ~ FC RESPONSES WITH a.sp.df

a.sp.df <- merge(a.responses.df, a.sp.df)

# make row names of a.sp.df be species names
rownames(a.sp.df) <- a.sp.df$Species

################################################ MODEL RESPONSES ~ TRAITS | PHYLO

#### Import phylo tree object 
a.trees <- read.nexus("Data/avian_complete_species.tre")

# calculate consensus tree
a.tre <- ls.consensus(a.trees, start=NULL, tol=1e-12, quiet=FALSE)

# arrange a.sp.df for generating comp phylo autocorrelation object
a.sp.df$Species <- gsub(".", "_", a.sp.df$Species, fixed = T)

# set row names to species names in format that matches tree object
rownames(a.sp.df) <- paste(a.sp.df$Species)
a.sp.df$Species <- as.factor(a.sp.df$Species)

# generate comparative phylogenetic autocorrelation object
a.comp <- comparative.data(a.tre, a.sp.df, Species, vcv=TRUE, force.root = T)

# BIRD MODEL
a_glmm.pgls <- pgls(a.beta ~ BM + Nest + Habitat + Diet + BiogeoPlas + TopSpatialScale, a.comp )

# explore model simplification
a_pgls_mods.table <- dredge(a_glmm.pgls, beta = c("none", "sd", "partial.sd"), 
                            evaluate = TRUE,
                            rank = "AICc", fixed = NULL)

write.csv(a_pgls_mods.table[1:2,], "a_pgls_mods.table.csv")

# so, we included all the traits. Biogeographic plasticity, diet and spatial scale were the only predictors here. 
# AIC selection suggests that the most parsimonious model is one of biogeographic plasticity alone, as diet doesn't improve model fit to the data


# Plot Diet and BiogeoPlast together (or any categorical and continuous variable)
intercepts <- c(coef(a_glmm.pgls2)["(Intercept)"] ,
                coef(a_glmm.pgls2)["(Intercept)"] + coef(a_glmm.pgls2)["DietG_H"],
                coef(a_glmm.pgls2)["(Intercept)"] + coef(a_glmm.pgls2)["DietINS"],
                coef(a_glmm.pgls2)["(Intercept)"] + coef(a_glmm.pgls2)["DietNEC"] )

lines.df <- data.frame(intercepts = intercepts,
                       slopes = rep(coef(a_glmm.pgls)["BiogeoPlas"], 4),
                       Diet = levels(a.comp$data$Diet))

qplot(x = BiogeoPlas, y = a.beta, color = Diet, data = a.comp$data) + 
  geom_abline(aes(intercept = intercepts, 
                  slope = slopes, 
                  color = Diet), data = lines.df)

a.comp 

summary(a_glmm.pgls)


# generate all possible models and rank by AIC

a_pgls_mods.table <- dredge(a_glmm.pgls, beta = c("none", "sd", "partial.sd"), 
                            evaluate = TRUE,
                            rank = "AICc", fixed = NULL)


###### PLOTTING ######
library(ggplot2)
library(reshape2)

# set text sizes
axis.text <- 15
axis.title <- 17
legend.text <- 15
legend.title <- 16.5

sp.scales <- c("200", "400", "800", "1000", "3000")
sp.scales <-  factor(sp.scales)
sp.scales <- factor(sp.scales, levels(sp.scales)[c(2,4,5,1,3)])

a.sp.counts <- c(
  b200_count_a_tm, 
  b400_count_a_tm,
  b800_count_a_tm,
  b1000_count_a_tm,
  b3000_count_a_tm
)

b.sp.counts <- c(
  b200_count_b_tm, 
  b400_count_b_tm,
  b800_count_b_tm,
  b1000_count_b_tm,
  b3000_count_b_tm
)

ab.sp.counts <- data.frame(a.sp.counts, b.sp.counts, sp.scales)

ab.sp.counts.melts <- melt(ab.sp.counts, id.vars = 'sp.scales')

# set colors for plotting
red <- "#F8766D"
blue <- "#00BFC4"


### barplot of top spatial scales
sp.scales <- ggplot(ab.sp.counts.melts, 
                    aes(x = sp.scales, y = value, fill = variable)) +
  geom_bar(stat = 'identity', position='dodge') +
  scale_fill_manual("Model System",
                    labels = c("Birds", "Dung Beetles"),
                    values=c(red, blue)) +
  xlab("Spatial Scale (m)") +
  theme(axis.text.x=element_text(size=axis.text), 
        axis.text.y=element_text(size=axis.text),
        axis.title.x = element_text(size = axis.title),
        axis.title.y = element_text(size = axis.title),
        legend.text = element_text(size = legend.text),
        legend.title = element_text(size = legend.title),
        legend.position = "right", 
        legend.direction = "vertical",
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black")) +
  ylab("Number of Species")

tiff(file = "../Figures/spatial_scales.tiff", res = 400, width = 9, height = 5,
     units = "in", compression = "lzw")
plot(sp.scales)

dev.off()

#### histograms of disturbance responses


# BEETLES 
# BIRDS
# create histogram of responses


tiff(file = "../Figures/b_wl_hist.tiff", res = 400, width = 8, height = 5.75,
     units = "in", compression = "lzw")
b.winner_loser.hist <- hist(b.responses.df$b.beta,
                            main = NULL, 
                            xlab = "Sensitivity to native forest loss",
                            ylab = NULL,
                            cex.lab = 1.7,
                            cex.axis = 1.7,
                            ylim = c(0.3,10),
                            xlim = c(-0.53, 0.2),
                            col = blue,
                            breaks = 35,
                            border = "gray") 
abline(v = 0, col = "black")
minor.tick(nx=10, ny=5, tick.ratio=0.3)
text(0, 9,  paste("Winners               Losers   "), cex=1.7)
dev.off()





# BIRDS
# create histogram of responses
tiff(file = "../Figures/a_wl_hist.tiff", res = 400, width = 8, height = 5.75,
     units = "in", compression = "lzw")
a.winner_loser.hist <- hist(a.responses.df$a.beta, 
                            main = NULL,
                            xlab = "Sensitivity to native forest loss",
                            ylab = NULL,
                            cex.lab = 1.7,
                            cex.axis = 1.7,
                            ylim = c(0.5, 15),
                            xlim = c(-0.53, 0.2),
                            col = red,
                            breaks = 35,
                            border = "gray") 

abline(v = 0, col = "black")
minor.tick(nx=10, ny=5, tick.ratio=0.3)
text(0,14,  paste("Winners               Losers   "), cex=1.7)

dev.off()


## create scatterplot of avian disturbance responses ~ biogeo plasticity

tiff(file = "Figures/a_biogeo_scatter.tiff", res = 400, width = 8, 
     height = 5.75,
     units = "in", compression = "lzw")
a.biogeo.scatter <- ggplot(a.sp.df, aes(x = BiogeoPlas, y = a.beta)) + 
  geom_point(stat = "identity", col = red, size = 3) + 
  geom_smooth(method=lm,   # Add linear regression line
              se=T) +
  scale_y_continuous(breaks=seq(-0.5,0.2, .1)) +
  scale_x_continuous(breaks=seq(0, 15, 2)) +
  xlab("Biogeographic plasticity") +
  theme(axis.text.x=element_text(size=axis.text), 
        axis.text.y=element_text(size=axis.text),
        axis.title.x = element_text(size = axis.title),
        axis.title.y = element_text(size = axis.title),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(), 
        axis.line = element_line(colour = "black")) +
  ylab("Sensitivity to native forest loss")

a.biogeo.scatter

dev.off()

## create boxplot for beetle hab pref

# subset data for plot
b.sp.df.hab <- b.sp.df[!(is.na(b.sp.df$Hab)),]

# start graphics driver
tiff(file = "Figures/b_habitat_boxplot.tiff", res = 400, width = 8, 
     height = 6.25,
     units = "in", compression = "lzw")

b.hab.boxplot <- ggplot(b.sp.df.hab, aes(x = Hab, y = b.beta)) + 
  stat_boxplot(geom ='errorbar', width = 0.25) +
  geom_boxplot(fill = blue) +
  geom_hline(yintercept = 0, lty = "dashed", col = "grey") +
  scale_x_discrete(labels = c("Forest Specialists", "Generalists")) +
  scale_y_continuous(breaks=seq(-0.5,0.2, 0.05)) +
  xlab("Habitat Preference") +
  theme(axis.text.x=element_text(size=14), 
        axis.text.y = element_text(size=axis.text, margin = margin(l = 8)),
        axis.title.x = element_text(size = axis.title, margin = margin(t=20)),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.title.y = element_text(size = axis.title)) +
  ylab("Disturbance Sensitivity")

b.hab.boxplot

dev.off()

red <- "#F8766D"
blue <- "#00BFC4"





## create boxplot for bird diet
a.sp.df$Diet

# subset data for plot
a.sp.df.diet <- a.sp.df[!(is.na(a.sp.df$Diet)),]
a.sp.df.diet <- a.sp.df[!(a.sp.df$a.beta< -0.3),]

# start graphics driver
tiff(file = "Figures/a_diet_boxplot.tiff", res = 400, width = 8, 
     height = 6.25,
     units = "in", compression = "lzw")
a_diet_boxplot <- ggplot(a.sp.df.diet, aes(x = Diet, y = a.beta)) + 
  stat_boxplot(geom ='errorbar', width = 0.25) +
  geom_boxplot(fill = red) +
  geom_hline(yintercept = 0, lty = "dashed", col = "grey") +
  scale_x_discrete(labels = c("Frugivore" , "Granivore/Herbivore", "Insectivore", "Nectivore")) +
  scale_y_continuous(breaks=seq(-0.5,0.2, 0.05)) +
  xlab("Diet") +
  theme(axis.text.x=element_text(size=14), 
        axis.text.y = element_text(size=axis.text, margin = margin(l = 8)),
        axis.title.x = element_text(size = axis.title, margin = margin(t=20)),
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        axis.title.y = element_text(size = axis.title)) +
  ylab("Disturbance Sensitivity")
a_diet_boxplot
dev.off()






red <- "#F8766D"
blue <- "#00BFC4"


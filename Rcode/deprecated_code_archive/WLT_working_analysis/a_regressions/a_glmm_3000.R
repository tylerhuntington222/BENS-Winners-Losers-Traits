#####################################################################################################################################
################################## REGRESSION ANALYSIS: TRAIT CORRELATES OF RESPONSES TO FOREST COVER IN BIRDS ###############
#####################################################################################################################################

# OVERVIEW:

# Step 1:
# GLMM: N ~ FC | Landscape
# (AVIAN ABUNDANCES SUMMED AT THE POINT LEVEL)
# (FC CALCULATED AT 3000m BUFFER AROUND POINTS)

# Step 2:
# PGLS: Response ~ Traits | Phylo

#####################################################################################################################################

# set working directory
setwd("/Users/nicholslab/Dropbox/Interface/Dados/Non-Geospatial/Winners_Losers_Traits/Analysis/Data")

list.files()

# load libraries

library(phytools)
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


# load beetle abundance and traits dataframe
a.df <- readRDS ( "a.df.rds" )

names(a.df)
options(max.print=3000000)
################################## MANIPULATE DATA FOR MODELING ################################## 

# subset a.df for species that were captured 3 or more times
a.df <- subset ( a.df , a.df$N.tot > 3 )

# subset a.df for species that were sighted at at three distinct points
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
a.sp.df <- ddply ( a.df, .(Species , species.genus, Nest , Biogeo , BiogeoPlas , BM , 
                           Fecund, DietPlas, Diet, Habitat, HabPlas, LS.total, Pts.total), 
                   summarize, N.tot = mean(N.tot))



################################################ Plot N ~ FC @ 3000m for all each beetle species

# make plot matrix of individual beetle species N ~ FC
par(mfrow=c(2,3))
par (pty = "s")
par(mar=c(2, 2,  2,  2))

for (i in a.sp.df$Species) {
  
  # create df with entries as all captures of a particular species
  a.par.df <- subset (a.df , a.df$Species == i) # omitting/including zeros doesn't seem to affect model
  
  # plot N ~ FC for particular species
  plot(a.par.df$perFC_3000, jitter(a.par.df$N, 0.5), main = i)
  
  print (i)
  
}

######################### EXPLORATORY ANALYSIS OF VARIABLE DISTRIBUTIONS

table(a.sp.df$Nest:a.sp.df$HabPlas)

hist (a.sp.df$BM)
hist(log(a.sp.df$BM))

plot(a.sp.df$Nest)

# Check contingency table and interaction between Diet & Nest
table (a.sp.df$Diet, a.sp.df$Nest)
chisq.test (a.sp.df$Diet, a.sp.df$Nest)

# Check contingency table and interaction between Diel & Biogeo
table (a.sp.df$Habitat, a.sp.df$Biogeo)
chisq.test (a.sp.df$Habitat, a.sp.df$Biogeo)

# Check contingency table and interaction between Diel & BM
table (a.sp.df$Diel, b.logBM.bins)
chisq.test (a.sp.df$Diel, b.logBM.bins)

# Check contingency table and interaction between Diel & Hab
table (a.sp.df$Diel, a.sp.df$Hab)
chisq.test (a.sp.df$Diel, a.sp.df$Hab)

# Check contingency table and interaction between Nest & Biogeo
table (a.sp.df$Nest, a.sp.df$Biogeo)
chisq.test (a.sp.df$Nest, a.sp.df$Biogeo)

# Check contingency table and interaction between Nest & Diet
table (a.sp.df$Fecund, a.sp.df$DietPlas)
chisq.test (a.sp.df$Fecund, a.sp.df$BiogeoPlas)

hist(a.sp.df$BiogeoPlas)

# Check contingency table and interaction between Nest & BM
table (a.sp.df$Diel, b.logBM.bins)
chisq.test (a.sp.df$Nest, b.logBM.bins)

# Check contingency table and interaction between Nest & Hab
table (a.sp.df$Nest, a.sp.df$HabPlas)
chisq.test (a.sp.df$Nest, a.sp.df$HabPlas)

plot(a.sp.df$BiogeoPlas, a.sp.df$HabPlas)

# Check contingency table and interaction between Biogeo & Diet
table (a.sp.df$Biogeo, a.sp.df$Diet)
chisq.test (a.sp.df$Biogeo, a.sp.df$Diet)

# Check contingency table and interaction between Biogeo & BM
table (a.sp.df$Biogeo, b.logBM.bins)
chisq.test (a.sp.df$Biogeo, b.logBM.bins )

# Check contingency table and interaction between Biogeo & Hab
table (a.sp.df$Biogeo, a.sp.df$Hab)
chisq.test (a.sp.df$Biogeo, a.sp.df$Hab)

# Check contingency table and interaction between Diet & BM
table (a.sp.df$Diet, b.logBM.bins)
chisq.test (a.sp.df$Diet, b.logBM.bins)

# Check contingency table and interaction between Diet & Hab
table (a.sp.df$Diet, a.sp.df$Hab)
chisq.test (a.sp.df$Diet, a.sp.df$Hab)

# Check contingency table and interaction between BM & Hab
table (b.logBM.bins, a.sp.df$Hab)
chisq.test (b.logBM.bins, a.sp.df$Hab)



#############################

# set plot matrix
par(mfrow=c(3,2))

## Diel
# plot histogram of beetle species by Diel
plot (a.sp.df$Diel, main = "Diel")
summary (a.sp.df$Diel)

## Nest
# plot histogram of beetle species by nesting strategy
plot (a.sp.df$Nest, main = "Nest")
summary (a.sp.df$Nest)

## Diet
# plot histogram of beetle species by diet 
plot (a.sp.df$Diet, main = "Diet")
summary (a.sp.df$Diet)

## Diet Plasticity
# plot histogram of beetle species by diet 
plot (a.sp.df$DietPlas, main = "Diet Plasticity")
summary (a.sp.df$DietPlas)

## Body Mass
# plot histogram of beetle species by logBM
hist (log(a.sp.df$BM), main = "log(BM)")
summary (a.sp.df$BM)

## Fecundity
# plot histogram of beetle species by logBM
hist (log(a.sp.df$Fecund), main = "Fecund")
summary (a.sp.df$Fecund)

## Biogeo association
# plot histogram of avian species by Biogeo
plot (a.sp.df$Biogeo, main = "Biogeo")
summary (a.sp.df$Biogeo)


## Biogeo association
# plot histogram of avian species by Biogeo
plot (a.sp.df$BiogeoPlas, main = "Biogeo Plasticity")
summary (a.sp.df$BiogeoPlas)

## Habitat preference
# plot histogram of avia species by Hab Assoc
plot (a.sp.df$Habitat, main = "Hab Assoc")
summary (a.sp.df$Habitat)

# plot histogram of avian species by Hab plasticity
plot (a.sp.df$HabPlase, main = "Habitat Plasticity")
summary (a.sp.df$HabPlas)




################################################ CALCULATE N ~ FC RESPONSES FOR ALL BEETLE SPECIES

# create blank beetle FC response df
a.3000respFC.df <- data.frame(
  Species = factor(), 
  Intercept = double(),
  Inter.error = double(),
  Inter.t.val = double(),
  Inter.p.val = double(), 
  a.slope_3000 = double(),
  a.slope_3000.error = double (),
  a.slope_3000.t.val = double(),
  a.slope_3000.p.val = double()
)

# # remove species with insufficient values to run
a.sp.df = a.sp.df[!a.sp.df$Species == "Calli.ameth", ]
a.sp.df = a.sp.df[!a.sp.df$Species == "Gral.varia", ]
a.sp.df = a.sp.df[!a.sp.df$Species == "Pseroco.dec", ]
a.sp.df = a.sp.df[!a.sp.df$Species == "Trogon.rufus", ]


# check that it worked
a.sp.df$Species = factor(a.sp.df$Species) #weirdly important step. 
levels(a.sp.df$Species)

# set levels of species var
levels(a.3000respFC.df$Species) <- unique (a.sp.df$Species)

# initialize row 1 as first slot of b.respFC.df to be filled
row = 1

# iterate response calculation over all beetle species and store in b.respFC.df
for (i in a.sp.df$Species) {
  
  # create df with entries as all captures of a particular (par) species
  a.par.df <- subset (a.df , a.df$Species == i) 
  
  # # scale both longitude and forest cover to the mean
  a.par.df$Long.mscale = scale(a.par.df$Long, center = TRUE, scale = FALSE)
  a.par.df$perFC_3000.mscale = scale(a.par.df$perFC_3000, center = TRUE, scale = FALSE)
  
  # model N ~ FC for particular species
  a.glmm.par <- glmer (N ~ perFC_3000.mscale  + Long.mscale + (1|Landscape), data = a.par.df, family = poisson)
  
  
  # add row to response df with model output for particular species
  a.3000respFC.df[row, 1] <- (i)
  a.3000respFC.df[row,2:5] <- summary (a.glmm.par)$coefficients[1,]
  a.3000respFC.df[row,6:9] <- summary (a.glmm.par)$coefficients[2,]
  
  # increase row counter
  row <- row + 1
}


a.3000respFC.df


# 
# # round p-values in response df
a.3000respFC.df$a.slope_3000.p.val <- round(a.3000respFC.df$a.slope_3000.p.val, 3)
# 
# # create histogram of responses
# hist (a.3000respFC.df$a.slope_3000, breaks = 20)
# 
# # subset for species with significant a.slope_3000 coefficients
# a.sigResp <- subset (a.3000respFC.df, a.3000respFC.df$a.slope_3000.p.val < 0.05)
# 
# hist (a.sigResp$a.slope_3000, breaks = 15)
# hist (a.sigResp$a.slope_3000.p.val, breaks = 15)
# 
# 
# ### Export summary of betas generated by this modeling approach
# 
# # subset a.3000respFC.df for a.slope_3000s and p-values
# a.3000respFC.df_betas <- subset(a.3000respFC.df, select = c("Species", "a.slope_3000", "a.slope_3000.p.val"))
# 
# # rename columns to specify modeling approach number
# colnames(b.respFC.df_betas) = c("Species", "a_glmm_3000_a.slope_3000", "b_glmm_3000_p.val")
# 
# #  save as RDS object
# # saveRDS(b.respFC.df_betas, "b.Mod1_betas.df")
# 
# ################################################ MERGE N ~ FC RESPONSES WITH a.sp.df
# 
a.sp.df <- merge(a.3000respFC.df, a.sp.df)


#### Create preliminary traits GLM without accounting for phylo
a_trait_resp_glm <- glm(a.slope_3000 ~ BM + Biogeo + BiogeoPlas + Fecund + DietPlas + Diet + Habitat + HabPlas, data = a.sp.df)


options('na.action' = na.fail)
a_3000_trait_resp_glm.table <- dredge(a_trait_resp_glm, beta = c("none", "sd", "partial.sd"), evaluate = TRUE,
                                      rank = "AICc", fixed = NULL)

names(a.sp.df)

# ################################################ MODEL RESPONSES ~ TRAITS | PHYLO
# 
#### Import phylo tree object 
a.trees <- read.nexus("avian_complete_species.tre")

# calculate consensus tree
a.tre <- ls.consensus(a.trees, start=NULL, tol=1e-12, quiet=FALSE)

# alternative function for calculating consensus tree -- VERY long runtime
#a.tre <- averageTree(a.trees, start=NULL, tol=1e-12, quiet=FALSE)

# arrange a.sp.df for generating comp phylo autocorrelation object
a.sp.df$species.genus <- gsub(".", "_", a.sp.df$species.genus, fixed = T)

# set row names to species names in format that matches tree object
rownames(a.sp.df) <- paste(a.sp.df$species.genus)

# generate comparative phylogenetic autocorrelation object
a.comp_3000 <- comparative.data(a.tre, a.sp.df, species.genus, vcv=TRUE, force.root = T)

# build pgls model
a_glmm_3000.pgls <- pgls(a.slope_3000 ~ BM + Diet + DietPlas + Nest + Biogeo + BiogeoPlas + Fecund + Habitat + HabPlas, a.comp_3000)

# generate all possible models and rank by AIC
a_pgls_3000_mods.table <- dredge(a_glmm_3000.pgls, beta = c("none", "sd", "partial.sd"), evaluate = TRUE,
                                 rank = "AICc", fixed = NULL)

print (a_pgls_3000_mods.table)

summary(a_pgls_3000_mods.table)


# make plots of relationships of interest
plot( a.sp.df$HabPlas, a.sp.df$a.slope_3000)

plot (a.sp.df$BiogeoPlas, a.sp.df$a.slope_3000)

plot (a.sp.df$Diet, a.sp.df$a.slope_3000)



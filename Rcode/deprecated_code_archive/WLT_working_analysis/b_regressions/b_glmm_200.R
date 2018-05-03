#####################################################################################################################################
################################## REGRESSION ANALYSIS: TRAIT CORRELATES OF RESPONSES TO FOREST COVER IN DUNG BEETLES ###############
#####################################################################################################################################

# OVERVIEW:

# Step 1:
# GLMM: N ~ FC | Landscape
# (BEETLE ABUNDANCES SUMMED AT THE POINT LEVEL)
# (FC CALCULATED AT 200m BUFFER AROUND POINTS)

# Step 2:
# PGLS: Response ~ Traits | Phylo

#####################################################################################################################################

# set working directory
setwd("/Users/tyler/Dropbox/Interface/Dados/Non-Geospatial/Winners_Losers_Traits/Analysis/Data")

list.files()


# load libraries
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
b.df <- readRDS ( "b.df.rds" )

names(b.df)
################################## MANIPULATE DATA FOR MODELING ################################## 

# create df with unique entry for each beetle species and associated attributes
b.sp.df <- ddply ( b.df, .(Species , Diel , Nest , Biogeo , BM , Hab , Diet, Pts.total, LS.total, N.tot), summarize, mean(N.tot))

# subset for species that were captured at at three distinct points
b.sp.df <- subset ( b.sp.df , b.sp.df$Pts.total >= 3 )

# subset for species that were captured in at least two different landscapes
b.sp.df <- subset ( b.sp.df , b.sp.df$LS.total >= 2 )

################################################ Plot N ~ FC for all each beetle species

# make plot matrix of individual beetle species N ~ FC

for (i in b.sp.df$Species) {
  
  # create df with entries as all captures of a particular species
  b.par.df <- subset (b.df , b.df$Species == i) # omitting/including zeros doesn't seem to affect model
  
  # plot N ~ FC for particular species
  plot(b.par.df$perFC_200, b.par.df$N, main = i)
  
  print (i)
  
}

################################################ CALCULATE N ~ FC RESPONSES FOR ALL BEETLE SPECIES

# create blank beetle FC response df
b.200respFC.df <- data.frame(
  Species = factor(), 
  Intercept = double(),
  Inter.error = double(),
  Inter.t.val = double(),
  Inter.p.val = double(), 
  b.slope_200 = double(),
  b.slope_200.error = double (),
  b.slope_200.t.val = double(),
  b.slope_200.p.val = double()
)

# remove species with insufficient values to run
b.sp.df = b.sp.df[!b.sp.df$Species == "Canthon.aff.semiopacus", ]
b.sp.df = b.sp.df[!b.sp.df$Species == "Deltochilum.dentipes", ]
b.sp.df = b.sp.df[!b.sp.df$Species == "Uroxys.sp1", ]
b.sp.df = b.sp.df[!b.sp.df$Species == "Dichotomius.depressicollis", ]
b.sp.df = b.sp.df[!b.sp.df$Species == "Scybalocanthon.nigriceps", ]

# # remove species with fewer than 2 captures at a point
# b.sp.df = b.sp.df[!b.sp.df$Species == "Dichotomius.sp5", ]
# b.sp.df = b.sp.df[!b.sp.df$Species == "Uroxys.sp4", ]
# b.sp.df = b.sp.df[!b.sp.df$Species == "Coprophanaeus.cerberus", ]



# check that it worked
b.sp.df$Species = factor(b.sp.df$Species) #weirdly important step. 
levels(b.sp.df$Species)


# set levels of species var

levels(b.200respFC.df$Species) <- unique (b.sp.df$Species)

# initialize row 1 as first slot of b.200respFC.df to be filled
row = 1

# iterate response calculation over all beetle species and store in b.200respFC.df
for (i in b.sp.df$Species) {
  
  # create df with entries as all captures of a particular (par) species
  b.par.df <- subset (b.df , b.df$Species == i) 
  
  # scale both longitude and forest cover to the mean
  b.par.df$Long.mscale = scale(b.par.df$Long, center = TRUE, scale = FALSE)
  b.par.df$perFC_200.mscale = scale(b.par.df$perFC_200, center = TRUE, scale = FALSE)
  
  # model N ~ FC for particular species
  b.glmm.par <- glmer (N ~ perFC_200.mscale  + Long.mscale + (1|Landscape), data = b.par.df, family = poisson)
  
  # add row to response df with model output for particular species
  b.200respFC.df[row, 1] <- (i)
  b.200respFC.df[row,2:5] <- summary (b.glmm.par)$coefficients[1,]
  b.200respFC.df[row,6:9] <- summary (b.glmm.par)$coefficients[2,]
  
  # increase row counter
  row <- row + 1
}

b.200respFC.df

# round p-values in response df
b.200respFC.df$b.slope_200.p.val <- round(b.200respFC.df$b.slope_200.p.val, 3)

# create histogram of responses
hist (b.200respFC.df$b.slope_200, main = "Distribution of Beetle Responses to Forest Cover",
      xlab = "Response Coefficient",
      ylab = "Frequency",
      breaks = 10) 
abline(v = 0, col = "red")


# subset for species with significant b.slope_200 coefficients
b.sigResp <- subset (b.200respFC.df, b.200respFC.df$b.slope_200.p.val < 0.05)

hist (b.sigResp$b.slope_200, breaks = 15)
hist (b.sigResp$b.slope_200.p.val, breaks = 15)


### Export summary of betas generated by this modeling approach

# subset b.200respFC.df for b.slope_200s and p-values
b.200respFC.df_betas <- subset(b.200respFC.df, select = c("Species", "b.slope_200", "b.slope_200.p.val"))

# rename columns to specify modeling approach number
colnames(b.200respFC.df_betas) = c("Species", "b_glmm_200_b.slope_200", "b_glmm_200_p.val")

#  save as RDS object
# saveRDS(b.200respFC.df_betas, "b.Mod1_betas.df")

################################################ MERGE N ~ FC RESPONSES WITH b.sp.df

b.sp.df <- merge(b.200respFC.df, b.sp.df)

# make row names of b.sp.df be species names
rownames(b.sp.df) <- b.sp.df$Species

################################################ MODEL RESPONSES ~ TRAITS | PHYLO

#### Import phylo tree object 
b.tre <- read.tree("b.sp.tre")
b.tre$node.label <- NULL

# generate comparative phylogenetic autocorrelation object
b.comp_200p <- comparative.data(b.tre, b.sp.df, Species, vcv=TRUE)

# build pgls model
b_glmm_200.pgls <- pgls(b.slope_200 ~ BM + Diel + Nest + BM*Diel, b.comp_200p)

# + Hab  --- generates error
# + Biogeo -- generates error


# generate all possible models and rank by AIC

b_pgls_200_mods.table <- dredge(b_glmm_200.pgls, beta = c("none", "sd", "partial.sd"), evaluate = TRUE,
                                 rank = "AICc", fixed = NULL)

print (b_pgls_200_mods.table)

summary(b_glmm_200.pgls)



plot(a.df$perFC_800, a.df$perFC_800)

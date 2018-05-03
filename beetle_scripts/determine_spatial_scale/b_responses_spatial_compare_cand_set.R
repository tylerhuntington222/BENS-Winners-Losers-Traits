################################################################################
# b_responses_spatial_compare_cand_set.R

# PURPOSE: To compare beetle sensitivity responses (Abundance ~ Forest Cover) 
#          at different buffer sizes around forest capture points using a 
#          (200m, 400m, 1000m, 3000m). "top model" approach.
#          Generates a sensitivity response model
#          for each beetle species at each spatial scale separately, then
#          creates a histogram of the species/spatial-scale combinations
#          resulted in a model with AIC < 2.

# INPUTS: clean beetle dataset outputted by `b_df_construct.R` to the 
#         `Intermediates/` directory.as a binary (.rds) file.

# OUTPUTS: Histogram plot outputted to `Figures/Exploratory/` directory.
#          Bin labels are spatial scales and bin counts are number of beetle
#          species for which modeling their abundance as a function of 
#          forest cover (controlling for landscape) resulted in a GLM
#          with AIC < 2.
            

# AUTHOR: Tyler Huntington

################################################################################

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

# Set WLT Analysis directory as working directory for this script.
# Filepaths are relative, so running the below code
# without modification on any machine should set the wd properly.
tryCatch (
  {
  this.dir <- dirname(parent.frame(2)$ofile)
  },
  error = function(x) {
    this.dir <<- dirname(rstudioapi::getActiveDocumentContext()$path)
  }
)
setwd(this.dir)
setwd("../../..")


########################### READ IN & ORGANIZE DATA ############################
# load pre-processed beetle data
b.df <- readRDS("Intermediates/b.df.rds")

# create df with unique entry for each beetle species and associated attributes
b.sp.df <- ddply ( b.df, .(Species , Diel , Nest , Biogeo , 
                           BM , Hab , Diet, N.tot, Pts.total, LS.total), 
                   summarize, 
                   mean(N.tot))

# subset for species captured more than three times
b.sp.df <- subset(b.sp.df, b.sp.df$Pts.total >= 3)

# subset for species that were captured in at least two different landscapes
b.sp.df <- subset ( b.sp.df , b.sp.df$LS.total >= 2 )

### Run N ~ FC multivariate regressions for each species using FC vals at each buffer size as IVs

# initialize blank df for top spatial models per species
a_spatial_mod.df <- data.frame("Species" = factor(), 
                               "Top_Spatial_Mod" = factor())

# initialize model counter vars
b1000_count_b <- 0
b800_count_b <- 0
b3000_count_b <- 0
b400_count_b <- 0
b200_count_b <- 0

for (i in b.sp.df$Species) {
  
  # initialize species row to be inserted in top model df
  a_species_row <- NULL
  
  # create df with entries as all captures of a particular species
  b.par.df <- subset (b.df , b.df$Species == i) 
  
  # regress N ~ FC for particular species at multiple spatial scales
  b.lm.par <- glmer (b.par.df$N ~ b.par.df$perFC_200
                     + b.par.df$perFC_400
                     + b.par.df$perFC_800
                     + b.par.df$perFC_1000
                     + b.par.df$perFC_3000
                     + (1|b.par.df$Landscape), family = poisson)
  
  options(na.action = "na.fail")
  b.par.spatial.mod.comp <- dredge(b.lm.par, beta = c("none", "sd", "partial.sd"), evaluate = TRUE,
                                   m.lim = c(1, 1), rank = "AICc", fixed = NULL)
  
  b.par.top.models <- get.models(b.par.spatial.mod.comp, subset = delta < 2)
  
  cand_set_codes <- rownames (summary(b.par.top.models))
  
  # KEY FOR MODEL NUMBERS IN DREDGE TABLE
  # Model 2 = 1000m buffer
  # Model 17 = 800m buffer
  # Model 5 = 3000m buffer
  # Model 9 = 400m buffer
  # Model 3 = 200m buffer
  
  
  # use model code to assign top model for particular species
  if ( "2" %in% cand_set_codes) {
    b1000_count_b <- b1000_count_b + 1
  }  
  if ("17" %in% cand_set_codes) {
    b800_count_b <- b800_count_b + 1
  } 
  if ("5" %in% cand_set_codes) {
    b3000_count_b <- b3000_count_b + 1
  } 
  if ("9" %in% cand_set_codes) {
    b400_count_b <- b400_count_b + 1
  } 
  if ("3" %in% cand_set_codes) {
    b200_count_b <- b200_count_b + 1
  }
  
}

##################### ANALYZE RESULTS OF SPATIAL SCALE COMPS ##################
# print counts as a table
print(data.frame(b200_count_b, b400_count_b, 
                 b800_count_b, b1000_count_b, b3000_count_b))

# plot distribution of top spatial models across beetle species
jpeg("Figures/Exploratory/b_spatial_scale_comp_hist_cand_set_approach.jpg")
par(mar = c(5, 5, 5, 5) + 0.1)
counts = c(b200_count_b, b400_count_b, b800_count_b, 
           b1000_count_b, b3000_count_b)
labels = c("200", "400", "800", "1000", "3000")
barplot <- barplot(counts, space = 0.25, names.arg = labels, ylim = c(0 , 25),
                   main = "Beetles: Candidate Set Approach\n for Determining Optimal Spatial Scale",
                   xlab = "Buffer Radius (m)",
                   ylab = "# Species w/ Model in Candidate Set")
text(barplot, counts, labels = counts, pos = 3)
dev.off()




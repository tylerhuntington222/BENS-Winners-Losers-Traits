################################################################################
# a_responses_spatial_compare_top_mod.R

# PURPOSE: To compare avian sensitivity responses (Abundance ~ Forest Cover) 
#          at different buffer sizes around forest capture points 
#          (200m, 400m, 1000m, 3000m) using a "candidate set" approach. 
#          Generates a sensitivity response model
#          for each avian species at each spatial scale separately, then
#          creates a histogram of spatial scales which best modeled 
#          (i.e. yielded the lowest AIC model) for each species. 

# INPUTS: clean avian dataset outputted by `b_df_construct.R` to the 
#         `Intermediates/` directory.as a binary (.rds) file.

# OUTPUTS: Histogram plot outputted to `Figures/Exploratory/` directory.
#          Bin labels are spatial scales and bin counts are number of avian
#          species for which the corresponding spatial scale yielded the
#          GLM of N~FC with the lowest AIC value (i.e. the "top model").
            

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
                           Fecund, DietPlas, Diet, Habitat, HabPlas, 
                           LS.total, Pts.total), 
                   summarize, N.tot = mean(N.tot))

# Run N ~ FC multivariate regressions for each species using FC vals 
# at each buffer size as IVs

# initialize blank df for top spatial models per species
a_spatial_mod.df <- data.frame("Species" = factor(), 
                               "Top_Spatial_Mod" = factor())

# initialize model counter vars
b1000_count_a_tm <- 0
b800_count_a_tm <- 0
b3000_count_a_tm <- 0
b400_count_a_tm <- 0
b200_count_a_tm <- 0

# iterate over species

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
  a.par.spatial.mod.comp <- dredge(a.lm.par, 
                                   beta = c("none", "sd", "partial.sd"), 
                                   evaluate = TRUE,
                         m.lim = c(1, 1), rank = "AICc", fixed = NULL)

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
    b1000_count_a_tm <- b1000_count_a_tm + 1
    
  } else if (top_model_code == "17") {
    top_mod <- 800
    b800_count_a_tm <- b800_count_a_tm + 1
    
  } else if (top_model_code == "5") {
    top_mod <- 3000
    b3000_count_a_tm <- b3000_count_a_tm + 1
    
  } else if (top_model_code == "9") {
    top_mod <- 400
    b400_count_a_tm <- b400_count_a_tm + 1
    
  } else if (top_model_code == "3") {
    top_mod <- 200
    b200_count_a_tm <- b200_count_a_tm + 1
  }
  
  a_species_row <- data.frame(factor(i), factor(top_mod))
  a_spatial_mod.df <- rbind (a_spatial_mod.df, a_species_row)
  print(a_spatial_mod.df)
}

# rectify column names of a_spatial_mod.df
colnames(a_spatial_mod.df) <- c("Species", "Top_Spatial_Mod")

############ ANALYZE RESULTS OF SPATIAL SCALE COMPS ###########################

# plot distribution of top spatial models across beetle species
jpeg("Figures/Exploratory/a_spatial_scale_comp_top_model_approach.jpeg")
par(mar = c(5, 5, 5, 5) + 0.1)
counts = c(b200_count_a_tm, b400_count_a_tm, 
           b800_count_a_tm, b1000_count_a_tm, b3000_count_a_tm)
labels = c("200", "400", "800", "1000", "3000")
barplot <- barplot(counts, space = 0.25, names.arg = labels, ylim = c(0 ,50),
                   main = "Birds: Top Model Approach\n for Determining Spatial Scale",
                   xlab = "Buffer Radius (m)",
                   ylab = "# Species w/ Model in Candidate Set")
text(barplot, counts, labels = counts, pos = 3)
dev.off()

table(a_spatial_mod.df$Top_Spatial_Mod)







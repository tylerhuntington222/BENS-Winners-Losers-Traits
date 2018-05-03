################################################################################
# b_responses_spatial_compare_top_mod.R

# PURPOSE: To compare beetle sensitivity responses (Abundance ~ Forest Cover) 
#          at different buffer sizes around forest capture points 
#          (200m, 400m, 1000m, 3000m) using a "candidate set" approach. 
#          Generates a sensitivity response model
#          for each beetle species at each spatial scale separately, then
#          creates a histogram of spatial scales which best modeled 
#          (i.e. yielded the lowest AIC model) for each species. 

# INPUTS: clean beetle dataset outputted by `b_df_construct.R` to the 
#         `Intermediates/` directory.as a binary (.rds) file.

# OUTPUTS: Histogram plot outputted to `Figures/Exploratory/` directory.
#          Bin labels are spatial scales and bin counts are number of beetle
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


############## READ IN & ORGANIZE DATA ########################################

# load base beetle dataframe
b.df <- readRDS("Intermediates/b.df.rds")

# create df with unique entry for each beetle species and associated attributes
b.sp.df <- ddply ( b.df, .(Species , Diel , Nest , Biogeo , BM , Hab , 
                           Diet, N.tot), summarize, mean(N.tot))

# subset relevant columns
b.sp.df <- b.sp.df[,1:8]

# make into 1 column df with only species names
b.sp.df <- data.frame(unique(b.df$Species))
colnames(b.sp.df) <- ("Species")


# Run N ~ FC multivariate regressions for each species 
#using FC vals at each buffer size as IVs

# initialize blank df for top spatial models per species
b_spatial_mod.df <- data.frame("Species" = factor(), 
                               "Top_Spatial_Mod" = factor())


# initialize model counter vars
b1000_count_b_tm <- 0
b800_count_b_tm <- 0
b3000_count_b_tm <- 0
b400_count_b_tm <- 0
b200_count_b_tm <- 0

for (i in b.sp.df$Species) {
  
  # initialize species row to be inserted in top model df
  b_species_row <- NULL
  
  # create df with entries as all captures of a particular species
  b.par.df <- subset (b.df , b.df$Species == i) 
  
  # regress N ~ FC for particular species at multiple spatial scales
  tryCatch({
    b.lm.par <- glmer (b.par.df$N ~ b.par.df$perFC_200
                       + b.par.df$perFC_400
                       + b.par.df$perFC_800
                       + b.par.df$perFC_1000
                       + b.par.df$perFC_3000
                       + (1|b.par.df$Landscape), family = poisson)
    options(na.action = "na.fail")
    b.par.spatial.mod.comp <- dredge(b.lm.par, 
                                     beta = c("none", "sd", "partial.sd"), 
                                     evaluate = TRUE,
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
    b_spatial_mod.df <- rbind (b_spatial_mod.df, b_species_row)
    print(b_spatial_mod.df)
    }, error = function(e) {
      print(i)
      print(e)
    }
  )
}

# rectify column names of b_spatial_mod.df
colnames(b_spatial_mod.df) <- c("Species", "Top_Spatial_Mod")

################ ANALYZE RESULTS OF SPATIAL SCALE COMPS #######################

# print table of results
table(b_spatial_mod.df$Top_Spatial_Mod)

# plot distribution of top spatial models across beetle species
jpeg("Figures/Exploratory/b_spatial_scale_comp_hist_top_model_approach.jpeg")
par(mar = c(4, 5, 3, 3) + 0.1)
counts = c(b200_count_b_tm, b400_count_b_tm, b800_count_b_tm, 
           b1000_count_b_tm, b3000_count_b_tm)
labels = c("200", "400", "800", "1000", "3000")
barplot <- barplot(counts, space = 0.25, names.arg = labels, ylim = c(0 , 15),
                   main = "Beetles: Top Model Approach\n for Determining Spatial Scale",
                   xlab = "Buffer Radius (m)",
                   ylab = "# Species w/ Spatial Scale as Top Model")
text(barplot, counts, labels = counts, pos = 3)
dev.off()








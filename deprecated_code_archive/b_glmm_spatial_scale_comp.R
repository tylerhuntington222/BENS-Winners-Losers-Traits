#####################################################################################################################################
################################## BEETLE N~FC REGRESSION MODELS COMPARISON #########################################################
#####################################################################################################################################

# PURPOSE: compare beetle responses (N.point ~ FC) at different buffer sizes around forest points (200m, 400m, 1000m, 3000m )

# AUTHORS: Tyler Huntington, Elizabeth Nichols, Larissa Boesing, Jean Paul Metzger


# set working directory
setwd("/Users/nicholslab/Dropbox/Interface/Dados/Non-Geospatial/Winners_Losers_Traits/Analysis/Data")

# load base beetle dataframe
b.df <- readRDS("b.df.rds")

head(b.df)

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

# create df with unique entry for each beetle species and associated attributes
b.sp.df <- ddply ( b.df, .(Species , Diel , Nest , Biogeo , BM , Hab , Diet, N.tot), summarize, mean(N.tot))

# subset relevant columns
b.sp.df <- b.sp.df[,1:8]

# make into 1 column df with only species names
b.sp.df <- data.frame(unique(b.df$Species))
colnames(b.sp.df) <- ("Species")

### add perFC columns with values from each buffer size 200m, 400m, 1000m, 3000m

# import buffer FC data
fc_pt_buffer.df <- read.csv("Landscape_composition_pitfalls.csv", header = T)
buffer_scales.vec <- c(200, 400, 1000, 3000)

# loop over spatial scales and add FC vals to b.df
for (i in buffer_scales.vec) {
  
  # initialize blank df for particular buffer size
  fc_par_buff.df = NULL
  
  # subset frag stats for buffer size
  fc_par_buff.df <- subset (fc_pt_buffer.df, fc_pt_buffer.df$Buffer_m == i)
  
  # subset for forest cover
  fc_par_buff.df <- subset (fc_par_buff.df, fc_par_buff.df$Class == 'forest')
  
  # collapse traps to point level, averaging trap buffer FC values
  fc_par_buff.df <- ddply (fc_par_buff.df, .(Landscape, Point, Code_min), 
                        summarize, perFC = mean(Per_area))
  
  colnames(fc_par_buff.df) <- c( "Landscape", "Point" , "TrapID" , paste("pt_", i ,"m_perFC", sep = ""))
  
  # subset for only TrapID and buffer FC cols
  fc_par_buff.df <- subset (fc_par_buff.df, select = c( "TrapID", (paste("pt_", i ,"m_perFC", sep = ""))))
  
  # merge particular buffer FC data with b.df
  b.df <- merge(b.df, fc_par_buff.df, by = "TrapID")
  
}

### Run N ~ FC multivariate regressions for each species using FC vals at each buffer size as IVs

# initialize blank df for top spatial models per species
b_spatial_mod.df <- data.frame("Species" = factor(), 
                               "Top_Spatial_Mod" = factor())

for (i in b.sp.df$Species) {
  
  # initialize species row to be inserted in top model df
  b_species_row <- NULL
  
  # create df with entries as all captures of a particular species
  b.par.df <- subset (b.df , b.df$Species == i) 
  
  # regress N ~ FC for particular species at multiple spatial scales
  b.glmm.par <- glmer (b.par.df$N ~ b.par.df$pt_200m_perFC
                     + b.par.df$pt_400m_perFC
                     + b.par.df$pt_1000m_perFC
                     + b.par.df$pt_3000m_perFC
                     + (1|b.par.df$Landscape), family = poisson)
  
  
  
  options(na.action = "na.fail")
  b.par.spatial.mod.comp <- dredge(b.lm.par, beta = c("none", "sd", "partial.sd"), evaluate = TRUE,
                         m.lim = c(1,1), rank = "AICc", fixed = NULL)

  b.par.top.model <- get.models(b.par.spatial.mod.comp, subset = delta == 0)
  
  top_model_code <- rownames (summary(b.par.top.model))
  
  # KEY FOR MODELS IN DREDGE TABLE
  # Model 2 = 1000m buffer
  # Model 3 = 200m buffer
  # Model 5 = 3000m buffer
  # Model 9 = 400m buffer
  
  # use model code to assign top model for particular species
  if (top_model_code == "2") {
    top_mod <- "1000"
  } else if (top_model_code == "3") {
    top_mod <- "200"
  } else if (top_model_code == "5") {
    top_mod <- "3000"
  } else if (top_model_code == "9") {
    top_mod <- "400"
  }
  
  b_species_row <- data.frame(factor(i), factor(top_mod))
  b_spatial_mod.df <- rbind (b_spatial_mod.df, b_species_row)

}

b_spatial_mod.df$Top_Spatial_Mod <- factor(b_spatial_mod.df)

# rectify column names of b_spatial_mod.df
colnames(b_spatial_mod.df) <- c("Species", "Top_Spatial_Mod")

################################################ ANALYZE RESULTS OF SPATIAL SCALE COMPARISON

# plot distribution of top spatial models across beetle species
plot(b_spatial_mod.df$Top_Spatial_Mod, ylab = "∆AIC Top Model Counts (by species)", 
     xlab = "Point-Level Buffer Radius for Calculating FC (m)")

table(b_spatial_mod.df$Top_Spatial_Mod)





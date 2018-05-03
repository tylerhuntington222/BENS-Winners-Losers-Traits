#####################################################################################################################################
################################## BEETLE N~FC REGRESSION MODELS COMPARISON #########################################################
#####################################################################################################################################

# PURPOSE: compare the three different modeling approaches used to generate betas (i.e. abundance responses to declining FC)

# MODEL SUMMARIES
# Model 1: N (Point) ~ FC (200m_pt_buffer) 
# Model 2: N (Point) ~ FC (Landscape)
# Model 3: N (Landscape) ~ FC (Landscape)


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

###################################

###### READ IN DATA

b.Mod1_betas.df <- readRDS("b.Mod1_betas.df")

b.Mod2_betas.df <- readRDS("b.Mod2_betas.df")

b.Mod3_betas.df <- readRDS("b.Mod3_betas.df")


###### MERGE MODEL OUTPUTS INTO SINGLE DF

b_mod_comp.df <- merge(b.Mod1_betas.df, b.Mod2_betas.df)
b_mod_comp.df <- merge(b_mod_comp.df, b.Mod3_betas.df)


###### SYNTHESIZE ALL THREE MODEL OUTPUTS

# calculate mean slope coefficient for each species and append column to b_mod_comp.df
b_mod_comp.df$avg_Slope <- rowMeans(b_mod_comp.df[,c(2,4,6)])
                                   
# calculate mean slope p-value for each species and append column to b_mod_comp.df
b_mod_comp.df$avg_pVal <- round(rowMeans(b_mod_comp.df[,c(3,5,7)]), 3)

b_mod_comp.df[order(b_mod_comp.df$avg_pVal),]


###### MAKE COMPARATIVE PLOTS

# make histograms of slopes for all three models

hist(b_mod_comp.df$Mod1_Slope, main = "Distribution of Model 1 Slope Coefficients")
abline(v=0, col = "red") 

hist(b_mod_comp.df$Mod2_Slope, main = "Distribution of Model 2 Slope Coefficients")
abline(v=0, col = "red")

hist(b_mod_comp.df$Mod3_Slope, main = "Distribution of Model 3 Slope Coefficients")
abline(v=0, col = "red")

# make histograms of p-values

hist(b_mod_comp.df$Mod1_p.val)
hist(b_mod_comp.df$Mod2_p.val)
hist(b_mod_comp.df$Mod3_p.val)



# make boxplots of slopes from all three models
boxplot(b_mod_comp.df$Mod1_Slope, b_mod_comp.df$Mod2_Slope, b_mod_comp.df$Mod3_Slope,
        main = "Magnitude of slope coefficients", ylab = "Slope Coefficient", 
        names = c("Mod1" , "Mod2", "Mod3"),
        ylim = c(-0.25, 0.15))
        abline(h=0, col = "red")

# make boxplots of p-vals from all three models
boxplot(b_mod_comp.df$Mod1_p.val, b_mod_comp.df$Mod2_p.val, b_mod_comp.df$Mod3_p.val,
        main = "Significance levels of slopes", ylab = "p-value", 
        names = c("Mod1" , "Mod2", "Mod3"))

# determine proportion of p-values that are significant from each modeling approach
prop_sig1 = (length(b_mod_comp.df$Mod1_p.val[b_mod_comp.df$Mod1_p.val<0.05]))/length(b_mod_comp.df$Mod1_p.val)
prop_sig2 = (length(b_mod_comp.df$Mod2_p.val[b_mod_comp.df$Mod2_p.val<0.05]))/length(b_mod_comp.df$Mod2_p.val)
prop_sig3 = (length(b_mod_comp.df$Mod3_p.val[b_mod_comp.df$Mod3_p.val<0.05]))/length(b_mod_comp.df$Mod3_p.val)
prop_mods_sig = c(round(prop_sig1, 2) , round (prop_sig2, 2), round (prop_sig3, 2))

# create barplot of proportion of p-values < 0.05 for each model
plot.new()
orig.par <- par()
orig.par$cex.lab


par (orig.par)

mod_labels <- c("(1) N(pt)~FC(200m)", "(2) N(pt)~FC(LS)", "(3) N(LS) ~ FC(LS)")
barplot <- barplot(prop_mods_sig, 13, axes = TRUE, axisnames = TRUE, names.arg = mod_labels,
                   main = "Comparison of significant betas yielded by each model", 
                   ylab = "Proportion of Slope p-values < 0.05",
                   xlab = "Model",
                   ylim = c(0,1),
                   cex.sub = 0.8)

text(barplot, prop_mods_sig, labels = prop_mods_sig, pos = 3)
 









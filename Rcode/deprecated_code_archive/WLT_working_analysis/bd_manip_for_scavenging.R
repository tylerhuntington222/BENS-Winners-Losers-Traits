#########################

# Authors: Nichols, Alarcon, Metzger, Langhans. DRAFT-DO NOT CIRCULATE
# Goal: evaluate dung beetle diversity thresholds over forest gradient
# 1. Clean data
# 2. Join with cleaned data to calculate TD, FD, PD per trap
# 3. Evaluate trends over gradient of three metrics

#########################

#########################
# Set directory         
#########################

# setwd("~/Dropbox/Grants/Active/Interface/Dados/Non-Geospatial/Biodiversity/Analysis/Data") #liz
setwd("~/Dropbox/Interface/Dados/Non-Geospatial/Biodiversity/Analysis/Data") #Kelley

# Specify directory for results and figures
# results = "~/Dropbox/Grants/Active/Interface/Dados/Non-Geospatial/Biodiversity/Analysis/Data/results" #liz
# figures = ""~/Dropbox/Grants/Active/Interface/Dados/Non-Geospatial/Biodiversity/Analysis/Data/figures" #liz

results = "~/Dropbox/Interface/Dados/Non-Geospatial/Biodiversity/Analysis/Data/results" #Kelley
figures = "~/Dropbox/Interface/Dados/Non-Geospatial/Biodiversity/Analysis/Data/figures" #Kelley

#########################
# Load libraries   
#########################

library(rtf)
library(DataCombine)
library(R2admb)
library(glmmADMB)
library(MuMIn)
library(lme4)
library(FD)
library(vegan)
library(nlme)
library(picante)
library(stringr)
library(MASS)
library("lme4")
library("coda")
library("reshape")
require("lattice")
library("plotrix")
library("dataframes2xls")
library("gdata")
library("vegan")
library("fields")
library("Matrix")  
library("abind")  
library("doBy")
library("gdistance")
library("sp")
library("xlsx")
library("reshape2")
library(languageR) 
library(car)  
library(ape)  
library(geiger)
library(phytools)
library(ade4)
library(betareg)
library(stats)
library(lmerTest)
library(perm)
library(graphics)
library(plotrix)
library(MASS)  
library(arm)  
library(xtable)  
library(fields)
library(Hmisc) #interferes with plyr, is it necessary for Lari code?
library(abind)  
library(gtools)
library(gmodels) 
library(geiger)
library(phangorn)
library(gridExtra)
library(gridGraphics)
library(ggplot2)
library(plyr)

options(stringsAsFactors=FALSE)

#############################
# Some basic functions    
#############################

SE = function (x)  sd(x)/sqrt(length(x))
Uniq = function(x) length(unique(x)) 
N <- function (x) sum(x,na.rm = TRUE)
S = function (x) length (x[x>0])
Sd = function (x) { p<- x/sum(x)
1-sum(p^2)}  

#function that moves item in list as you specify
moveme <- function (invec, movecommand) {
  movecommand <- lapply(strsplit(strsplit(movecommand, ";")[[1]], 
                                 ",|\\s+"), function(x) x[x != ""])
  movelist <- lapply(movecommand, function(x) {
    Where <- x[which(x %in% c("before", "after", "first", 
                              "last")):length(x)]
    ToMove <- setdiff(x, Where)
    list(ToMove, Where)
  })
  myVec <- invec
  for (i in seq_along(movelist)) {
    temp <- setdiff(myVec, movelist[[i]][[1]])
    A <- movelist[[i]][[2]][1]
    if (A %in% c("before", "after")) {
      ba <- movelist[[i]][[2]][2]
      if (A == "before") {
        after <- match(ba, temp) - 1
      }
      else if (A == "after") {
        after <- match(ba, temp)
      }
    }
    else if (A == "first") {
      after <- 0
    }
    else if (A == "last") {
      after <- length(myVec)
    }
    myVec <- append(temp, values = movelist[[i]][[1]], after = after)
  }
  myVec
}

# this function can add new branches of specified length at specified nodes
bind.tip<-function(tree,tip.label,edge.length=NULL,where=NULL){
  if(is.null(where)) where<-length(tree$tip)+1
  tip<-list(edge=matrix(c(2,1),1,2),
            tip.label=tip.label,
            edge.length=edge.length,
            Nnode=1)
  class(tip)<-"phylo"
  obj<-bind.tree(tree,tip,where=where)
  return(obj)
}

# source: httFS://github.com/willpearse/phyloGenerator/blob/master/merging.R, 
# takes a phylogeny with genera and adds species within those genera as polytomies.
# Edited to get rid of a default step which cut genera branch lenghts in half, and then gave species branch lengths of 1/2 genera branch length
# now it maintains genera branch lengths, and gives species branches lengths of 1.


#Replace a given list of genera in a phylogeny with polytomes of their species
#Needs a tree, and vectors of genera and species
#NOTE: vectors should be of same length, e.g. c('Quercus', 'Quercus', 'Phocoena') and c('ilex', 'robur', 'phocoena')
make.composite.with.polytomies.2 <- function(tree, genera, species, max.genus.age=NA){
  #Functions#
  #Binds a clade into a phylogeny, replacing a given tip, and making it ultrametric	
  #A hack to fix the ultrametric nature, but APE is a pile of shite so it's all I can do	
  bind.ultrametric.by.replacement <- function(backbone, donor, replacing.tip.label, donor.length=NA){	
    #Find the species we're replacing	
    bind.point <- which(backbone$tip.label == replacing.tip.label)	
    #Bind *badly* the two trees together	
    backbone <- bind.tree(backbone, donor, where=bind.point)	
    #Now we've bound, where's the genus located?	
    which.tip <- which(backbone$tip.label == donor$tip.label[1])	
    #What node is above the genus?	
    which.node <- backbone$edge[which(backbone$edge[,2] == which.tip),1]	
    #Where is that node in the edge.list (and thus edge.length list)	
    which.edge <- which(backbone$edge[,2] == which.node)	
    #What length is that branch?	
    tip.length <- backbone$edge.length[which.edge]
    if(is.na(donor.length)){
      #DON'T REPLACE WITH TIP.LENGTH/2
      backbone$edge.length[which.edge] <- tip.length
    } else {
      #Correct for the manual adjustment...
      backbone$edge.length[which.edge] <- tip.length - donor.length/2
    }
    #Return the ultrametric tree!	
    return(backbone)	
  }	
  #Make a polytomy from a given species list, optionally with a tip.length	
  make.polytomy <- function(species, tip.length=NA){	
    d.f <- data.frame(spp=factor(species))	
    polytomy <- as.phylo.formula(~spp, data=d.f)	
    if(!is.na(tip.length)) polytomy$edge.length <- rep(tip.length, length(species))	
    return(polytomy)	
  }	
  #Find the unique branch length of a tip on a given phlogeny	
  find.unique.branch.length <- function(tree, tip){	
    #What number of tip.label is the tip	
    which.tip <- which(tree$tip.label == tip)	
    #What edge is the tip.label on	
    which.edge <- which(tree$edge[,2] == which.tip)	
    #What's the edge.length for that edge?	
    tip.length <- tree$edge.length[which.edge]	
    #Return that edge.length	
    return(tip.length)	
  }	
  
  #Code#
  #The genera and species vectors *must* be characters, so convert
  genera<-as.character(genera);species<-as.character(species)
  #Go along all the genera to be replaced with polytomes
  for(genus in unique(genera)){
    #Make a vector of the species that will go into this particular polytomy
    species.to.bind <- species[genera == genus]
    #If the length of that vector is 1...
    if(length(species.to.bind) == 1){
      #...we don't need to bind anything in, so just change the tip's name
      tree$tip.label[tree$tip.label == genus] <- species.to.bind
    } else {
      #Other, find the branch length leading to the particular genus
      tip.length <- find.unique.branch.length(tree, genus)
      #Don't warn about edge issues unless we need to think about it...
      edge.warning <- NA
      #Set the genus to the correct age
      if(!is.na(max.genus.age)){
        if(max.genus.age*2 < tip.length){
          tip.length <- min(tip.length, max.genus.age*2)
          edge.warning <- tip.length
        }
      }
      #Make the polytomy of species, with branch lengths appropriate for that genus
      polytomy <- make.polytomy(species.to.bind, tip.length = 1)
      #Bind that polytomy in, ultrametrially, to that phylogeny
      tree <- bind.ultrametric.by.replacement(tree, polytomy, genus, edge.warning)
    }
  }
  #Return the tree
  return(tree)
}

# Funtion to add text to the corner of a plot
Corner_text <- function(text, location="topright"){
  legend(location,legend=text, bty ="n", pch=NA, cex = 3) 
}

# Banks-Leite 2014 function, calculate preference for FC given input presence/absence species/site matrix, and list of % forest cover for each site in order
weight.avg <- function(community, pref){
  weights <- numeric(ncol(community))
  for(i in seq(ncol(community))){
    x <- community[,i] * pref
    weights[i] <- mean(x[x!=0])
  }
  return(weights)
}

# Function to convert a base graphic to a grob
grab_grob <- function(){
  grid.echo()
  grid.grab()
}

#################
# Set data
#################

# taxonomic diversity
df = read.xls("Biodiversity_27_09_16.xlsx", sheet = 2, header = TRUE)

# functional diversity
# Abundance data comes from df
traits.df = read.xls("Species_traits_23_09_16.xlsx", sheet = 1, header = TRUE, stringsAsFactors = FALSE) 

# phylogenetic diversity
# Trait data, used to get list of our genera of interest comes from traits.df
# Phylogenetic tree supplied by Monaghan
base.phy <- read.nexus("Monaghan_Scarabaeinae_r8s.tre_02.txt")
# Data table from Monaghan to match code in tree to genus, species
code.df = read.xlsx("Monaghan_dungbeetles.mrb.xlsx", sheetIndex = 1, header = TRUE)

# buffer data
buffer.df <- read.csv("Landscape_composition_pitfalls.csv", header = TRUE)

#landscape configuration data
land.config <- read.csv("Fragstats.csv", header = T)

# set seed--so all random number generations will be reproducible
set.seed(200) 


######################################### Manip Taxonomic ############################################################

##########################
# convert to long format
##########################

df_melt = melt(df, id.vars=1:16, measure.vars=17:78, variable_name="code")
names(df_melt)[match("value", names(df_melt))]="N.trap" 

###############################################
# get rid of observations that are 0 bc missing
###############################################

df_melt1 <- subset(df_melt, df_melt$missing==0)

#########################################
# get rid of outliers 
########################################

# 329-0-M--23% FC
df_melt2 <- subset(df_melt1, !(df_melt1$Landscape==329 & df_melt1$Point==0))

##########################################
# Separate based on habitat use classification
#############################################

# set up df with just our species and habitat use
hab.assoc.1 <- subset(traits.df, traits.df$Reference=="Nichols_unpub_2016")
hab.assoc.2 <- subset(hab.assoc.1, hab.assoc.1$Classification=="FS" | hab.assoc.1$Classification=="NFS")
hab.assoc.3 <- NULL
hab.assoc.3$code <- hab.assoc.2$Code
hab.assoc.3$Classification <- hab.assoc.2$Classification
hab.assoc <- as.data.frame(hab.assoc.3)

# set up df
# all = false here gets rid of morpho species
hab_exp1 = merge(df_melt2, hab.assoc, by = c("code"), all = FALSE)

# subset for only forest captures
hab_exp2 <- subset(hab_exp1, hab_exp1$Habitat=="M")

#separate by habitat use
df.group.1 = subset(hab_exp2, hab_exp2$Classification=="FS")  ; dim(df.group.1) ; df.group.1$code = factor(df.group.1$code)
df.group.2 = subset(hab_exp2, hab_exp2$Classification=="NFS")  ; dim(df.group.2) ; df.group.2$code = factor(df.group.2$code)

# summarize mean richness and abundance per point per landscape, per habitat use subset
# for each species within each point, sums up N over all traps, then creates a presence/absence column per species per point

#FS
sum.df.group.1 = ddply(df.group.1, .(Landscape, Habitat, Point, code), summarize, 
                       N = sum(N.trap))
sum.df.group.1$presence = ifelse(sum.df.group.1$N>=1, 1, 0)
#sums up presence/absence for all species/point [s] and sums up N for all species per point [N]
sum.df.group.1.2 = ddply(sum.df.group.1, .(Landscape, Habitat, Point), summarize, 
                         N = sum(N), 
                         S = sum(presence))
# get rid of points with 0 species
sum.df.group.1.3 <- subset(sum.df.group.1.2, sum.df.group.1.2$S > 0)
#rename for easy use
TD.FS <- sum.df.group.1.3
names(TD.FS)[names(TD.FS)=="S"] <- "TD"

# NFS
sum.df.group.2 = ddply(df.group.2, .(Landscape, Habitat, Point, code), summarize, 
                       N = sum(N.trap))
sum.df.group.2$presence = ifelse(sum.df.group.2$N>=1, 1, 0)
sum.df.group.2.2 = ddply(sum.df.group.2, .(Landscape, Habitat, Point), summarize, 
                         N = sum(N), 
                         S = sum(presence))
#get rid of points with 0 species
sum.df.group.2.3 <- subset(sum.df.group.2.2, sum.df.group.2.2$S > 0)
#rename for easy use
TD.NFS <- sum.df.group.2.3
names(TD.NFS)[names(TD.NFS)=="S"] <- "TD"

# Total
# subset so only contains FS and NFS species
hab_exp3 <- subset(hab_exp2, !hab_exp2$Classification=="na")
sum.df = ddply(hab_exp3, .(Landscape, Habitat, Point, code), summarize, 
               N = sum(N.trap)) 
sum.df$presence = ifelse(sum.df$N>=1, 1, 0)
sum.df2 = ddply(sum.df, .(Landscape, Habitat, Point), summarize, 
                N = sum(N), 
                S = sum(presence))
# get rid of points with 0 species
sum.df3 <- subset(sum.df2, sum.df2$S > 0)
#rename for easy use
TD.tot <- sum.df3
names(TD.tot)[names(TD.tot)=="S"] <- "TD"

# extract a list of all species analyzed--FS and NFS with non-zero captures
spec.keep.sum <- ddply(sum.df, .(code), summarize, N = sum(N))
spec.keep.sum.1 <- subset(spec.keep.sum, spec.keep.sum$N > 0 ) # this is a list with all species analyzed and how many of each were captured
spec.keep <- spec.keep.sum.1$code

############################################################## Manip Functional ####################################################
# 
# ###################################
# # Create abundance data input for FD
# ####################################
# 
# # work from dataframe with 0s and outliers removed
# abund.df_melt1 <- df_melt2
# 
# # summarize number of individuals of each species trapped per point
# abund.sum.df = ddply(abund.df_melt1, .(Landscape, Habitat, Point, code), summarize, N = sum(N.trap))
# 
# # subset for only Forest captures
# abund.sum.df2 <- subset(abund.sum.df, abund.sum.df$Habitat=="M")
# 
# # create a single ID column
# abund.sum.df2$ID <- paste(abund.sum.df2$Landscape, abund.sum.df2$Point, abund.sum.df2$Habitat, sep="-")
# # This will print an ID column that lists Landscape-Point-Habitat
# 
# # get rid of old ID columns
# abund.sum.df2$Landscape <- NULL
# abund.sum.df2$Point <- NULL
# abund.sum.df2$Habitat <- NULL
# 
# #get rid of points where 0 species were collected
# #create presence/absence column
# abund.sum.df2$presence = ifelse(abund.sum.df2$N>=1, 1, 0)
# #sum per point to get species richness per point
# richness.df = ddply(abund.sum.df2, .(ID), summarize, S = sum(presence))
# #create df with points to keep
# richness.df2 <- subset(richness.df, richness.df$S > 0)
# #merge with original dataframe to only preserve species with enough species
# abund.sum.df3 = merge(abund.sum.df2,richness.df2,  by = c("ID"), all = FALSE)
# #remove extra columns
# abund.sum.df3$presence <- NULL
# abund.sum.df3$S <- NULL
# #check that there are fewer points
# length(unique(abund.sum.df3$ID))
# length(unique(abund.sum.df2$ID))
# 
# # convert data back into wide format
# abund.final = dcast(abund.sum.df3, ID ~ code, value.var = "N")
# 
# # try to get species names alphabetical
# #column names alphabetical
# abund.final.alph <- abund.final[,order(colnames(abund.final))]
# #move ID column to front
# abun.temp <- abund.final.alph[moveme(names(abund.final.alph), "ID first")]
# 
# #########################################
# # Create separate PS, FS abundance inputs
# #########################################
# 
# abund.sum.df.habuse <- merge(abund.sum.df2, hab.assoc, by = c("code"), all = FALSE)
# abund.sum.df.NFS <- subset(abund.sum.df.habuse, abund.sum.df.habuse$Classification=="NFS")
# abund.sum.df.FS <- subset(abund.sum.df.habuse, abund.sum.df.habuse$Classification=="FS")
# 
# #sum per point to get species richness per point
# richness.NFS = ddply(abund.sum.df.NFS, .(ID), summarize, S = sum(presence))
# # get rid of points with 0  species
# #create df with points to keep
# richness.NFS2 <- subset(richness.NFS, richness.NFS$S > 0)
# #merge with original dataframe to only preserve points with enough species
# abund.sum.NFS2 = merge(abund.sum.df.NFS,richness.NFS2,  by = c("ID"), all = FALSE)
# #remove extra columns
# abund.sum.NFS2$presence <- NULL
# abund.sum.NFS2$S <- NULL
# abund.sum.NFS2$Classification <- NULL
# 
# #sum per point to get species richness per point
# richness.FS = ddply(abund.sum.df.FS, .(ID), summarize, S = sum(presence))
# # get rid of points with 0  species
# #create df with points to keep
# richness.FS2 <- subset(richness.FS, richness.FS$S > 0)
# #merge with original dataframe to only preserve species with enough species
# abund.sum.FS2 = merge(abund.sum.df.FS,richness.FS2,  by = c("ID"), all = FALSE)
# #remove extra columns
# abund.sum.FS2$presence <- NULL
# abund.sum.FS2$S <- NULL
# abund.sum.FS2$Classification <- NULL
# 
# length(unique(abund.sum.NFS2$ID)) #87
# length(unique(abund.sum.FS2$ID)) #84
# # compared to number of forest points included for pooled species == 93
# 
# # convert data back into wide format
# abund.final.NFS = dcast(abund.sum.NFS2, ID ~ code, value.var = "N")
# # try to get species names alphabetical
# #column names alphabetical
# abund.final.NFS.alph <- abund.final.NFS[,order(colnames(abund.final.NFS))]
# #move ID column to front
# abun.temp.NFS <- abund.final.NFS.alph[moveme(names(abund.final.NFS.alph), "ID first")]
# #make ID rowname
# rownames(abun.temp.NFS) <- abun.temp.NFS$ID
# abun.temp.NFS$ID <- NULL
# 
# # convert data back into wide format
# abund.final.FS = dcast(abund.sum.FS2, ID ~ code, value.var = "N")
# # try to get species names alphabetical
# #column names alphabetical
# abund.final.FS.alph <- abund.final.FS[,order(colnames(abund.final.FS))]
# #move ID column to front
# abun.temp.FS <- abund.final.FS.alph[moveme(names(abund.final.FS.alph), "ID first")]
# #make ID rowname
# rownames(abun.temp.FS) <- abun.temp.FS$ID
# abun.temp.FS$ID <- NULL
# 
# ##########################
# # Create trait input for FD
# ##########################
# 
# # Subset for only Nichols data
# traits.df1 <- subset(traits.df, traits.df$Reference=="Nichols_unpub_2016")
# 
# # Get rid of undesired traits and columns
# traits.df1$Biogeo_assoc <- NULL
# traits.df1$ID <- NULL
# traits.df1$Reference <- NULL
# traits.df1$Genus <- NULL
# traits.df1$spp <- NULL
# traits.df1$latin_bin <- NULL
# traits.df1$Bod_leng_mm <- NULL
# traits.df1$Bod_width_mm<- NULL
# traits.df1$Bod_leng_source <-NULL
# traits.df1$Diet <- NULL
# traits.df1$Colour <- NULL
# traits.df1$BM1 <- NULL
# traits.df1$BM2 <- NULL
# traits.df1$BM3 <- NULL
# traits.df1$BM4 <- NULL
# traits.df1$BM5 <- NULL
# traits.df1$Hab_pref <-NULL
# traits.df1$Classification <-NULL
# traits.df1$Lat <- NULL
# traits.df1$Long <- NULL
# traits.df1$Sources <- NULL
# traits.df1$Biogeo_prov_source <- NULL
# traits.df1$Diel_source <- NULL
# traits.df1$Diet_source <- NULL
# traits.df1$BM_source <- NULL
# traits.df1$Habitat_spec_source <- NULL
# traits.df1$Nesting_source <- NULL
# 
# # give it a shorter name
# trait.temp <- traits.df1
# 
# ####################################
# # Clean up abundance and trait input
# ####################################
# 
# #make species names column into row names in trait.temp, and same thing for ID in abun.temp.pre
# rownames(trait.temp) <- trait.temp$Code
# rownames(abun.temp) <- abun.temp$ID
# 
# trait.temp$Code <- NULL
# abun.temp$ID <- NULL
# 
# # match species names in both files
# colnames(abun.temp) = rownames(trait.temp)
# identical(rownames(trait.temp),colnames(abun.temp))
# 
# ##################################
# # Making a dendrogram of FD
# ##################################
# 
# # Took code for hclust and making a dendrogram from Cadotte 2009, which purports to use the same method as petchy 2002.
# # dist.ktab code is my own, recommended by Edwards 2013, as it allows input of multiple variable types.
# 
# # 1. Create distance matrix
# 
# # make df. into ktab file (a collection of several dataframes)
# # first, separate my dataframe of traits into 1 dataframe per trait
# kdf1 <- trait.temp[,-(2:3), drop = FALSE]
# kdf2 <- trait.temp[,-(1:2), drop = FALSE]
# kdf3.temp <- trait.temp[,-(1), drop = FALSE]
# kdf3 <- kdf3.temp[,-(2), drop = FALSE]
# # turn into list
# my.list <- list(kdf1, kdf2, kdf3)
# # turn into ktab
# kta <- ktab.list.df(my.list)
# 
# # caculate distance matrix
# # this method uses Gower distance NOT euclidean distance as used in Petchey 2002, Cadotte 2009
# # order of traits in ktab is: Diel, BM, Nest
# # create a vector specifying variable type in order, n = nominal, q = quantitative
# type <- c("N", "Q", "N")
# # calculate using dist.ktab, scale by range
# ktabdist <- dist.ktab(kta, type = type, option = "scaledBYrange")
# 
# # 2. Do hierarchical clustering as per Petchey & Gaston
# #"average" clustering method is UPGMA, used by Petchey, Edwards, and Cadotte
# 
# trait.clust<-hclust(ktabdist, method="average")
# 
# # 3. Transform trait.clust into phylo object
# clust.phy <- as.phylo(trait.clust)
# 
# # only keep species that are NFS and FS captured in forest.
# 
# # extract a list of NFS species
# NFS.spec <- subset(hab.assoc, hab.assoc$Classification=="NFS")
# NFS.spec.n <- NFS.spec$code #21
# 
# # create a list of NFS species that were captured more than 0 times in forest
# NFS.spec.sum <- ddply(abund.sum.NFS2, .(code), summarize, N.trap.tot = sum(N))
# NFS.spec.sum.0 <- subset(NFS.spec.sum, NFS.spec.sum$N.trap.tot==0)$code
# NFS.spec.n.2 <- subset(NFS.spec.n, !(NFS.spec.n %in% NFS.spec.sum.0))
# 
# # extract a list of FS species
# FS.spec <- subset(hab.assoc, hab.assoc$Classification=="FS")
# FS.spec.n <- FS.spec$code #17
# 
# # create a list of FS species that were captured more than 0 times in forest
# FS.spec.sum <- ddply(abund.sum.FS2, .(code), summarize, N.trap.tot = sum(N))
# FS.spec.sum.0 <- subset(FS.spec.sum, FS.spec.sum$N.trap.tot==0)$code
# FS.spec.n.2 <- subset(FS.spec.n, !(FS.spec.n %in% FS.spec.sum.0))
# 
# # trim all tips from dendrogram not in those 2 lists
# tiplabs <- clust.phy$tip.label
# tip.remove.1 <- subset(tiplabs, !(tiplabs %in% NFS.spec.n.2))
# tip.remove.2 <- subset(tip.remove.1, !(tip.remove.1 %in% FS.spec.n.2))
# clust.phy <- drop.tip(clust.phy, tip.remove.2, trim.internal = TRUE)
# 
# # save a PDF
# filename <- paste(figures, "Nichols_2016_dendrogram.pdf", sep = "/")
# pdf(file = filename, width = 10, height = 15)
# 
# par(2,2,2,5)
# 
# plot.phylo(clust.phy, show.tip.label = F, x.lim = c(0,0.7))
# tiplabels(text = c("Canthon amabilis", "Canthon aff. angularis", "Canthon ibarragrassoi", "Canthon aff. luctuosus",
#             "Canthon podagricus", "Canthon septemmaculatus", "Canthon virens chalybaeus",
#             "Chalcocopris hesperus", "Coprophaneus bellicosus", "Coprophanaeus cerberus", "Coprophanaeus saphirinus",
#             "Deltochilum brasilense", "Deltochilum dentipes", "Deltochilum furcatum", "Deltochilum morbillosum",
#             "Deltochilum rubripenne", "Dichotomius assifer", "Dichotomius bos", "Dichotomius aff. carbonarius",
#             "Dichotomius depressicollis", "Dichotomius aff. fissus", "Dichotomius mormon", "Dichotomius aff. carbonarius sp.2",
#             "Dichotomius quadrinodosus", "Eurysternus cyanescens", "Eurysternus francinae",
#             "Eurysternus hirtellus", "Eurysternus inflexus", "Eurysternus parallelus", "Ontherus sulcator",
#             "Phanaeus dejeani", "Phanaeus splendidulus", "Pseudocanthon aff. felix", "Scybalocanthon nigriceps", "Sylvicanthon aff. foveiventris"),
#           tip = c(1:35), adj = -.015, family = "serif", cex = 1.5, frame = "n", font = 3)
# 
# dev.off()
# 
# # save a tiff
# filename <- paste(figures, "Nichols_2016_dendrogram.tiff", sep = "/")
# tiff(filename = filename, width = 2550, height = 3000, compression = c("none"), res = 400)
# 
# plot.phylo(clust.phy, show.tip.label = F, x.lim = c(0,0.7))
# tiplabels(text = c("Canthon amabilis", "Canthon aff. angularis", "Canthon ibarragrassoi", "Canthon aff. luctuosus",
#             "Canthon podagricus", "Canthon septemmaculatus", "Canthon virens chalybaeus",
#             "Chalcocopris hesperus", "Coprophaneus bellicosus", "Coprophanaeus cerberus", "Coprophanaeus saphirinus",
#             "Deltochilum brasilense", "Deltochilum dentipes", "Deltochilum furcatum", "Deltochilum morbillosum",
#             "Deltochilum rubripenne", "Dichotomius assifer", "Dichotomius bos", "Dichotomius aff. carbonarius",
#             "Dichotomius depressicollis", "Dichotomius aff. fissus", "Dichotomius mormon", "Dichotomius aff. carbonarius sp.2",
#             "Dichotomius quadrinodosus", "Eurysternus cyanescens", "Eurysternus francinae",
#             "Eurysternus hirtellus", "Eurysternus inflexus", "Eurysternus parallelus", "Ontherus sulcator",
#             "Phanaeus dejeani", "Phanaeus splendidulus", "Pseudocanthon aff. felix", "Scybalocanthon nigriceps", "Sylvicanthon aff. foveiventris"),
#           tip = c(1:35), adj = -.015, family = "serif", cex = 1, frame = "n", font = 3)
# 
# dev.off()
# 
# #########################
# # Create NFS only dendrogram
# ########################
# 
# # Prunes clust.phy of any species that aren't NFS
# 
# # subset tip labels to create a list of only non-NFS species
# tiplabs <- clust.phy$tip.label
# tips.remove <- subset(tiplabs, !(tiplabs %in% NFS.spec.n.2))
# 
# # trim our tree to get rid of non-NFS species
# NFS.dend.phy <- drop.tip(clust.phy, tips.remove, trim.internal = TRUE)
# plot.phylo(NFS.dend.phy)
# 
# #########################
# # Create FS only dendrogram
# ########################
# 
# # Prunes clust.phy of any species that aren't FS
# 
# # subset tip labels to create a list of only non-FS species
# tiplabs <- clust.phy$tip.label
# tips.remove <- subset(tiplabs, !(tiplabs %in% FS.spec.n.2))
# 
# # trim our tree to get rid of non-FS species
# FS.dend.phy <- drop.tip(clust.phy, tips.remove, trim.internal = TRUE)
# plot.phylo(FS.dend.phy)
# 
# ###########################
# # calculate FDses for NFS
# ###########################
# 
# #Using abundance
# FD.ses.NFS <- ses.pd(abun.temp.NFS, NFS.dend.phy, null.model = "frequency", runs = 999, iterations = 1000)
# # gets an error message about being unable to calculate FD if there is a single species, that is ok.
# 
# # normalize FDses output so that it is bound between 0 and 1
# x <- FD.ses.NFS$pd.obs.z
# FD.ses.NFS$FDses_norm = (x-min(x))/(max(x)-min(x))
# 
# # write a temp output, read input to save time and calculation--just read it in for now instead of re-calculating.
# write.csv(FD.ses.NFS, "FD.ses.NFS.csv")

# read in temp FD.ses file to save time
FD.ses.NFS <- read.csv("FD.ses.NFS.csv", header = TRUE)
rownames(FD.ses.NFS) <- FD.ses.NFS$X
FD.ses.NFS$X <- NULL

# ###########################
# # calculate FDses for FS
# ###########################
# 
# #Using abundance
# FD.ses.FS <- ses.pd(abun.temp.FS, FS.dend.phy, null.model = "frequency", runs = 999, iterations = 1000)
# # gets an error message about being unable to calculate PD if there is a single species, that is ok.
# 
# # normalize FDses output so that it is bound between 0 and 1
# x <- FD.ses.FS$pd.obs.z
# FD.ses.FS$FDses_norm = (x-min(x))/(max(x)-min(x))
# 
# # write a temp output, read input to save time and calculation--just read it in for now instead of re-calculating.
# write.csv(FD.ses.FS, "FD.ses.FS.csv")

# read in temp FD.ses file to save time
FD.ses.FS <- read.csv("FD.ses.FS.csv", header = TRUE)
rownames(FD.ses.FS) <- FD.ses.FS$X
FD.ses.FS$X <- NULL

# ###################################### Manip Phylogenetic ####################################################
# 
# # Step 1: Create phylogenetic tree
# # Step 2: Calculate PD
# 
# ################## Step 1: Create phylogenetic tree
# 
# ########################
# # Inspect Monaghan tree
# ########################
# 
# # plot the tree
# plot.phylo(base.phy, cex = .3)
# 
# # Does it have branch lengths?
# base.phy$edge.length
# # yes
# 
# #####################################
# # Replace code names with genus names
# #####################################
# 
# # Print list of all tip labels
# # Species are represented by a code made of abbreviations and BMNH
# # We will match to genus using BMNH
# base.phy$tip.labels
# 
# # Clean up Monaghan input file
# # get rid of all columns but genus and BMNH
# code.df$Sequence <- NULL
# code.df$Matrix.and.Tree.label <- NULL
# code.df$species <- NULL
# code.df$Genus.and.species <-NULL
# 
# # Create a vector from the tip labels with the Monaghan code truncated so that there is only BMNH
# splitlabel <- strsplit(base.phy$tip.label, "_")
# new.code <- sapply(splitlabel, '[', 1)
# # this is still in same order as the tip labels
# 
# # use this BMNH only vector to replace tip labels
# base.phy$tip.label <- new.code
# 
# #replace BMNH with genus by matching with Monaghan dataframe
# base.phy$tip.label <- code.df[[2]][match(base.phy$tip.label, code.df[[1]])]
# 
# #this leaves us with one NA at position 157--examining the original data file to look for the BMNH value 676985 reveals that it was probably a typo, because the code column inclues 676985 in the code, while the BMNH column has 673985. Regardless, this is for genus Helictopleurus, not of concern to us.
# 
# ################################
# # Prune the tree to only include genera of interest
# ###############################
# 
# # extract a list of our study genera from the traits database
# traits.df1 <- subset(traits.df, traits.df$Reference=="Nichols_unpub_2016")
# genera <- unique(traits.df1$Genus)
# 
# # subset tip labels to create a list of only genera not in our study
# tiplabs <- base.phy$tip.label
# tips.remove <- subset(tiplabs, !(tiplabs %in% genera))
# 
# # trim our tree to get rid of genera not in our study
# pruned.phy <- drop.tip(base.phy, tips.remove, trim.internal = TRUE)
# plot.phylo(pruned.phy)
# 
# # check only genera in our study included
# pruned.phy$tip.label
# 
# ##############################
# # Prune the tree so that there is only one tip for each genus
# #############################
# 
# # From the step before we have a tree with multiple repetitions of each genus
# # Here we get rid of specified tips so there is only one tip for each genus name
# # Some were easy to collapse down to a node, for others informed decisions were made about which group to keep based on which position they should be in
# # Initially I keep one extra Dichotomius to use for node
# 
# pruned.phy2 <- drop.tip(pruned.phy, c(1:3, 10:12, 15:19, 22:28, 30:32, 35, 37:64) , trim.internal = TRUE)
# plot.phylo(pruned.phy2)
# 
# # Make sure branch lengths changed
# edgelabels(pruned.phy2$edge.length)
# plot.phylo(pruned.phy)
# edgelabels(pruned.phy$edge.length)
# #yes, a comparison of these reveals that branch lengths changed between figures to reflect a new hierarchy.
# 
# ##############################
# # Add in genera missing from Monaghan
# #############################
# 
# # Here we add in genera from our study missing from Monaghan
# # Informed decisions were made about where to place genera
# # Branch lengths were made to match sister genera of those inserted
# 
# plot.phylo(pruned.phy2)
# 
# # Visualize nodes so know where to insert
# # Note that node #s change each time tree does, so this is necessary at every step
# nodelabels()
# 
# #Visualize edge lengths to match with sister genera
# edgelabels(pruned.phy2$edge.length)
# 
# 
# btree.1 <- bind.tip(pruned.phy2, "Paracanthon", edge.length = 23.355309, where = 26)
# plot.phylo(btree.1)
# nodelabels()
# 
# btree.2 <- bind.tip(btree.1, "Psuedocanthon", edge.length = 23.355309, where = 27)
# plot.phylo(btree.2)
# nodelabels()
# 
# btree.3 <- bind.tip(btree.2, "Sylvicanthon", edge.length = 23.355309, where = 28)
# plot.phylo(btree.3)
# nodelabels()
# edgelabels(btree.3$edge.length)
# 
# btree.4 <- bind.tip(btree.3, "Sulcophanaeus", edge.length = 35.966992, where = 27)
# plot.phylo(btree.4)
# nodelabels()
# edgelabels(btree.4$edge.length)
# 
# btree.5 <- bind.tip(btree.4, "Chalcocopris", edge.length = 17.027549, where = 25)
# plot.phylo(btree.5)
# nodelabels()
# 
# #Drop extra Dichotomius
# btree.5$tip.label
# graft.phy <- drop.tip(btree.5, c(2) , trim.internal = TRUE)
# plot.phylo(graft.phy)
# 
# #IMPORTANT NOTE: the tip.label vector no longer represents tips in the order they appeared on the phylogeny, but rather in the order in which they were added.
# 
# #####################
# # Attaching species names
# #####################
# 
# # Now we want to attach species names onto each genus.
# # Species branches will have a unit length of 1, we can make no assumptions.
# # Basically, we want to create a polytomy for species with branch length 1 at each genus tip.
# 
# #Creating a new tree with polytomies at genus tips
# 
# # list of genera for all of our species, in order with repeating genus names
# genera.re <- traits.df1$Genus
# # list of the names that we want on the tips of our new tree, matching with genera
# species.re <- traits.df1$Code
# 
# # create a tree with species polytomies attached to genus tips, species branches have length one
# species.phy <- make.composite.with.polytomies.2(graft.phy, genera.re, species.re)
# 
# # plot it
# plot.phylo(species.phy)
# # visualize branch lengths
# edgelabels(species.phy$edge.length)
# 
# # compare with genus phylogeny
# plot.phylo(graft.phy)
# edgelabels(graft.phy$edge.length)
# 
# # problem: while polytomies have been added with branch length one, this code simply renames genus branches for singleton species
# # We can manually add one onto the branch length of all singleton species to be consistent
# 
# # visualize the plot, first with the number of each branch and then the length of each branch.
# # this way, for each singleton species, we know which element to replace in the edge.length vector, and we can replace it with the corresponding branch length + 1.
# plot.phylo(species.phy)
# edgelabels()
# edgelabels(species.phy$edge.length)
# 
# #ontherus sulcator
# species.phy$edge.length[80] = 66.842983
# #trichillium
# species.phy$edge.length[79] = 64.218118
# #scybalocanthon
# species.phy$edge.length[59] = 24.355309
# #pseudocanthon
# species.phy$edge.length[46] = 24.355309
# #sylvicanthon
# species.phy$edge.length[45] = 24.355309
# #sulcophanaeus
# species.phy$edge.length[35] = 36.966992
# #chalcocopris
# species.phy$edge.length[7] = 18.027549
# 
# # check that branch lengths have been changed
# plot.phylo(species.phy)
# 
# # pseudocanthon still only has its genus name and not species name. This was true before I manipulated the polytomy code as well. However, the code name is in the "species" input vecotr.
# # replace manually for now
# species.phy$tip.label
# species.phy$tip.label[4] = "Pseudocanthon.aff.felix"
# plot.phylo(species.phy)
# 
# # only keep species that are NFS and FS captured in forest.
# tiplabs <- species.phy$tip.label
# tip.remove.1 <- subset(tiplabs, !(tiplabs %in% NFS.spec.n.2))
# tip.remove.2 <- subset(tip.remove.1, !(tip.remove.1 %in% FS.spec.n.2))
# species.phy <- drop.tip(species.phy, tip.remove.2, trim.internal = TRUE)
# 
# # save a PDF
# filename <- paste(figures, "Nichols_2016_tree.pdf", sep = "/")
# pdf(file = filename, width = 8, height = 15)
# 
# par(2,2,2,5)
# 
# plot.phylo(species.phy, show.tip.label = F, x.lim = c(0,204.4628))
# tiplabels(text = c("Scybalocanthon nigriceps", "Ontherus sulcator", "Pseudocanthon aff. felix", "Sylvicanthon aff. foveiventris",
#                    "Chalcocopris hesperus", "Canthon amabilis", "Canthon aff. angularis", "Canthon ibarragrassoi",
#                    "Canthon aff. luctuosus", "Canthon podagricus", "Canthon septemmaculatus", "Canthon virens chalybaeus",
#              "Coprophaneus bellicosus", "Coprophanaeus cerberus", "Coprophanaeus saphirinus",
#             "Deltochilum brasilense", "Deltochilum dentipes", "Deltochilum furcatum", "Deltochilum morbillosum",
#             "Deltochilum rubripenne", "Dichotomius assifer", "Dichotomius bos", "Dichotomius aff. carbonarius",
#             "Dichotomius depressicollis", "Dichotomius aff. fissus", "Dichotomius mormon", "Dichotomius aff. carbonarius sp.2",
#             "Dichotomius quadrinodosus", "Eurysternus cyanescens", "Eurysternus francinae",
#             "Eurysternus hirtellus", "Eurysternus inflexus", "Eurysternus parallelus",
#             "Phanaeus dejeani", "Phanaeus splendidulus"),
#           tip = c(1:35), adj = -.015, family = "serif", cex = 1.5, frame = "n", font = 3)
# 
# dev.off()
# 
# # Print a tiff file
# filename <- paste(figures, "Nichols_2016_tree.tiff", sep = "/")
# tiff(filename = filename, width = 2700, height = 3750, compression = c("none"), res = 400)
# 
# plot.phylo(species.phy, show.tip.label = F, x.lim = c(0,204.4628))
# tiplabels(text = c("Scybalocanthon nigriceps", "Ontherus sulcator", "Pseudocanthon aff. felix", "Sylvicanthon aff. foveiventris",
#                    "Chalcocopris hesperus", "Canthon amabilis", "Canthon aff. angularis", "Canthon ibarragrassoi",
#                    "Canthon aff. luctuosus", "Canthon podagricus", "Canthon septemmaculatus", "Canthon virens chalybaeus",
#              "Coprophaneus bellicosus", "Coprophanaeus cerberus", "Coprophanaeus saphirinus",
#             "Deltochilum brasilense", "Deltochilum dentipes", "Deltochilum furcatum", "Deltochilum morbillosum",
#             "Deltochilum rubripenne", "Dichotomius assifer", "Dichotomius bos", "Dichotomius aff. carbonarius",
#             "Dichotomius depressicollis", "Dichotomius aff. fissus", "Dichotomius mormon", "Dichotomius aff. carbonarius sp.2",
#             "Dichotomius quadrinodosus", "Eurysternus cyanescens", "Eurysternus francinae",
#             "Eurysternus hirtellus", "Eurysternus inflexus", "Eurysternus parallelus",
#             "Phanaeus dejeani", "Phanaeus splendidulus"),
#           tip = c(1:35), adj = -.015, family = "serif", cex = 1.5, frame = "n", font = 3)
# 
# dev.off()
# 
# #write it out
# filename <- paste(results, "Nichols_2016_tree.tre", sep = "/")
# write.tree(species.phy, file = filename, digits = 8)
# 
# #########################
# # Create NFS only phlogeny
# ########################
# 
# # Prunes species.phy of any species that aren't NFS
# # Use NFS.spec.n.2, a list of all NFS species captured in forest calculated in FD manip
# 
# # subset tip labels to create a list of only non-NFS species
# tiplabs <- species.phy$tip.label
# tips.remove <- subset(tiplabs, !(tiplabs %in% NFS.spec.n.2))
# 
# # trim our tree to get rid of non-NFS species
# NFS.phylo.phy <- drop.tip(species.phy, tips.remove, trim.internal = TRUE)
# plot.phylo(NFS.phylo.phy)
# 
# NFS.phylo.phy$tip.label
# 
# #########################
# # Create FS only phylogeny
# ########################
# 
# # Prunes species.phy of any species that aren't FS
# # Use FS.spec.n.2, a list of all FS species captured in forest calculated in FD manip
# 
# # subset tip labels to create a list of only non-FS species
# tiplabs <- species.phy$tip.label
# tips.remove <- subset(tiplabs, !(tiplabs %in% FS.spec.n.2))
# 
# # trim our tree to get rid of non-FS species
# FS.phylo.phy <- drop.tip(species.phy, tips.remove, trim.internal = TRUE)
# plot.phylo(FS.phylo.phy)
# 
# FS.phylo.phy$tip.label
# 
# 
# ################################ Step 2: Calculate PD
# 
# ####################################
# # calculation of PDses for NFS only
# ####################################
# 
# #Using abundance
# PD.ses.NFS <- ses.pd(abun.temp.NFS, NFS.phylo.phy, null.model = "frequency", runs = 999, iterations = 1000)
# 
# # normalize FDses output so that it is bound between 0 and 1
# x <- PD.ses.NFS$pd.obs.z
# PD.ses.NFS$PDses_norm = (x-min(x))/(max(x)-min(x))
# 
# # write a temp output, read input to save time and calculation--just read it in for now instead of re-calculating.
# write.csv(PD.ses.NFS, "PD.ses.NFS.csv")

# read in temp FD.ses file to save time
PD.ses.NFS <- read.csv("PD.ses.NFS.csv", header = TRUE)
rownames(PD.ses.NFS) <- PD.ses.NFS$X
PD.ses.NFS$X <- NULL 

# ###########################
# # calculate PDses for FS
# ###########################
# 
# #Using abundance
# PD.ses.FS <- ses.pd(abun.temp.FS, FS.phylo.phy, null.model = "frequency", runs = 999, iterations = 1000)
# 
# # normalize FDses output so that it is bound between 0 and 1
# x <- PD.ses.FS$pd.obs.z
# PD.ses.FS$PDses_norm = (x-min(x))/(max(x)-min(x))
# 
# # write a temp output, read input to save time and calculation--just read it in for now instead of re-calculating.
# write.csv(PD.ses.FS, "PD.ses.FS.csv")

# read in temp FD.ses file to save time
PD.ses.FS <- read.csv("PD.ses.FS.csv", header = TRUE)
rownames(PD.ses.FS) <- PD.ses.FS$X
PD.ses.FS$X <- NULL 

############################################################### Analysis #####################################################################

#################################### Step 1. AIC landscape scale selection ###############

##############
# Set up data
#############

# create FC input data for the models, summarized on the point scale
#subset buffer.df to forest
buffer.forest <- subset(buffer.df, buffer.df$Class=="forest")
#ddply to get rid of trap
buffer.forest.point <- ddply(buffer.forest, .(Landscape, Point, Buffer_m), summarize, FC=mean(Per_area))

# merge TD with forest cover data
TD.FS.forest <- merge(TD.FS, buffer.forest.point, by = c("Landscape", "Point"), all = FALSE)
TD.NFS.forest <- merge(TD.NFS, buffer.forest.point, by = c("Landscape", "Point"), all = FALSE)

####################
# Linear FS scale models
###################

# poisson regression
# modeling N vs. FC at each possible buffer scale, with LS as a random effect

# convert LS to a factor for models
TD.FS.forest$Landscape.f <- factor(TD.FS.forest$Landscape)
is.factor(TD.FS.forest$Landscape.f)

m200   = glmmadmb(formula = N ~ FC + (1|Landscape.f), data = subset(TD.FS.forest, TD.FS.forest$Buffer_m=="200"), family = "poisson")
m400   = glmmadmb(formula = N ~ FC + (1|Landscape.f), data = subset(TD.FS.forest, TD.FS.forest$Buffer_m=="400"), family = "poisson")
m800   = glmmadmb(formula = N ~ FC + (1|Landscape.f), data = subset(TD.FS.forest, TD.FS.forest$Buffer_m=="800"), family = "poisson")
m2000  = glmmadmb(formula = N ~ FC + (1|Landscape.f), data = subset(TD.FS.forest, TD.FS.forest$Buffer_m=="2000"), family = "poisson")
m3000  = glmmadmb(formula = N ~ FC + (1|Landscape.f), data = subset(TD.FS.forest, TD.FS.forest$Buffer_m=="3000"), family = "poisson")

# all converge

# Arrange AIC values in table
ms <- model.sel(m200, m400, m800, m2000, m3000)

# clean up, add, and order columns
ms$Buffer = rownames(ms)
rownames(ms) <- NULL
ms$`(Intercept)` <- NULL
ms$family <- NULL
# extracting p values for FC
ms$pvals <- c(data.frame(coef(summary(m200)))[2,4],data.frame(coef(summary(m400)))[2,4],data.frame(coef(summary(m800)))[2,4],data.frame(coef(summary(m2000)))[2,4],data.frame(coef(summary(m3000)))[2,4])
ms$Hab_use <- c("FS", " ", " "," ", " ")

table <- ms[moveme(names(ms), "Hab_use first; Buffer after Hab_use; pvals last; FC before pvals")]

aic.names = c("Class", "Buffer","k", "LL","AICc","delta","Weight", "FC", "p-values")
colnames(table) <- aic.names
table[,4:9] <-round(table[,4:9], 3)
table.TD.FS = table[with(table, order(delta)),]

#####################
# Linear NFS scale models
####################

# poisson regressions
# modeling N vs. FC at each possible buffer scale, with LS as a random effect

# convert LS to a factor for models
TD.NFS.forest$Landscape.f <- factor(TD.NFS.forest$Landscape)
is.factor(TD.NFS.forest$Landscape.f)

m200   = glmmadmb(formula = N ~ FC + (1|Landscape.f), data = subset(TD.NFS.forest, TD.NFS.forest$Buffer_m=="200"), family = "poisson")
m400   = glmmadmb(formula = N ~ FC + (1|Landscape.f), data = subset(TD.NFS.forest, TD.NFS.forest$Buffer_m=="400"), family = "poisson")
m800   = glmmadmb(formula = N ~ FC + (1|Landscape.f), data = subset(TD.NFS.forest, TD.NFS.forest$Buffer_m=="800"), family = "poisson")
m2000  = glmmadmb(formula = N ~ FC + (1|Landscape.f), data = subset(TD.NFS.forest, TD.NFS.forest$Buffer_m=="2000"), family = "poisson")
m3000  = glmmadmb(formula = N ~ FC + (1|Landscape.f), data = subset(TD.NFS.forest, TD.NFS.forest$Buffer_m=="3000"), family = "poisson")


# all converge

# Arrange AIC values in table
ms <- model.sel(m200, m400, m800, m2000, m3000)

# clean up, add, and order columns
ms$Buffer = rownames(ms)
rownames(ms) <- NULL
ms$`(Intercept)` <- NULL
ms$family <- NULL
# extracting p values for FC
ms$pvals <- c(data.frame(coef(summary(m200)))[2,4],data.frame(coef(summary(m400)))[2,4],data.frame(coef(summary(m800)))[2,4],data.frame(coef(summary(m2000)))[2,4],data.frame(coef(summary(m3000)))[2,4])
ms$Hab_use <- c("NFS", " ", " "," ", " ")

table <- ms[moveme(names(ms), "Hab_use first; Buffer after Hab_use; pvals last; FC before pvals")]

aic.names = c("Class", "Buffer","k", "LL","AICc","delta","Weight", "FC", "p-values")
colnames(table) <- aic.names
table[,4:9] <-round(table[,4:9], 3)
table.TD.NFS = table[with(table, order(delta)),]

######################################
# Direct export AIC table for FS and NFS
######################################

table.1 <- rbind(table.TD.FS, table.TD.NFS)

filename <- paste(results, "AIC_table_abundance.doc", sep = "/")
AIC_table <- RTF(file = filename)
setFontSize(AIC_table, 12)
addParagraph(AIC_table, "Table 1. AIC selection table for FS and NFS models \n")
addTable(AIC_table, as.data.frame(table.1))
done(AIC_table)

######################################### Step 2. Diversity Models #####################################

#############################
# Clean up forest buffer data to model with
##############################

# use FC input data for the models, summarized on the point scale

# extract 200m (landscape-scale) FC
buffer.forest.200 <- subset(buffer.forest.point, buffer.forest.point$Buffer_m==200)
buffer.forest.200$FC_LS <- buffer.forest.200$FC
buffer.forest.200$FC <- NULL
buffer.forest.200$Buffer_m <- NULL

##########################
# Get lat long data for each point
##########################

df_melt2.M <- subset(df_melt2, df_melt2$Habitat=="M")

# make lat/long columns numeric
df_melt2.M$Lat.n <- as.numeric(df_melt2.M$Lat)
df_melt2.M$Long.n <- as.numeric(df_melt2.M$Long)

# there are some NAs in lat/long which will mess up calculation, remove them
df_melt2.M.2 <- na.omit(df_melt2.M)

# for each point, lat and long are the average lat and long of all traps
point.coords <- ddply(df_melt2.M.2, .(Landscape, Point), summarize, 
                      lat = mean(Lat.n), 
                      long = mean(Long.n))

#############################
# clean up fragstats configuration data
#############################

# ddply per point
land.data.pt <- ddply(land.config, .(Landscape, Point), summarize,
                      TE = mean(TE))

#################################
# merge them all together into one geospatial dataframe
#############################

geos.1 <- merge(land.data.pt, point.coords, by =c("Landscape", "Point"))
geos.final <- merge(geos.1, buffer.forest.200, by =c("Landscape", "Point"))

##########################
# Merge with all biodiversity metrics
###############################

# TD--raw richness, raw abundance
TD.FS.test <- NULL
TD.FS.test$Landscape <- TD.FS$Landscape
TD.FS.test$Point <- TD.FS$Point
TD.FS.test$TD_FS_raw <- TD.FS$TD
TD.FS.test$TD_FS_abund <- TD.FS$N
TD.FS.test.df <- as.data.frame(TD.FS.test)

TD.NFS.test <- NULL
TD.NFS.test$Landscape <- TD.NFS$Landscape
TD.NFS.test$Point <- TD.NFS$Point
TD.NFS.test$TD_NFS_raw <- TD.NFS$TD
TD.NFS.test$TD_NFS_abund <- TD.NFS$N
TD.NFS.test.df <- as.data.frame(TD.NFS.test)

# FD--use normalized ses that is standardized
FD.ses.NFS$ID <- rownames(FD.ses.NFS)
splitID <- strsplit(FD.ses.NFS$ID, "-")
landscape <- sapply(splitID, '[', 1)
point <- sapply(splitID, '[', 2)
habitat <- sapply(splitID, '[', 3)
FD.ses.NFS$Landscape <- landscape
FD.ses.NFS$Point <- point
FD.ses.NFS$Habitat <- habitat
FD.ses.NFS$FDses_stand = (FD.ses.NFS$FDses_norm * (length(FD.ses.NFS$FDses_norm) - 1) + 0.5) / length(FD.ses.NFS$FDses_norm)
FD.test.NFS <- NULL
FD.test.NFS$Landscape <- FD.ses.NFS$Landscape
FD.test.NFS$Point <- FD.ses.NFS$Point
FD.test.NFS$FD_NFS <- FD.ses.NFS$FDses_stand
FD.test.NFS.df <- as.data.frame(FD.test.NFS)

FD.ses.FS$ID <- rownames(FD.ses.FS)
splitID <- strsplit(FD.ses.FS$ID, "-")
landscape <- sapply(splitID, '[', 1)
point <- sapply(splitID, '[', 2)
habitat <- sapply(splitID, '[', 3)
FD.ses.FS$Landscape <- landscape
FD.ses.FS$Point <- point
FD.ses.FS$Habitat <- habitat
FD.ses.FS$FDses_stand = (FD.ses.FS$FDses_norm * (length(FD.ses.FS$FDses_norm) - 1) + 0.5) / length(FD.ses.FS$FDses_norm)
FD.test.FS <- NULL
FD.test.FS$Landscape <- FD.ses.FS$Landscape
FD.test.FS$Point <- FD.ses.FS$Point
FD.test.FS$FD_FS <- FD.ses.FS$FDses_stand
FD.test.FS.df <- as.data.frame(FD.test.FS)

# PD--use normalized ses that is standardized
PD.ses.NFS$ID <- rownames(PD.ses.NFS)
splitID <- strsplit(PD.ses.NFS$ID, "-")
landscape <- sapply(splitID, '[', 1)
point <- sapply(splitID, '[', 2)
habitat <- sapply(splitID, '[', 3)
PD.ses.NFS$Landscape <- landscape
PD.ses.NFS$Point <- point
PD.ses.NFS$Habitat <- habitat
PD.ses.NFS$PDses_stand = (PD.ses.NFS$PDses_norm * (length(PD.ses.NFS$PDses_norm) - 1) + 0.5) / length(PD.ses.NFS$PDses_norm)
PD.test.NFS <- NULL
PD.test.NFS$Landscape <- PD.ses.NFS$Landscape
PD.test.NFS$Point <- PD.ses.NFS$Point
PD.test.NFS$PD_NFS <- PD.ses.NFS$PDses_stand
PD.test.NFS.df <- as.data.frame(PD.test.NFS)

PD.ses.FS$ID <- rownames(PD.ses.FS)
splitID <- strsplit(PD.ses.FS$ID, "-")
landscape <- sapply(splitID, '[', 1)
point <- sapply(splitID, '[', 2)
habitat <- sapply(splitID, '[', 3)
PD.ses.FS$Landscape <- landscape
PD.ses.FS$Point <- point
PD.ses.FS$Habitat <- habitat
PD.ses.FS$PDses_stand = (PD.ses.FS$PDses_norm * (length(PD.ses.FS$PDses_norm) - 1) + 0.5) / length(PD.ses.FS$PDses_norm)
PD.test.FS <- NULL
PD.test.FS$Landscape <- PD.ses.FS$Landscape
PD.test.FS$Point <- PD.ses.FS$Point
PD.test.FS$PD_FS <- PD.ses.FS$PDses_stand
PD.test.FS.df <- as.data.frame(PD.test.FS)

# merge
# all.x here is important, it makes sure that all sites are preserved even when certain metrics lack values at them. For example, at some sites FS species were captured by no NFS species, so all biodiversity metrics for NFS are NA.
model.df.1 <- merge(geos.final, TD.NFS.test.df, by=c("Landscape", "Point"), all.x = TRUE) 
model.df.1.5 <- merge(model.df.1, TD.FS.test.df, by=c("Landscape", "Point"), all.x = TRUE )
model.df.2 <- merge(model.df.1.5, FD.test.NFS.df, by=c("Landscape", "Point"), all.x = TRUE)
model.df.3 <- merge(model.df.2, FD.test.FS.df, by=c("Landscape", "Point"), all.x = TRUE)
model.df.4 <- merge(model.df.3, PD.test.NFS.df, by=c("Landscape", "Point"), all.x = TRUE)
model.df.5 <- merge(model.df.4, PD.test.FS.df, by=c("Landscape", "Point"), all.x = TRUE) # 335-2 is included here even though no species were captured.

# remove 335-2
model.df.final <- subset(model.df.5, !(model.df.5$Landscape=="335" & model.df.5$Point=="2"))

# create a column in model.df.final with ID
model.df.final$ID <- paste(model.df.final$Landscape, "-", model.df.final$Point, "-M", sep = "")

#extract list of points in final model
points.final <- model.df.final$ID

# subset for each non-TD metric, to deal with NAs bc they can't be fed into models
model.df.final.NFS <- subset(model.df.final, !(model.df.final$FD_NFS=="NA"))
model.df.final.FS <- subset(model.df.final, !(model.df.final$FD_FS=="NA"))


#######################################
# Models
################################

# modeling biodiversity metrics vs. FC, long, TE with LS and a random effect
# beta mixed models for FD and PD
# poisson mixed models for N and TD

# convert LS to a factor
model.df.final$Landscape.f <- as.factor(model.df.final$Landscape)
model.df.final.FS$Landscape.f <- as.factor(model.df.final.FS$Landscape)
model.df.final.NFS$Landscape.f <- as.factor(model.df.final.NFS$Landscape)

# check for correlation between indp. vars using a scatterplot matrix
pairs(model.df.final[,3:9])

# load a later version of glmm ADMB to solve convergence problems
install.packages("glmmADMB_0.7.tar.gz",repos=NULL,type="source")
library(glmmADMB) 

# everything is beta distributed, so we will use glmmADMB
TD.raw.FS.mod    = glmmadmb(formula = TD_FS_raw     ~ FC_LS + long +  TE + (1|Landscape.f), data = model.df.final.FS,  family = "poisson") # all get this error message: Estimated covariance matrix may not be positive definite
TD.abund.FS.mod  = glmmadmb(formula = TD_FS_abund   ~ FC_LS + long +  TE + (1|Landscape.f), data = model.df.final.FS,  family = "poisson")
TD.raw.NFS.mod   = glmmadmb(formula = TD_NFS_raw    ~ FC_LS + long +  TE + (1|Landscape.f), data = model.df.final.NFS, family = "poisson")
TD.abund.NFS.mod = glmmadmb(formula = TD_NFS_abund  ~ FC_LS + long +  TE + (1|Landscape.f), data = model.df.final.NFS, family = "poisson")
FD.FS.mod        = glmmadmb(formula = FD_FS         ~ FC_LS + long +  TE + (1|Landscape.f), data = model.df.final.FS,  family = "beta", link = "logit") 
FD.NFS.mod       = glmmadmb(formula = FD_NFS        ~ FC_LS + long +  TE + (1|Landscape.f), data = model.df.final.NFS, family = "beta", link = "logit")
PD.FS.mod        = glmmadmb(formula = PD_FS         ~ FC_LS + long +  TE + (1|Landscape.f), data = model.df.final.FS,  family = "beta", link = "logit") 
PD.NFS.mod       = glmmadmb(formula = PD_NFS        ~ FC_LS + long +  TE + (1|Landscape.f), data = model.df.final.NFS, family = "beta", link = "logit")

# summarize model results
summary(TD.raw.FS.mod)
summary(TD.abund.FS.mod)
summary(TD.raw.NFS.mod)
summary(TD.abund.NFS.mod)
summary(FD.FS.mod)
summary(FD.NFS.mod)
summary(PD.FS.mod)
summary(PD.NFS.mod)

############################################
# write out a table of model coefficients
#########################################

# extract model coefficients
TD.raw.FS.sum <- data.frame(coef(summary(TD.raw.FS.mod)))
TD.raw.FS.sum$model <- "TD.raw_FS"
TD.raw.FS.sum$var <- rownames(TD.raw.FS.sum)

TD.raw.NFS.sum <- data.frame(coef(summary(TD.raw.NFS.mod)))
TD.raw.NFS.sum$model <- "TD.raw_NFS"
TD.raw.NFS.sum$var <- rownames(TD.raw.NFS.sum)

TD.abund.FS.sum <- data.frame(coef(summary(TD.abund.FS.mod)))
TD.abund.FS.sum$model <- "TD.abund_FS"
TD.abund.FS.sum$var <- rownames(TD.abund.FS.sum)

TD.abund.NFS.sum <- data.frame(coef(summary(TD.abund.NFS.mod)))
TD.abund.NFS.sum$model <- "TD.abund_NFS"
TD.abund.NFS.sum$var <- rownames(TD.abund.NFS.sum)

FD.FS.sum <- data.frame(coef(summary(FD.FS.mod)))
FD.FS.sum$model <- "FD_FS"
FD.FS.sum$var <- rownames(FD.FS.sum)

FD.NFS.sum <- data.frame(coef(summary(FD.NFS.mod)))
FD.NFS.sum$model <- "FD_NFS"
FD.NFS.sum$var <- rownames(FD.NFS.sum)

PD.FS.sum <- data.frame(coef(summary(PD.FS.mod)))
PD.FS.sum$model <- "PD_FS"
PD.FS.sum$var <- rownames(PD.FS.sum)

PD.NFS.sum <- data.frame(coef(summary(PD.NFS.mod)))
PD.NFS.sum$model <- "PD_NFS"
PD.NFS.sum$var <- rownames(PD.NFS.sum)

# rbind all models together
mod.sum.0 <- rbind(TD.abund.FS.sum, TD.abund.NFS.sum)
mod.sum.0.1 <- rbind(mod.sum.0, TD.raw.FS.sum)
mod.sum.0.2 <- rbind(mod.sum.0.1, TD.raw.NFS.sum)
mod.sum.2 <- rbind(mod.sum.0.2, FD.FS.sum)
mod.sum.3 <- rbind(mod.sum.2, FD.NFS.sum)
mod.sum.4 <- rbind(mod.sum.3, PD.FS.sum)
mod.sum.5 <- rbind(mod.sum.4, PD.NFS.sum)

# clean up unnecessary columns and rows intercept, add label columns
rownames(mod.sum.5) <- NULL
mod.sum.5$Std..Error <- NULL
mod.sum.5$z.value <- NULL

mod.sum.6 <- subset(mod.sum.5, !(mod.sum.5$var=="(Intercept)"))

mod.sum.6$metric <- c("N","", "", "", "", "", "TD", "", "", "", "", "", "FD", "", "", "", "", "","PD", "", "", "", "", "")
mod.sum.6$hab_use <- c("FS", "", "", "NFS", "","","FS", "", "", "NFS", "","","FS","", "", "NFS", "","","FS", "", "", "NFS", "","")
mod.sum.6$model <- NULL

# replace variable abbreviations with full names
mod.sum.7 <- lapply(mod.sum.6, function(x) {
  gsub("FC_LS", "Forest Cover", x)
})
mod.sum.8 <- lapply(mod.sum.7, function(x) {
  gsub("long", "Longitude", x)
})
mod.sum.9 <- lapply(mod.sum.8, function(x) {
  gsub("TE", "Total Edge", x)
})

mod.sum.9.df <- as.data.frame(mod.sum.9)

# rearrange and rename columns
table <- mod.sum.9.df[moveme(names(mod.sum.9.df), "metric first; hab_use after metric; var after hab_use")]
aic.names = c("Metric", "Classification","Variable", "Coefficient","p-value")
colnames(table) <- aic.names

# round p-vals and coefficients
table[,5] <-formatC(round(as.numeric(table[,5]), 3), format = 'f', digits = 3)
table[,4] <-signif(as.numeric(table[,4]), digits = 3)
table$Coefficient <- format(table$Coefficient, scientific = TRUE)

# print it out
filename <- paste(results, "Model_coefficients_table.doc", sep = "/")
AIC_table <- RTF(file = filename)
setFontSize(AIC_table, 12)
addParagraph(AIC_table, "Table 1. Model coefficients for all mixed models \n")
addTable(AIC_table, as.data.frame(table))
done(AIC_table)


############################################################
# make some plots of each indv var vs. metric for each metric
########################################################

filename <- paste(figures, "TD.FS.mod.pdf", sep = "/")
pdf(file = filename, width = 12, height = 12)

par(mfrow=c(2,2))
par(margin(t = 2, r = 2, b = 4, l = 6))

plot(model.df.final$FC_LS, model.df.final$TD_FS_raw, xlab = "Landscape-scale forest cover", ylab = "TD", cex.lab = 1.5, pch = 19)
plot(model.df.final$long, model.df.final$TD_FS_raw, xlab = "Longitude", ylab = "TD", cex.lab = 1.5, pch = 19)
plot(model.df.final$TE, model.df.final$TD_FS_raw, xlab = "Total Edge", ylab = "TD", cex.lab = 1.5, pch = 19)

dev.off()


filename <- paste(figures, "TD.NFS.mod.pdf", sep = "/")
pdf(file = filename, width = 12, height = 12)

par(mfrow=c(2,2))
par(margin(t = 2, r = 2, b = 4, l = 6))

plot(model.df.final$FC_LS, model.df.final$TD_NFS_raw, xlab = "Landscape-scale forest cover", ylab = "TD", cex.lab = 1.5, pch = 19)
plot(model.df.final$long, model.df.final$TD_NFS_raw, xlab = "Longitude", ylab = "TD", cex.lab = 1.5, pch = 19)
plot(model.df.final$TE, model.df.final$TD_NFS_raw, xlab = "Total Edge", ylab = "TD", cex.lab = 1.5, pch = 19)

dev.off()


filename <- paste(figures, "N.FS.mod.pdf", sep = "/")
pdf(file = filename, width = 12, height = 12)

par(mfrow=c(2,2))
par(margin(t = 2, r = 2, b = 4, l = 6))

plot(model.df.final$FC_LS, model.df.final$TD_FS_abund, xlab = "Landscape-scale forest cover", ylab = "N", cex.lab = 1.5, pch = 19)
plot(model.df.final$long, model.df.final$TD_FS_abund, xlab = "Longitude", ylab = "N", cex.lab = 1.5, pch = 19)
plot(model.df.final$TE, model.df.final$TD_FS_abund, xlab = "Total Edge", ylab = "N", cex.lab = 1.5, pch = 19)

dev.off()


filename <- paste(figures, "N.NFS.mod.pdf", sep = "/")
pdf(file = filename, width = 12, height = 12)

par(mfrow=c(2,2))
par(margin(t = 2, r = 2, b = 4, l = 6))

plot(model.df.final$FC_LS, model.df.final$TD_NFS_abund, xlab = "Landscape-scale forest cover", ylab = "N", cex.lab = 1.5, pch = 19)
plot(model.df.final$long, model.df.final$TD_NFS_abund, xlab = "Longitude", ylab = "N", cex.lab = 1.5, pch = 19)
plot(model.df.final$TE, model.df.final$TD_NFS_abund, xlab = "Total Edge", ylab = "N", cex.lab = 1.5, pch = 19)

dev.off()


filename <- paste(figures, "FD.NFS.mod.pdf", sep = "/")
pdf(file = filename, width = 12, height = 12)

par(mfrow=c(2,2))
par(margin(t = 2, r = 2, b = 4, l = 6))

plot(model.df.final$FC_LS, model.df.final$FD_NFS, xlab = "Landscape-scale forest cover", ylab = "FD", cex.lab = 1.5, pch = 19)
plot(model.df.final$long, model.df.final$FD_NFS, xlab = "Longitude", ylab = "FD", cex.lab = 1.5, pch = 19)
plot(model.df.final$TE, model.df.final$FD_NFS, xlab = "Total Edge", ylab = "FD", cex.lab = 1.5, pch = 19)

dev.off()


filename <- paste(figures, "FD.FS.mod.pdf", sep = "/")
pdf(file = filename, width = 12, height = 12)

par(mfrow=c(2,2))
par(margin(t = 2, r = 2, b = 4, l = 6))

plot(model.df.final$FC_LS, model.df.final$FD_FS, xlab = "Landscape-scale forest cover", ylab = "FD", cex.lab = 1.5, pch = 19)
plot(model.df.final$long, model.df.final$FD_FS, xlab = "Longitude", ylab = "FD", cex.lab = 1.5, pch = 19)
plot(model.df.final$TE, model.df.final$FD_FS, xlab = "Total Edge", ylab = "FD", cex.lab = 1.5, pch = 19)

dev.off()


filename <- paste(figures, "PD.NFS.mod.pdf", sep = "/")
pdf(file = filename, width = 12, height = 12)

par(mfrow=c(2,2))
par(margin(t = 2, r = 2, b = 4, l = 6))

plot(model.df.final$FC_LS, model.df.final$PD_NFS, xlab = "Landscape-scale forest cover", ylab = "PD", cex.lab = 1.5, pch = 19)
plot(model.df.final$long, model.df.final$PD_NFS, xlab = "Longitude", ylab = "PD", cex.lab = 1.5, pch = 19)
plot(model.df.final$TE, model.df.final$PD_NFS, xlab = "Total Edge", ylab = "PD", cex.lab = 1.5, pch = 19)

dev.off()


filename <- paste(figures, "PD.FS.mod.pdf", sep = "/")
pdf(file = filename, width = 12, height = 12)

par(mfrow=c(2,2))
par(margin(t = 2, r = 2, b = 4, l = 6))

plot(model.df.final$FC_LS, model.df.final$PD_FS, xlab = "Landscape-scale forest cover", ylab = "PD", cex.lab = 1.5, pch = 19)
plot(model.df.final$long, model.df.final$PD_FS, xlab = "Longitude", ylab = "PD", cex.lab = 1.5, pch = 19)
plot(model.df.final$TE, model.df.final$PD_FS, xlab = "Total Edge", ylab = "PD", cex.lab = 1.5, pch = 19)

dev.off()

####################################
# generate summary of modeling metrics
####################################

# N FS
mean(model.df.final.FS$TD_FS_abund) #11.91667
sd(model.df.final.FS$TD_FS_abund) #11.55939
min(model.df.final.FS$TD_FS_abund) #1
max(model.df.final.FS$TD_FS_abund) #56

# N NFS
mean(model.df.final.NFS$TD_NFS_abund) #18.35632
sd(model.df.final.NFS$TD_NFS_abund) #22.86664
min(model.df.final.NFS$TD_NFS_abund) #1
max(model.df.final.NFS$TD_NFS_abund) #133

# TD FS
mean(model.df.final.FS$TD_FS_raw) #3.190476
sd(model.df.final.FS$TD_FS_raw) #1.724914
min(model.df.final.FS$TD_FS_raw) #1
max(model.df.final.FS$TD_FS_raw) #8

# TD NFS
mean(model.df.final.NFS$TD_NFS_raw) #3.517241
sd(model.df.final.NFS$TD_NFS_raw) #2.106811
min(model.df.final.NFS$TD_NFS_raw) #1
max(model.df.final.NFS$TD_NFS_raw) #8

# FD FS
mean(model.df.final.FS$FD_FS) #0.4166199
sd(model.df.final.FS$FD_FS) #0.2793247
min(model.df.final.FS$FD_FS) #0.005952381
max(model.df.final.FS$FD_FS) # 0.9940476

# FD NFS
mean(model.df.final.NFS$FD_NFS) #0.3640388
sd(model.df.final.NFS$FD_NFS) #0.2522606
min(model.df.final.NFS$FD_NFS) #0.005747126
max(model.df.final.NFS$FD_NFS) #0.9942529

# PD FS
mean(model.df.final.FS$PD_FS) #0.3863029
sd(model.df.final.FS$PD_FS) #0.2764124
min(model.df.final.FS$PD_FS) #0.005952381
max(model.df.final.FS$PD_FS) # 0.9940476

# PD NFS
mean(model.df.final.NFS$PD_NFS) #0.396529
sd(model.df.final.NFS$PD_NFS) #0.2681362
min(model.df.final.NFS$PD_NFS) #0.005747126
max(model.df.final.NFS$PD_NFS) #0.9942529

# FC
mean(model.df.final.NFS$FC) #64.32046
sd(model.df.final.NFS$FC) #17.3041
min(model.df.final.NFS$FC) #23.85039
max(model.df.final.NFS$FC) #96.40319

# TE
mean(model.df.final.NFS$TE) #1496.404
sd(model.df.final.NFS$TE) #469.946
min(model.df.final.NFS$TE) #909.0566
max(model.df.final.NFS$TE) #3833.574


################################ Step 3. Permutation tests E vs. W ###########################

##################
# subset master dataframes
###################

# east 148, 266, 291, 317
# west 215, 263, 282, 329, 333, 335, 359, 399
east.FS <- subset(model.df.final.FS, model.df.final.FS$Landscape=="148" | 
                    model.df.final.FS$Landscape== "266" | 
                    model.df.final.FS$Landscape== "291" | 
                    model.df.final.FS$Landscape== "317")
west.FS <- subset(model.df.final.FS, model.df.final.FS$Landscape=="215" | 
                    model.df.final.FS$Landscape== "263" | 
                    model.df.final.FS$Landscape== "282" | 
                    model.df.final.FS$Landscape== "329" |
                    model.df.final.FS$Landscape== "333" |
                    model.df.final.FS$Landscape== "335" |
                    model.df.final.FS$Landscape== "359" |
                    model.df.final.FS$Landscape== "399")

east.NFS <- subset(model.df.final.NFS, model.df.final.NFS$Landscape=="148" | 
                     model.df.final.NFS$Landscape== "266" | 
                     model.df.final.NFS$Landscape== "291" | 
                     model.df.final.NFS$Landscape== "317")
west.NFS <- subset(model.df.final.NFS, model.df.final.NFS$Landscape=="215" | 
                     model.df.final.NFS$Landscape== "263" | 
                     model.df.final.NFS$Landscape== "282" | 
                     model.df.final.NFS$Landscape== "329" |
                     model.df.final.NFS$Landscape== "333" |
                     model.df.final.NFS$Landscape== "335" |
                     model.df.final.NFS$Landscape== "359" |
                     model.df.final.NFS$Landscape== "399")

################
# Permutation tests
################

# TD FS
x <- east.FS$TD_FS_raw
y <- west.FS$TD_FS_raw
t1 <- permTS(x, y, alternative = c("two.sided"), method = "exact.mc", 
             control = permControl(tsmethod = "abs"))
# p-value = 0.113

# TD NFS
x <- east.NFS$TD_NFS_raw
y <- west.NFS$TD_NFS_raw
t2 <- permTS(x, y, alternative = c("two.sided"), method = "exact.mc", 
             control = permControl(tsmethod = "abs"))
# p-value = 0.001

# N FS
x <- east.FS$TD_FS_abund
y <- west.FS$TD_FS_abund
t3 <- permTS(x, y, alternative = c("two.sided"), method = "exact.mc", 
             control = permControl(tsmethod = "abs"))
# p-value = 0.047

# N NFS
x <- east.NFS$TD_NFS_abund
y <- west.NFS$TD_NFS_abund
t4 <- permTS(x, y, alternative = c("two.sided"), method = "exact.mc", 
             control = permControl(tsmethod = "abs"))
# p-value = 0.002

# FD FS
x <- east.FS$FD_FS
y <- west.FS$FD_FS
t5 <- permTS(x, y, alternative = c("two.sided"), method = "exact.mc", 
             control = permControl(tsmethod = "abs"))
# p-value = 0.095

# FD NFS
x <- east.NFS$FD_NFS
y <- west.NFS$FD_NFS
t6 <- permTS(x, y, alternative = c("two.sided"), method = "exact.mc", 
             control = permControl(tsmethod = "abs"))
# p-value = 0.001

# PD FS
x <- east.FS$PD_FS
y <- west.FS$PD_FS
t7 <- permTS(x, y, alternative = c("two.sided"), method = "exact.mc", 
             control = permControl(tsmethod = "abs"))
# p-value = 0.138

# PD NFS
x <- east.NFS$PD_NFS
y <- west.NFS$PD_NFS
t8 <- permTS(x, y, alternative = c("two.sided"), method = "exact.mc", 
             control = permControl(tsmethod = "abs"))
# p-value = 0.001

########################################################
# write out a table of permutation test results
########################################################

# extract t-test coefficients into table
T.table <- NULL
T.table$metric <- c("N", "", "TD", "", "FD", "","PD", "")
T.table$class <- c("FS", "NFS","FS", "NFS","FS", "NFS","FS", "NFS")
T.table$stat <- c(t3$estimate, t4$estimate,t1$estimate, t2$estimate, t5$estimate, t6$estimate, t7$estimate,t8$estimate)
T.table$p <- c(t3$p.value, t4$p.value,t1$p.value, t2$p.value, t5$p.value, t6$p.value, t7$p.value,t8$p.value)

# make into df
T.table.df <- as.data.frame(T.table)

# rename columns
aic.names = c("Metric", "Classification","Statistic", "p-value")
colnames(T.table.df) <- aic.names

# round statistic and p-value
T.table.df[,3] <-formatC(round(T.table.df[,3], 2), format = 'f', digits = 2)
T.table.df[,4] <-formatC(round(T.table.df[,4], 3), format = 'f', digits = 3)

# print it out
filename <- paste(results, "Permutation_test_results.doc", sep = "/")
AIC_table <- RTF(file = filename)
setFontSize(AIC_table, 12)
addParagraph(AIC_table, "Table 1. Results of permutation tests testing the difference in each biodiversity metric between Eastern and Western sites for FS and NFS separately. The statistic represents Eastern minus Western sites for each metric, and repetitions represents the number of monte carlo repetitions.\n")
addTable(AIC_table, as.data.frame(T.table.df))
done(AIC_table)

######################################
# explore higher FS N values in west than east
######################################

par(mfrow=c(1,2))
hist(east.FS$TD_FS_abund, breaks = 20)
hist(west.FS$TD_FS_abund, breaks = 20) # really are a lot more high values

# who is that outlier?
subset(west.FS, west.FS$TD_FS_abund > 50)
# 263-0, abund is 56, really low LS level FC--one of the smallest 27.47394--5 spec.
# solve with species plots
# several species high there, esp. ibarragrassoi

mean(east.FS$TD_FS_abund) #8.208333
mean(west.FS$TD_FS_abund) #13.94643

####################################################################  Figures ########################################################### 

##############################################
# Reproduce Fahrig 2003 Figure 3
##############################################

# merge configuration with composition data
buffer.forest$X <- NULL
buffer.forest.200.2 <- subset(buffer.forest, buffer.forest$Buffer_m==200)
land.data <- merge(buffer.forest.200.2, land.config, by = c("Code"))


# ddply per ls
land.data.ls <- ddply(land.data, .(Landscape.x), summarize,
                      FC = mean(Per_area),
                      TE = mean(TE),
                      MPI = mean(MPI), 
                      NP = mean(NP))

#Set up plotting space and PDF command
filename <- paste(figures, "Fahrig_reproduction_landscape.pdf", sep = "/")
pdf(file = filename, width = 15, height = 5)

par(mfrow=c(1,3))

#Number of patches
plot(land.data.ls$FC, land.data.ls$NP., xlab = "Habitat Amount (%)", ylab = "Number of Patches", pch = 19)
#Mean patch size
plot(land.data.ls$FC, land.data.ls$MPI, xlab = "Habitat Amount (%)", ylab = "Mean Patch Isolation", pch =19)
#Total edge
plot(land.data.ls$FC, land.data.ls$TE, xlab = "Habitat Amount (%)", ylab = "Total Edge", pch =19)

dev.off()


#####################################################################
# Binned Banks-Leite matrix, species by site arranged by % FC
####################################################################

# Goal--reproduce the matrix from Banks-Leite 2014, showing species presence/absence at each site, with sites arranged by %FC, and a division between FS and NFS
# Key difference--we are going to bin sites by FC

# Create an input species/site matrix

# use abundance data input (outliers etc. removed)
abund.df_melt1 <- df_melt2

# summarize number of individuals of each species trapped per point
abund.sum.df = ddply(abund.df_melt1, .(Landscape, Habitat, Point, code), summarize, N = sum(N.trap))

# subset for only Forest captures
abund.sum.df2 <- subset(abund.sum.df, abund.sum.df$Habitat=="M")

# create a single ID column
abund.sum.df2$ID <- paste(abund.sum.df2$Landscape, abund.sum.df2$Point, abund.sum.df2$Habitat, sep="-")
# This will print an ID column that lists Landscape-Point-Habitat
head(abund.sum.df2)

# get rid of old ID columns
abund.sum.df2$Landscape <- NULL
abund.sum.df2$Point <- NULL
abund.sum.df2$Habitat <- NULL

# subset only for FS and NFS with non-0 captures
abund.subs <- abund.sum.df2[abund.sum.df2$code %in% spec.keep,] ######PROBLEM IS HERE--SPECIES.KEEP IS ALL.

#get rid of points of all points not used in final model 
abund.sum.df3 <- abund.subs[abund.subs$ID %in% points.final,] # 91 points

# cast back into matrix,  with sites as columns and spec rows.
PA.df = dcast(abund.sum.df3, code ~ ID, value.var = "N")
rownames(PA.df) <- PA.df$code
PA.df$code <- NULL

#convert to presence/absence
PA.df[PA.df>0] <- 1

# order df by FC
order.model.df.final <- model.df.final[order(model.df.final$FC_LS),]
order.sites <- order.model.df.final$ID

PA.order.df <- PA.df[,order.sites]


# use Banks-Leite 2014 method of ordering species by forest pref.

# input needed---
# list of % FC in order
FC_LS <- order.model.df.final$FC_LS # 85 sites
# presence/absence matrix with sites in % FC order
PA.spec.site.df <- t(PA.order.df) # 85 sites

# calculate preference of each species for % FC
weight.avg(PA.spec.site.df, FC_LS)

# Create a df that appends forest prefence weights to species names
forest.pref <- NULL
forest.pref$code <- colnames(PA.spec.site.df)
forest.pref$weights <- weight.avg(PA.spec.site.df, FC_LS)
forest.pref.df <- as.data.frame(forest.pref)

# Pick binning scale
# FC_LS has a list of all FC values from modeling df
FC_LS # Values range from 23-96

# let's try visualizing as hist to see what lvls bin at
hist(FC_LS, breaks = c(seq(20,100,4)), right = FALSE)
# all bins still filled!

# bin at 4% increase in FC

# go back to original abundance df (or PA), bin sites by FC, convert to PA
PA.test <- t(PA.df) # make a df with spec as columns, sites as rows
PA.test.df <- as.data.frame(PA.test)

# attach FC column
#make mini df with just FC and ID
LS.FC <- NULL
LS.FC$ID <- model.df.final$ID
LS.FC$FC <- model.df.final$FC_LS
LS.FC.df <- as.data.frame(LS.FC)
PA.test.df$ID <- rownames(PA.test.df)
PA.test.df2 <- merge(PA.test.df, LS.FC.df, by=c("ID"))

# create a column that lists which bin each site should be in
# bin by 4s 
PA.test.df2$FC_bin <- cut(PA.test.df2$FC, breaks=c(seq(20,100,4)), right=FALSE)
PA.test.df2$FC <- NULL
PA.test.df2$ID <- NULL

# sum each species column by bins
PA.bin <- aggregate(.~FC_bin, PA.test.df2, sum)
# convert to P/A
PA.bin[PA.bin>0] <- 1

# arrange species by forest pref, append MA column as before--copied code
rownames(PA.bin) <- PA.bin$FC_bin
PA.bin$FC_bin<- NULL
PA.bin.2 <- t(PA.bin)
PA.bin.2.df <- as.data.frame(PA.bin.2)

# add an empty bin at the beginning for graphing ease
PA.bin.2.df <- cbind(x1 = 0, PA.bin.2.df)

# subset PA.bin.2 (sites=col, spec=rows) by FS, NFS
PA.bin.2.df$code <- rownames(PA.bin.2.df)
PA.bin.order.df <- merge(PA.bin.2.df, hab.assoc, by=c("code"))
PA.FS <- subset(PA.bin.order.df, PA.bin.order.df$Classification=="FS")
PA.NFS <- subset(PA.bin.order.df, PA.bin.order.df$Classification=="NFS")

# merge weights with each subset
PA.FS.weight <- merge(PA.FS, forest.pref.df, by=c("code") )
PA.NFS.weight <- merge(PA.NFS, forest.pref.df, by=c("code") )

# order each subset by weight
PA.FS.order <- PA.FS.weight[order(-PA.FS.weight$weights),]
PA.NFS.order <- PA.NFS.weight[order(-PA.NFS.weight$weights),]

# while subsetted apart, change all 1s to 2s in PA.NFS so plot in diff colors (this may mess up ordering posthoc)
PA.NFS.order[PA.NFS.order==1] <- 2

# rbind together, clean up column names
PA.oo.bin.df <- rbind(PA.FS.order, PA.NFS.order)
rownames(PA.oo.bin.df) <- PA.oo.bin.df$code

# attach column with endemic to MA
MA.end <- hab.assoc.2[,-(1:5),drop=FALSE]
MA.end <- MA.end[,-(2:3), drop=FALSE]
MA.end <- MA.end[,-(3:24), drop=FALSE]
MA.end$end <- ifelse(MA.end$Biogeo_assoc=="MA", 3, 4)
# order in same order as species in dataframe
MA.end.subs <- subset(MA.end, MA.end$Code %in% PA.oo.bin.df$code)
MA.end.subs.order <- MA.end.subs[match(PA.oo.bin.df$code, MA.end.subs$Code),]

PA.oo.bin.df$end <- MA.end.subs.order$end

# clean up columns
PA.oo.bin.df$code <- NULL
PA.oo.bin.df$Classification <- NULL
PA.oo.bin.df$weights <- NULL

# flip columns and rows so plots correctly
PA.oo.bin.flip <- t(PA.oo.bin.df)
PA.oo.bin.m <- as.matrix(PA.oo.bin.flip)

# reverse order of columns
PA.oo.bin.m <- PA.oo.bin.m[,ncol(PA.oo.bin.m):1]

# save different name so can plot later
PA.oo.bin.m.FC <- PA.oo.bin.m

# write a PDF
filename <- paste(figures, "Matrix_FC.pdf", sep = "/")
pdf(file = filename, width = 8, height = 6)

par(mar = c(4,5,5,3))
image(PA.oo.bin.m.FC, col= c("white","dodgerblue", "brown1", "black", "gray"), axes=FALSE, main="Forest cover", ylab="Species")
axis(1, at=c(0.024, .262, .5, .738, .976), labels = c(20,40,60,80,100)) 

dev.off()

##################################
# Binned Banks-Leite matrix for longitude
#################################

# Find range of longitude values
min(model.df.final$long) #-46.51152
max(model.df.final$long) #-46.01835

# Visualize to pick binning scale

# bins at .01
hist(model.df.final$long, breaks = c(seq(-46.52,-46.01,.01)), right = FALSE)
# 6 breaks in the data

# go back to PA.df, bin sites by long, convert to PA
PA.test <- t(PA.df) # make a df with spec as columns, sites as rows
PA.test.df <- as.data.frame(PA.test)

# attach long column
#make mini df with just long and ID
LS.long <- NULL
LS.long$ID <- model.df.final$ID
LS.long$long <- model.df.final$long
LS.long.df <- as.data.frame(LS.long)
PA.test.df$ID <- rownames(PA.test.df)
PA.test.df2 <- merge(PA.test.df, LS.long.df, by=c("ID"))

# create a column that lists which bin each site should be in
PA.test.df2$long_bin <- cut(PA.test.df2$long, breaks=c(seq(-46.52,-46.01,.01)), right=FALSE) # bin by .01

# get rid of extra columns
PA.test.df2$long <- NULL
PA.test.df2$ID <- NULL

# sum each species column by bins
PA.bin <- aggregate(.~long_bin, PA.test.df2, sum)

# convert to P/A
PA.bin[PA.bin>0] <- 1

# clean up columns, transpose so sites are columns
rownames(PA.bin) <- PA.bin$long_bin
PA.bin$long_bin<- NULL
PA.bin.2 <- t(PA.bin)
PA.bin.2.df <- as.data.frame(PA.bin.2)

# if binning by .01, add in empty bins, including one at beginning for graphing ease (see list of empty bins below)
# [-46.47,-46.46) [-46.44,-46.43);  [-46.41,-46.4) [-46.37,-46.36);  [-46.33,-46.32) [-46.31,-46.3); [-46.28,-46.27) [-46.25,-46.24); [-46.17,-46.16) [-46.15,-46.14); [-46.14,-46.13) [-46.06,-46.05)
PA.bin.2.df <- cbind(x1 = 0, PA.bin.2.df)
target <- which(names(PA.bin.2.df) == '[-46.47,-46.46)')[1]
PA.bin.2.df <- cbind(PA.bin.2.df[,1:target,drop=F], data.frame("x2"=0), PA.bin.2.df[,(target+1):length(PA.bin.2.df),drop=F])
target <- which(names(PA.bin.2.df) == '[-46.41,-46.4)')[1]
PA.bin.2.df <- cbind(PA.bin.2.df[,1:target,drop=F], data.frame("x3"=0), PA.bin.2.df[,(target+1):length(PA.bin.2.df),drop=F])
target <- which(names(PA.bin.2.df) == '[-46.33,-46.32)')[1]
PA.bin.2.df <- cbind(PA.bin.2.df[,1:target,drop=F], data.frame("x4"=0), PA.bin.2.df[,(target+1):length(PA.bin.2.df),drop=F])
target <- which(names(PA.bin.2.df) == '[-46.28,-46.27)')[1]
PA.bin.2.df <- cbind(PA.bin.2.df[,1:target,drop=F], data.frame("x5"=0), PA.bin.2.df[,(target+1):length(PA.bin.2.df),drop=F])
target <- which(names(PA.bin.2.df) == '[-46.17,-46.16)')[1]
PA.bin.2.df <- cbind(PA.bin.2.df[,1:target,drop=F], data.frame("x6"=0), PA.bin.2.df[,(target+1):length(PA.bin.2.df),drop=F])
target <- which(names(PA.bin.2.df) == '[-46.14,-46.13)')[1]
PA.bin.2.df <- cbind(PA.bin.2.df[,1:target,drop=F], data.frame("x7"=0), PA.bin.2.df[,(target+1):length(PA.bin.2.df),drop=F])

# subset PA.bin.2 (sites=col, spec=rows) by FS, NFS
PA.bin.2.df$code <- rownames(PA.bin.2.df)
PA.bin.order.df <- merge(PA.bin.2.df, hab.assoc, by=c("code"))
PA.FS <- subset(PA.bin.order.df, PA.bin.order.df$Classification=="FS")
PA.NFS <- subset(PA.bin.order.df, PA.bin.order.df$Classification=="NFS")

# weight each species by long pref
# list of long in order of % FC
long <- order.model.df.final$long # 91 sites

# presence/absence matrix with sites in % FC order
PA.order.df$code <- NULL
PA.spec.site.df <- t(PA.order.df) # 91 sites

weight.avg(PA.spec.site.df, long)

# append long prefence data to species names
long.pref <- NULL
long.pref$code <- colnames(PA.spec.site.df)
long.pref$weights <- weight.avg(PA.spec.site.df, long)
long.pref.df <- as.data.frame(long.pref)

# merge weights with each subset
PA.FS.weight <- merge(PA.FS, long.pref.df, by=c("code") )
PA.NFS.weight <- merge(PA.NFS, long.pref.df, by=c("code") )

# order each subset by weight
PA.FS.order <- PA.FS.weight[order(-PA.FS.weight$weights),]
PA.NFS.order <- PA.NFS.weight[order(-PA.NFS.weight$weights),]

# while subsetted apart, change all 1s to 2s in PA.NFS so plot in diff colors 
PA.NFS.order[PA.NFS.order==1] <- 2

# rbind together, clean up column names
PA.oo.bin.df <- rbind(PA.FS.order, PA.NFS.order)
rownames(PA.oo.bin.df) <- PA.oo.bin.df$code

# attach column with endemic to MA
MA.end <- hab.assoc.2[,-(1:5),drop=FALSE]
MA.end <- MA.end[,-(2:3), drop=FALSE]
MA.end <- MA.end[,-(3:24), drop=FALSE]
MA.end$end <- ifelse(MA.end$Biogeo_assoc=="MA", 3, 4)
# order in same order as species in dataframe
MA.end.subs <- subset(MA.end, MA.end$Code %in% PA.oo.bin.df$code)
MA.end.subs.order <- MA.end.subs[match(PA.oo.bin.df$code, MA.end.subs$Code),]

PA.oo.bin.df$end <- MA.end.subs.order$end
PA.oo.bin.df$end2 <- MA.end.subs.order$end # do it twice so it shows up wider in figure

# clean up extra columns
PA.oo.bin.df$code <- NULL
PA.oo.bin.df$Classification <- NULL
PA.oo.bin.df$weights <- NULL

# flip columns and rows so plots correctly
PA.oo.bin.flip <- t(PA.oo.bin.df)
PA.oo.bin.m <- as.matrix(PA.oo.bin.flip)

# reverse order of columns
PA.oo.bin.m <- PA.oo.bin.m[,ncol(PA.oo.bin.m):1]

# write a PDF (note that axes and line breaks were customized)
filename <- paste(figures, "Matrix_long_01.pdf", sep = "/")
pdf(file = filename, width = 8, height = 6)

par(mar = c(4,5,5,3))
image(PA.oo.bin.m, col= c("white","dodgerblue", "brown1", "black", "gray"), axes=FALSE, main="Longitude", ylab="Species")
axis(1, at=c(0.0128, .1754, .2914, .4306, .5466, .7786, .8482, .965), labels = c(-46.52, -46.44, -46.37, -46.31, -46.25, -46.15, -46.06, -46.01)) 
axis.break(1,.1638,style="zigzag") 
axis.break(1,.2798,style="zigzag")
axis.break(1,.419,style="zigzag")
axis.break(1,.535,style="zigzag")
axis.break(1,.767,style="zigzag")
axis.break(1,.8366,style="zigzag")

dev.off()


###########################################################################
# Bar plots of each metric vs. FC and Long, with NFS and FS separated
###########################################################################

# bar plots for each metric vs. FC and Long, on the point scale
# NFS and FS plotted in the same window, FS on top with the Y axis flipped
# Note: weird spacing in y-axis labels is to make plots line up in master figure

# replace NAs in model.df.final with 0s for plotting
model.df.final[is.na(model.df.final)] <- 0

# set plotting theme
theme_set(theme_bw(base_size = 19))

#N vs. FC
ptest.1 <- ggplot(model.df.final, aes(FC_LS, TD_NFS_abund)) +
  geom_bar(stat = "identity", color = "brown1" ) +
  theme(axis.line.y=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(), axis.title.y=element_blank(), legend.position="none",
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank()) +
  ylim(0,150) + xlim(20,100) 

ptest.2 <- ggplot(model.df.final, aes(FC_LS, TD_FS_abund)) +
  geom_bar(stat = "identity", color = "dodgerblue" ) +
  scale_y_reverse(lim= c(60, 0), breaks = c(60, 30, 0), labels = c(60, 30, 0)) + # reverses y axis
  xlim(20,100) +
  theme(axis.line=element_blank(),axis.text.x=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(), axis.title.y=element_blank(), legend.position="none",
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank())   

pp1 <- arrangeGrob(rbind(ggplotGrob(ptest.2), ggplotGrob(ptest.1), size = "last")) # stack two plots on top of each other, lining up common x axis


#TD vs. FC
ptest.3 <- ggplot(model.df.final, aes(FC_LS, TD_NFS_raw)) +
  geom_bar(stat = "identity", color = "brown1" ) +
  theme(axis.line=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(), axis.title.y=element_blank(), legend.position="none",
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank()) + xlim(20,100) +
  scale_y_continuous(breaks=c(0, 4, 8), labels = c("    0 ","    4 ", "    8 "))

ptest.4 <- ggplot(model.df.final, aes(FC_LS, TD_FS_raw)) +
  geom_bar(stat = "identity", color = "dodgerblue" ) +
  scale_y_reverse(breaks=c(8, 4, 0), labels = c("    8 ","    4 ", "    0 ")) +
  theme(axis.line=element_blank(),axis.text.x=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(), axis.title.y=element_blank(), legend.position="none",
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank()) + xlim(20,100) 

pp2 <- arrangeGrob(rbind(ggplotGrob(ptest.4), ggplotGrob(ptest.3), size = "last"))


#FD vs. FC
ptest.5 <- ggplot(model.df.final, aes(FC_LS, FD_NFS)) +
  geom_bar(stat = "identity", color = "brown1" ) +
  theme(axis.line=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(), axis.title.y=element_blank(), legend.position="none",
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank()) + xlim(20,100) +
  scale_y_continuous(breaks=c(0, .5, 1), labels = c(" 0"," 0.5", " 1"))

ptest.6 <- ggplot(model.df.final, aes(FC_LS, FD_FS)) +
  geom_bar(stat = "identity", color = "dodgerblue" ) +
  scale_y_reverse(breaks=c(1, .5, 0), labels = c(" 1"," 0.5", " 0")) +
  theme(axis.line=element_blank(),axis.text.x=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(), axis.title.y=element_blank(), legend.position="none",
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank()) + xlim(20,100)

pp3 <- arrangeGrob(rbind(ggplotGrob(ptest.6), ggplotGrob(ptest.5), size = "last"))


#PD vs. FC
ptest.7 <- ggplot(model.df.final, aes(FC_LS, PD_NFS)) +
  geom_bar(stat = "identity", color = "brown1" ) +
  theme(axis.line=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(), axis.title.y=element_blank(), legend.position="none",
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank()) + xlim(20,100) +
  scale_y_continuous(breaks=c(0, .5, 1), labels = c(" 0"," 0.5", " 1"))

ptest.8 <- ggplot(model.df.final, aes(FC_LS, PD_FS)) +
  geom_bar(stat = "identity", color = "dodgerblue" ) +
  scale_y_reverse(breaks=c(1, .5, 0), labels = c(" 1"," 0.5", " 0")) +
  theme(axis.line=element_blank(),axis.text.x=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(), axis.title.y=element_blank(), legend.position="none",
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank()) + xlim(20,100)

pp4 <- arrangeGrob(rbind(ggplotGrob(ptest.8), ggplotGrob(ptest.7), size = "last"))


# Longitude
#N vs. long
ptest.9 <- ggplot(model.df.final, aes(long, TD_NFS_abund)) +
  geom_bar(stat = "identity", color = "brown1" ) +
  theme(axis.line=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(), axis.title.y=element_blank(), legend.position="none",
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank()) +
  ylim(0,150)

ptest.10 <- ggplot(model.df.final, aes(long, TD_FS_abund)) +
  geom_bar(stat = "identity", color = "dodgerblue" ) +
  scale_y_reverse(lim= c(60, 0), breaks = c(60,30, 0), labels = c(60, 30, 0)) + 
  theme(axis.line=element_blank(),axis.text.x=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(), axis.title.y=element_blank(), legend.position="none",
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank()) 

pp5 <- arrangeGrob(rbind(ggplotGrob(ptest.10), ggplotGrob(ptest.9), size = "last"))


#TD vs. long
ptest.11 <- ggplot(model.df.final, aes(long, TD_NFS_raw)) +
  geom_bar(stat = "identity", color = "brown1" ) +
  theme(axis.line=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(), axis.title.y=element_blank(), legend.position="none",
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank()) +
  scale_y_continuous(breaks=c(0, 4, 8), labels = c("    0 ","    4 ", "    8 "))

ptest.12 <- ggplot(model.df.final, aes(long, TD_FS_raw)) +
  geom_bar(stat = "identity", color = "dodgerblue" ) +
  scale_y_reverse(breaks=c(8, 4, 0), labels = c("    8 ","    4 ", "    0 ")) +
  theme(axis.line=element_blank(),axis.text.x=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(), axis.title.y=element_blank(), legend.position="none",
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank())

pp6 <- arrangeGrob(rbind(ggplotGrob(ptest.12), ggplotGrob(ptest.11), size = "last"))


#FD vs. long
ptest.13 <- ggplot(model.df.final, aes(long, FD_NFS)) +
  geom_bar(stat = "identity", color = "brown1" ) +
  theme(axis.line=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(), axis.title.y=element_blank(), legend.position="none",
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank()) +
  scale_y_continuous(breaks=c(0, .5, 1), labels = c(" 0"," 0.5", " 1"))

ptest.14 <- ggplot(model.df.final, aes(long, FD_FS)) +
  geom_bar(stat = "identity", color = "dodgerblue" ) +
  scale_y_reverse(breaks=c(1, .5, 0), labels = c(" 1"," 0.5", " 0")) +
  theme(axis.line=element_blank(),axis.text.x=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(), axis.title.y=element_blank(), legend.position="none",
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank()) 

pp7 <- arrangeGrob(rbind(ggplotGrob(ptest.14), ggplotGrob(ptest.13), size = "last"))


#PD vs. long
ptest.15 <- ggplot(model.df.final, aes(long, PD_NFS)) +
  geom_bar(stat = "identity", color = "brown1" ) +
  theme(axis.line=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(), axis.title.y=element_blank(), legend.position="none",
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank()) +
  scale_y_continuous(breaks=c(0, .5, 1), labels = c(" 0"," 0.5", " 1"))

ptest.16 <- ggplot(model.df.final, aes(long, PD_FS)) +
  geom_bar(stat = "identity", color = "dodgerblue" ) +
  scale_y_reverse(breaks=c(1, .5, 0), labels = c(" 1"," 0.5", " 0")) +
  theme(axis.line=element_blank(),axis.text.x=element_blank(),axis.ticks=element_blank(),
        axis.title.x=element_blank(), axis.title.y=element_blank(), legend.position="none",
        panel.background=element_blank(),panel.border=element_blank(),panel.grid.major=element_blank(),
        panel.grid.minor=element_blank(),plot.background=element_blank()) 

pp8 <- arrangeGrob(rbind(ggplotGrob(ptest.16), ggplotGrob(ptest.15), size = "last"))

# take a look at all of them stacked
grid.arrange(pp1, pp5, pp2, pp6, pp3, pp7, pp4, pp8, ncol=2)

############################################################
# Combine matrices and metric plots into master figure
#############################################################

# Goal-- 10 panel plot, with 1 column for FC, 1 for long
# each column has a matrix, then plots of the variable vs. N, TD, FD, and PD

# Problem--need to combine base graphics and ggplot graphics
# Can do this by converting base graphics to Grobs, the ggplot output

# clear image space so that image conversion doesn't mess up
dev.off()
# Convert FC matrix to grob
par(mar = c(4,1.3,1,0))
image(PA.oo.bin.m.FC, col= c("white","dodgerblue", "brown1", "black", "gray"), axes=FALSE)
axis(1, at=c(0.024, .262, .5, .738, .976), labels = c(20,40,60,80,100), cex.axis = 1.25, col.ticks = "white") 
g1 <- grab_grob()

# clear image space
dev.off()
# Convert long matrix to grob
par(mar = c(4,1.6,1,0))
image(PA.oo.bin.m, col= c("white","dodgerblue", "brown1", "black", "gray"), axes=FALSE)
axis(1, at=c(0.0128, .20324, .39368, .58412, .77456, .965), labels = c(-46.5, -46.4, -46.3, -46.2, -46.1, "-46.0"), cex.axis = 1.25, col.ticks = "white") 
axis.break(1,.1638,style="zigzag") 
axis.break(1,.2798,style="zigzag")
axis.break(1,.419,style="zigzag")
axis.break(1,.535,style="zigzag")
axis.break(1,.767,style="zigzag")
axis.break(1,.8366,style="zigzag")
g2 <- grab_grob()

# Make a blank grob
blank<-rectGrob(gp=gpar(col="white")) 

# write all text/headings
t1 <- textGrob("Forest Cover (%)", gp = gpar(fontsize = 24, fontface = "bold"))
t2 <- textGrob("Longitude", gp = gpar(fontsize = 24, fontface = "bold"), x = .55)

t3 <- textGrob("Species", gp = gpar(fontsize = 19), rot = 90)
t4 <- textGrob("N", gp = gpar(fontsize = 19), rot = 90)
t5 <- textGrob("TD", gp = gpar(fontsize = 19), rot = 90)
t6 <- textGrob("FD", gp = gpar(fontsize = 19), rot = 90)
t7 <- textGrob("PD", gp = gpar(fontsize = 19), rot = 90)

t8 <-  textGrob("A", gp = gpar(fontsize = 19, fontface = "bold"), x = .8, y = .95)
t9 <-  textGrob("C", gp = gpar(fontsize = 19, fontface = "bold"), x = .8, y = .95)
t10 <- textGrob("E", gp = gpar(fontsize = 19, fontface = "bold"), x = .8, y = .95)
t11 <- textGrob("G", gp = gpar(fontsize = 19, fontface = "bold"), x = .8, y = .95)
t12 <- textGrob("I", gp = gpar(fontsize = 19, fontface = "bold"), x = .8, y = .95)
t13 <- textGrob("B", gp = gpar(fontsize = 19, fontface = "bold"), x = .8, y = .95)
t14 <- textGrob("D", gp = gpar(fontsize = 19, fontface = "bold"), x = .8, y = .95)
t15 <- textGrob("F", gp = gpar(fontsize = 19, fontface = "bold"), x = .8, y = .95)
t16 <- textGrob("H", gp = gpar(fontsize = 19, fontface = "bold"), x = .8, y = .95)
t17 <- textGrob("J", gp = gpar(fontsize = 19, fontface = "bold"), x = .8, y = .95)


# draw master figure, save as pdf
filename <- paste(figures, "Master_fig.pdf", sep = "/") 
pdf(file = filename, width = 12, height = 20)

# arrange all grobs into a matrix so that cell widths can be controlled
# blank grob fills in empty space
# text is also a grob
# the numbers in the layout matrix correspond with the order of the elements in grid.arrange
grid.arrange(
  g1, pp1, pp2, pp3, pp4,
  g2, pp5, pp6, pp7, pp8,
  blank,
  t1, t2,
  t3, t4, t5, t6, t7,
  t8, t9, t10, t11, t12, t13, t14, t15, t16, t17,
  layout_matrix = cbind(c(14, 15, 16, 17, 18, 11), c(19, 20, 21, 22, 23, 11), c(1,2,3,4,5,12), c(11,11,11,11,11,11), c(24, 25, 26, 27, 28, 11), c(6,7,8,9,10,13), c(11,11,11,11,11,11)),
  widths = c(.04, .02, .46, .10, .02,.46, .02),
  heights = c(.19, .19, .19, .19, .19, .05)
)

dev.off()

# draw master figure, save as tiff-----this isn't working :( 
filename <- paste(figures, "Master_fig.tiff", sep = "/")
tiff(filename = filename, width = 6.5, height = 10, units = "in", res = 300)

grid.arrange(
  g1, pp1, pp2, pp3, pp4,
  g2, pp5, pp6, pp7, pp8,
  blank,
  t1, t2,
  t3, t4, t5, t6, t7,
  t8, t9, t10, t11, t12, t13, t14, t15, t16, t17,
  layout_matrix = cbind(c(14, 15, 16, 17, 18, 11), c(19, 20, 21, 22, 23, 11), c(1,2,3,4,5,12), c(11,11,11,11,11,11), c(24, 25, 26, 27, 28, 11), c(6,7,8,9,10,13), c(11,11,11,11,11,11)),
  widths = c(.04, .02, .45, .12, .02,.45, .02),
  heights = c(.19, .19, .19, .19, .19, .05)
)

dev.off()

###################################################
# investigate high abundance with high FC for NFS  # get rid of this later?
###############################################

# outlier N points NFS
subset(model.df.final.NFS, model.df.final.NFS$TD_NFS_abund>100)
# 263-4
# 399-18

# let's look at what species are present at those points
# this df has number of each NFS species summed at each point
sum.df.group.2 = ddply(df.group.2, .(Landscape, Habitat, Point, code), summarize, 
                       N = sum(N.trap))
# subset for points of interest
p263.4 <- subset(sum.df.group.2, sum.df.group.2$Landscape==263 & sum.df.group.2$Point==4) # 65 Canthon aff. angularis, 35 Eurysternus cyanescens
p399.18 <- subset(sum.df.group.2, sum.df.group.2$Landscape==399 & sum.df.group.2$Point==18) # 111 Eurysternus parallelus

# Canthon aff. angularis--no
# E cyanescens--yes, borderline + FC relationship
# E parallelus--no


# from NFS spec vs. abundance plots, these spec have + relationship with FC:

# E cyanescens--borderline
# Dichotomius aff. carbonarius sp. 2--borderline

########################################
# plot FC vs. abundance for each species
########################################

# this section only runs when FD section run

# merge df with captures of each species at each point with FC, sep for FS and NFS
# sum.df.group.1  FS
# sum.df.group.2  NFS
# buffer.forest.200  FC data
FS.spec.abund <- merge(sum.df.group.1, buffer.forest.200, by = c("Landscape", "Point")) 
NFS.spec.abund <- merge(sum.df.group.2, buffer.forest.200, by = c("Landscape", "Point"))

# only keep species used in final model
FS.species <- unique(FS.spec.abund$code)
NFS.species.pre <- unique(NFS.spec.abund$code) # 21 species, only 18 in NFS.spec.n.2
NFS.species <- NFS.species.pre[NFS.species.pre %in% NFS.spec.n.2]


filename <- paste(figures, "FS_N_v_FC.pdf", sep = "/")
pdf(file = filename, width = 10, height = 20)

par(mfrow=c(6, 3))

for ( i in FS.species){
  df <- subset(FS.spec.abund, FS.spec.abund$code==i)
  df <- subset(df, df$N>0)
  Nmax <- max(df$N)
  plot(N ~ FC_LS, data = df, pch = 19, xlab="", ylab="", ylim = c(0, Nmax + 3))
  mtext(i, side = 3, cex = 1.5, font = 2)
  if(i != "Coprophanaeus.bellicosus"){
    lm = lm(N ~ FC_LS, data = df)
    abline(coef(lm))
    p <- summary(lm)$coefficients[2,4]
    r <- summary(lm)$r.squared
    mtext(paste("p = ", round(p, 3), "r-squared = ", round(r, 3), ""), side = 1, line = 3)
  }
}

dev.off()



filename <- paste(figures, "NFS_N_v_FC.pdf", sep = "/")
pdf(file = filename, width = 10, height = 20)

par(mfrow=c(6, 3))

for ( i in NFS.species){
  df <- subset(NFS.spec.abund, NFS.spec.abund$code==i)
  df <- subset(df, df$N>0)
  Nmax <- max(df$N)
  plot(N ~ FC_LS, data = df, pch = 19, xlab="", ylab="", ylim = c(0, Nmax + 3))
  mtext(i, side = 3, cex = 1.5, font = 2)
  if(i != "Canthon.septemmaculatus" & i != "Dichotomius.bos" & i != "Canthon.podagricus"){
    lm = lm(N ~ FC_LS, data = df)
    abline(coef(lm))
    p <- summary(lm)$coefficients[2,4]
    r <- summary(lm)$r.squared
    mtext(paste("p = ", round(p, 3), "r-squared = ", round(r, 3), ""), side = 1, line = 3)
  }
}

dev.off()

###############################
# Individual species N vs. long
############################

FS.spec.abund.2 <- merge(FS.spec.abund, point.coords, by=c("Landscape", "Point"))
NFS.spec.abund.2 <- merge(NFS.spec.abund, point.coords, by=c("Landscape", "Point"))

filename <- paste(figures, "FS_N_v_Y.pdf", sep = "/")
pdf(file = filename, width = 10, height = 20)

par(mfrow=c(6, 3))

for ( i in FS.species){
  df <- subset(FS.spec.abund.2, FS.spec.abund.2$code==i)
  df <- subset(df, df$N>0)
  Nmax <- max(df$N)
  plot(N ~ long, data = df, pch = 19, xlab="", ylab="", ylim = c(0, Nmax + 3))
  mtext(i, side = 3, cex = 1.5, font = 2)
  if(i != "Coprophanaeus.bellicosus"){
    lm = lm(N ~ long, data = df)
    abline(coef(lm))
    p <- summary(lm)$coefficients[2,4]
    r <- summary(lm)$r.squared
    mtext(paste("p = ", round(p, 3), "r-squared = ", round(r, 3), ""), side = 1, line = 3)
  }
}

dev.off()



filename <- paste(figures, "NFS_N_v_Y.pdf", sep = "/")
pdf(file = filename, width = 10, height = 20)

par(mfrow=c(6, 3))

for ( i in NFS.species){
  df <- subset(NFS.spec.abund.2, NFS.spec.abund.2$code==i)
  df <- subset(df, df$N>0)
  Nmax <- max(df$N)
  plot(N ~ long, data = df, pch = 19, xlab="", ylab="", ylim = c(0, Nmax + 3))
  mtext(i, side = 3, cex = 1.5, font = 2)
  if(i != "Canthon.septemmaculatus" & i != "Dichotomius.bos" & i != "Canthon.podagricus"){
    lm = lm(N ~ long, data = df)
    abline(coef(lm))
    p <- summary(lm)$coefficients[2,4]
    r <- summary(lm)$r.squared
    mtext(paste("p = ", round(p, 3), "r-squared = ", round(r, 3), ""), side = 1, line = 3)
  }
}

dev.off()


###############################################
# SOM Table 1. Species captured, hab use, traits
###############################################

# this section only runs when FD section run

# Merge df with species used captures, habitat use
SOM.1.1 <- merge(spec.keep.sum.1, hab.assoc, by = c("code"), all = FALSE) 

# Merge that with traits
trait.temp$code <- rownames(trait.temp)
SOM.1.2 <- merge(SOM.1.1, trait.temp, by = c("code"), all = FALSE)

# Sort alphabetically 
SOM.1.3 <- SOM.1.2[order(SOM.1.2$code),]

# replace code with clean names
SOM.1.3$spec <- c("Chalcocopris hesperus", "Canthon ibarragrassoi", "Canthon amabilis", "Canthon aff. angularis", 
                  "Canthon virens chalybaeus", "Canthon septemmaculatus", "Canthon aff. luctuosus", "Canthon podagricus", 
                  "Copropahaneus bellicosus", "Coprophanaeus cerberus", "Coprophanaes saphirinus", 
                  "Deltochilum brasiliense", "Deltochilum dentipes", "Deltochilum furcatum", "Deltochilum morbillosum",
                  "Deltochilum rubripenne", "Dichotomius assifer", "Dichotomius bos", "Dichotomius aff. carbonarius", 
                  "Dichotomius depressicollis", "Dichotomius aff. fissus", "Dichotomius mormon", "Dichotomius quadrinodosus",
                  "Dichotomius aff. carbonarius sp.2", "Eurysternus cyanescens", "Eurysternus francinae", "Eurysternus hirtellus",
                  "Eurysternus infelxus", "Eurysternus parallelus", "Ontherus sulcator", "Pseudocanthon aff. felix", "Phanaeus dejeani", 
                  "Phanaeus splendidulus", "Scybalocanthon nigriceps", "Sylvicanthon aff. foveiventris")

# Redo alphabetical sorting
SOM.1.4 <- SOM.1.3[order(SOM.1.3$spec),]

# Clean up extra columns, order, rename
SOM.1.4$code <- NULL

table <- SOM.1.4[moveme(names(SOM.1.4), "spec first")]

aic.names = c("Species", "N", "Classification","Activity period", "Nesting strategy","Body mass (mg)")
colnames(table) <- aic.names

# Print table
filename <- paste(results, "SOM_Table_1.doc", sep = "/")
AIC_table <- RTF(file = filename)
setFontSize(AIC_table, 12)
addParagraph(AIC_table, "SOM Table 1. All species captured, their abundances, habitat uses, and traits. N stands for the total number of individuals of each species captured. Bodymass was calculated as average bodymass in mg of species measured. In total, 35 species from 11 genera were captured. Of those, 17 were forest specialists and 18 were non-forest specialists. (Abbreviations: habitat use, NFS = non-forest specialist, FS = forest specialist; activity period, D = diurnal and N = nocturnal; guild, R = roller (telocoprid), B = burrower (paracoprid), D = dweller (endocoprid)) \n")
addTable(AIC_table, as.data.frame(table))
done(AIC_table)

########################################
# calculating phylogenetic signal
#######################################

# note that this section only works if you run the section of code that makes the phylogeny

# we need a file with traits and names--trait.temp 
# we need a phylogeny--species.phy
# we need codes to match between them

# prune trait db to get rid of species not studied
spec <- species.phy$tip.label
trait.temp2 <- subset(trait.temp, trait.temp$code %in% spec)

# order in same order as tiplabels
trait.temp3 <- trait.temp2[species.phy$tip.label, ]

# Use Lari's code to calculate phylogenetic signal for categorical traits:

#from https://stat.ethz.ch/pipermail/r-sig-phylo/2011-March/001037.html

# This function tests for phylogenetic signal with categorical traits.
# It works similarly to phylo.signal, by randomizing the tip data and comparing the minimum number of character state changes with a null model.
# Minimum character state change is obtained with parsimony, and the syntax allows for different evolutionary models.
# To build a matrix of costs of character state transition, see Maddison & Maddison 2000. MacClade 4 Manual pp.69-72 (unordered parsimony is the default, cost=NULL).
# Note that this function corresponds to the "Fixed Tree, Character Randomly Reshuffled" model proposed in Maddison & Slatkin (1991) Evolution 45:1184.

'phylo.signal.disc' <- function (trait, phy, rep = 999, cost = NULL) 
{lev <- attributes(factor(trait))$levels
if (length(lev) == length(trait))
  stop("Sure this trait is categorical?")
if (is.null(cost)){
  cost1 <- 1 - diag (length(lev))
}
else {
  if (length(lev) != dim(cost) [1] )
    stop("Dimensions of the character state transition matrix 
         do not aggre with the number of levels")
  cost1 <- t(cost)
}
dimnames(cost1) <- list(lev, lev)
trait <-as.numeric(trait)
attributes(trait)$names <- phy$tip
NULL.MODEL <- matrix(NA, rep, 1)
obs <-  t ( data.frame ( trait ) )
obs <- phyDat (t (obs), type = "USER", 
               levels = attributes(factor(obs))$levels)
OBS <- parsimony(phy, obs)
for (i in 1:rep){
  null <- sample (as.numeric(trait))
  attributes(null)$names <- attributes(trait)$names
  null <- t(data.frame(null))
  null <-
    phyDat (t(null) , type = "USER", 
            levels = attributes(factor(null))$levels)
  NULL.MODEL[i,] <- parsimony(phy, null)
  P.value <- sum(OBS >=NULL.MODEL) / (rep + 1)
}
par(mfrow = c(1,2))
hist(NULL.MODEL, xlab = "Transitions.in.Randomizations", 
     xlim = c(min(c(min(NULL.MODEL, OBS-1))), max(NULL.MODEL)+1))
arrows(OBS, rep/10, OBS, 0, angle=20, col = "red", lwd = 4)
phy$tip.label <- rep (".", length(trait))
plot(phy, tip.col = trait + 10, cex = 250/length(trait), 
     font = 1)
title("Character states")
par(mfrow = c(1,1))
OUTPUT1 <- t(data.frame(Number.of.Levels = length(attributes				(factor(trait))$levels), 
                        Evolutionary.Transitions.Observed=OBS, Evolutionary.Transitions.Randomization.Median=median(NULL.MODEL), Evolutionary.Transitions.Randomization.Min=min(NULL.MODEL), Evolutionary.Transitions.Randomization.Max=max(NULL.MODEL),P.value))
if(is.null(cost)) {
  list(.Randomization.Results = OUTPUT1,.Levels = 
         lev,.Costs.of.character.state.transition.UNORDERED.PARSIMONY  =  
         t(cost1))
}
else {
  list(.Randomization.Results=OUTPUT1,.Levels= lev,.Costs.of.character.state.transition.FROM.ROW.TO.COL = t(cost1))
}
}


result<-phylo.signal.disc(trait.temp3$Diel,species.phy,cost=NULL)
summary(result) # p = .999

result<-phylo.signal.disc(trait.temp3$Nest,species.phy,cost=NULL)
summary(result) # p = .999


# Use bloomberg's k to calculate phylogenetic signal for continuous traits

bm <- trait.temp3$BM_mean_mg

result <- phylosig(species.phy, bm, method = "K", test = TRUE, nsim = 999) #species with NA (1) had to be dropped from tree
# $K
# [1] 0.03404962
# 
# $P
# [1] 0.1441441

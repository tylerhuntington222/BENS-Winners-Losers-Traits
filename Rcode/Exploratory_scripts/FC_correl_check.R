##### Check correlation between averaged LC and reported for landscape 3km transects

# load packages

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

# set working directory

setwd("/Users/Tyler/Dropbox/Interface/Dados/Non-Geospatial/Winners_Losers_Traits/Analysis/Data")

list.files()
# load data

pt_stats.df <- read.csv("Landscape_composition_pitfalls.csv" , header = T)

ls_stats.df <- read.xlsx("LS_summary_stats.xls", sheetIndex = 3)

ls_stats.df <- subset (ls_stats.df, select = c("Landscape" , "X3km_FC._new"))

ls_stats.df <- subset (ls_stats.df, ls_stats.df$Landscape != 220)

# subset pt_stats for 3km buffers
pt_stats.df <- subset (pt_stats.df, pt_stats.df$Buffer_m == 200)

pt_stats.df <- subset (pt_stats.df, pt_stats.df$Class == "forest")
summary (as.factor(pt_stats.df$Buffer_m))

head (pt_stats.df)

levels (pt_stats.df$Class)

hist(pt_stats.df$Per_area)

# subset relvant columns of pt_stats.df 
pt_stats.df = subset (pt_stats.df, select = c("Landscape", "Point", "Code_min" , "Per_area", "Trap_num"))

# summarize pt level data at LS level
pt_stats.df = ddply(.data = pt_stats.df, c("Landscape", "Point", "Code_min"), .fun = summarize, pt_meanFC = mean(Per_area))
pt_stats_av.df <- ddply(.data = pt_stats.df, c("Landscape"), .fun = summarize, ls_meanFC = mean(pt_meanFC))

pt_stats_av.df


sort(ls_stats.df$Landscape)
sort(pt_stats_av.df$Landscape)

sort(ls_stats.df$X3km_FC._new)
sort(pt_stats_av.df$ls_meanFC)

hist(ls_stats.df$X3km_FC._new, breaks = 6)
hist(pt_stats_av.df$ls_meanFC, breaks = 6)

# plot relationship
plot (sort(pt_stats_av.df$ls_meanFC), sort(ls_stats.df$X3km_FC._new))
abline(0,1)


pbuff_fc <- sort(pt_stats_av.df$ls_meanFC)
ls_FC <- sort(ls_stats.df$X3km_FC._new)

lm (pbuff_fc ~ ls_FC)
cor(pbuff_fc , ls_FC)

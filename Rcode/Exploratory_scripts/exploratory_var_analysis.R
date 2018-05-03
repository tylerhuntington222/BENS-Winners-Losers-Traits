#####################################################################################################################################
################################## EXPLORATORY ANALYSIS OF VARIABLE DISTRIBUTIONS IN a.df AND b.df ################################## 
#####################################################################################################################################

# set working directory
setwd("/Users/nicholslab/Dropbox/Interface/Dados/Non-Geospatial/Winners_Losers_Traits/Analysis/Data")

# load beetle and avian dataframes
a.df <- readRDS ( "a.df.rds" )
b.df <- readRDS ( "b.df.rds" )

head(b.df)


################################## EXPLORATORY ANALYSIS OF VARIABLE DISTRIBUTIONS IN b.df ################################## 

# create b.df variant with no species/point pairings where no capture occurre
b.df.caps <- subset ( b.df , b.df$Presence == 1)

ls.count.df <-ddply ( b.df.caps, .(Species , Landscape) , 
                         summarise, N.total = mean(N.tot))

ls.count.df <- data.frame (table (ls.count.df$Species , ls.count.df$Landscape))
colnames(ls.count.df) <- c("Species" , "Landscape", "Presence")

# create df with num of landscapes each species was cap'd in
b.sp.ls.count.df <-ddply ( ls.count.df, .(Species) , 
                      summarise, LS.total = sum(Presence))

# merge b.sp.ls.count.df with b.df
b.df <- merge (b.df, b.sp.ls.count.df)

# subset b.df for species that were cap'd in at least three landscapes
b.df <- subset ( b.df , LS.total >= 3 )

# create df with unique entry for each beetle species and associated attributes
b.sp.df <- ddply ( b.df, .(Species , Diel , Nest , Biogeo , BM , Hab , Diet, N.tot), summarize, mean(N.tot))
b.sp.df <- b.sp.df[,1:8]


### EXAMINE VARIABLES FOR NON-STATIONARITY

# take log BMs of beetle species and bin into four weight classes
b.logBM.bins <- .bincode (log(b.sp.df$BM) , breaks = c(0,2,4,6,8), right = T)

# Check contingency table and interaction between Diel & Nest
table (b.sp.df$Diel, b.sp.df$Nest)
chisq.test (b.sp.df$Diel, b.sp.df$Nest)

# Check contingency table and interaction between Diel & Biogeo
table (b.sp.df$Diel, b.sp.df$Biogeo)
chisq.test (b.sp.df$Diel, b.sp.df$Biogeo)

# Check contingency table and interaction between Diel & BM
table (b.sp.df$Diel, b.logBM.bins)
chisq.test (b.sp.df$Diel, b.logBM.bins)

# Check contingency table and interaction between Diel & Hab
table (b.sp.df$Diel, b.sp.df$Hab)
chisq.test (b.sp.df$Diel, b.sp.df$Hab)

# Check contingency table and interaction between Nest & Biogeo
table (b.sp.df$Nest, b.sp.df$Biogeo)
chisq.test (b.sp.df$Nest, b.sp.df$Biogeo)

# Check contingency table and interaction between Nest & Diet
table (b.sp.df$Nest, b.sp.df$Diet)
chisq.test (b.sp.df$Nest, b.sp.df$Diet)

# Check contingency table and interaction between Nest & BM
table (b.sp.df$Diel, b.logBM.bins)
chisq.test (b.sp.df$Nest, b.logBM.bins)

# Check contingency table and interaction between Nest & Hab
table (b.sp.df$Nest, b.sp.df$Hab)
chisq.test (b.sp.df$Nest, b.sp.df$Hab)

# Check contingency table and interaction between Biogeo & Diet
table (b.sp.df$Biogeo, b.sp.df$Diet)
chisq.test (b.sp.df$Biogeo, b.sp.df$Diet)

# Check contingency table and interaction between Biogeo & BM
table (b.sp.df$Biogeo, b.logBM.bins)
chisq.test (b.sp.df$Biogeo, b.logBM.bins )

# Check contingency table and interaction between Biogeo & Hab
table (b.sp.df$Biogeo, b.sp.df$Hab)
chisq.test (b.sp.df$Biogeo, b.sp.df$Hab)

# Check contingency table and interaction between Diet & BM
table (b.sp.df$Diet, b.logBM.bins)
chisq.test (b.sp.df$Diet, b.logBM.bins)

# Check contingency table and interaction between Diet & Hab
table (b.sp.df$Diet, b.sp.df$Hab)
chisq.test (b.sp.df$Diet, b.sp.df$Hab)

# Check contingency table and interaction between BM & Hab
table (b.logBM.bins, b.sp.df$Hab)
chisq.test (b.logBM.bins, b.sp.df$Hab)



#############################

# set plot matrix
par(mfrow=c(3,2))

## Diel
# plot histogram of beetle species by Diel
plot (b.sp.df$Diel, main = "Diel")
summary (b.sp.df$Diel)

## Nest
# plot histogram of beetle species by nesting strategy
plot (b.sp.df$Nest, main = "Nest")
summary (b.sp.df$Nest)

## Diet
# plot histogram of beetle species by diet 
plot (b.sp.df$Diet, main = "Diet")
summary (b.sp.df$Diet)

## Body Mass
# plot histogram of beetle species by logBM
hist (log(b.sp.df$BM), main = "log(BM)")
summary (b.sp.df$BM)

## Biogeo association
# plot histogram of beetle species by Biogeo
plot (b.sp.df$Biogeo, main = "Biogeo")
summary (b.sp.df$Biogeo)

## Habitat preference
# plot histogram of beetle species by Hab Assoc
plot (b.sp.df$Hab, main = "Hab Assoc")
summary (b.sp.df$Hab)




################################################ Plot N ~ FC for all each beetle species

# make plot matrix of individual beetle species N ~ FC
par(mfrow=c(2,3))
par (pty = "s")
par(mar=c(2, 2,  2,  2))

for (i in b.sp.df$Species) {
  
  # create df with entries as all captures of a particular species
  b.par.df <- subset (b.df , b.df$Species == i & b.df$Presence == 1)

  # plot N ~ FC for particular species
  plot(b.par.df$perFC, jitter(b.par.df$N, 5), main = i)

  print (i)
  
  }

################################################ CALCULATE N ~ FC RESPONSES FOR ALL BEETLE SPECIES

# create blank beetle FC response df
b.respFC.df <- data.frame(
                          Species = factor(), 
                          Intercept = double(),
                          Inter.error = double(),
                          Inter.t.val = double(),
                          Inter.p.val = double(), 
                          Slope = double (),
                          Slope.error = double (),
                          Slope.t.val = double(),
                          Slope.p.val = double()
                          ) 
# set levels of species var
levels(b.respFC.df$Species) <- unique (b.sp.df$Species)

# initialize row 1 as first slot of b.respFC.df to be filled
row = 1

# iterate response calculation over all beetle species and store in b.respFC.df
for (i in b.sp.df$Species) {
  
  # create df with entries as all captures of a particular species
  b.par.df <- subset (b.df , b.df$Species == i )
  
  # plot N ~ FC for particular species
  b.lm.par <- lm (b.par.df$N ~ b.par.df$perFC)
  
  print (summary (b.lm.par)$coefficients)
  
  b.respFC.df[row, 1] <- (i)
  b.respFC.df[row,2:5] <- summary (b.lm.par)$coefficients[1,]
  b.respFC.df[row,6:9] <- summary (b.lm.par)$coefficients[2,]
  
  row <- row + 1
}

b.respFC.df[1,1] = "hello"

print (b.respFC.df)
summary (b.respFC.df)

hist (b.respFC.df$Slope, breaks = 20)

# species with significant slope coefficients
b.sigResp <- subset (b.respFC.df, b.respFC.df$Slope.p.val < 0.05)

hist (b.sigResp$Slope, breaks = 15)
hist (b.sigResp$Slope.p.val, breaks = 15)

################################## EXPLORATORY ANALYSIS OF VARIABLE DISTRIBUTIONS IN a.df ################################## 

# create df with unique entry for each beetle species and associated attributes
a.spec.df <-ddply ( a.df, .(Species , Nest , Biogeo , BiogeoPlas , BM , Habitat, HabPlas, BM, Fecund, Diet ) , summarize , N.total = sum ( N ) )


# create df with entries = all captures of a particular species
A.par.df <- subset (b.df , b.df$Species == "Onthophagus.sp1" )

# plot N ~ FC for particular species
plot (b.par.df$N ~ b.par.df$FC)

a.df[a.df$Biogeo == "T"]
head (a.df)
unique (subset (a.df$Species, a.df$Biogeo == "T"))


summary (a.df$Biogeo)
which(a.spec.df$BiogeoPlas == 21)

count(a.spec.df)


hist(subset (a.spec.df$BiogeoPlas,a.spec.df$BiogeoPlas != 1))



levels(b.df$TrapID)

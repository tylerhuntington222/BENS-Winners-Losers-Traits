
######################
# GEE model
######################
comm <- matrix(runif(1896,1), nrow = 1)
colnames(comm) <- beetles$SPCODE4
phy <- prune.sample(comm, exp_phy)
beetles.clean = data.frame(beetles$ref, beetles$smod, beetles$Prov, beetles$logratio, beetles$ndiel, beetles$nguild, beetles$z.biomass)
colnames(beetles.clean) = c("ref", "smod", "prov",  "logratio", "ndiel", "nguild", "z.biomass")
beetles.clean <- data.frame(beetles.clean)
beetles.clean$nmod = factor(beetles.clean$smod)
attach(beetles.clean)

globalPC = compar.gee (logratio ~ nmod 
                       + z.biomass 
                       + ndiel 
                       + nguild 
                       + nmod * z.biomass 
                       + nmod * ndiel
                       + nmod * nguild
                       + ref,
                       phy= phy)				

print(globalPC)
#non PC now
globalNPC = glm (logratio ~ nmod 
                 + z.biomass 
                 + ndiel 
                 + nguild 
                 + nmod * z.biomass 
                 + nmod * ndiel
                 + nmod * nguild
                 + ref)

detach(beetles.clean)				
global.drop =  drop1(globalPC, scope, quiet = FALSE)

global.drop

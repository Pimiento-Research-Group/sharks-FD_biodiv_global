#########################################################################
####                                                                  ###
####        Example of how to use the Fuse function                   ###
####                                                                  ###
####  7/05/2020 Camille Albouy & Fabien leprieur                      ###
#########################################################################

### Set working direectory
setwd("/home/emh-albouy/Documents/Projets/Package_FD/Function_fab/")

### Packages loading
lib_vect<-c("ade4","cluster","geometry","reshape2", "vegan", "ape")
sapply (lib_vect, library,character.only=TRUE) 

### Source function
source("FUSE_function_end.R")

###-------------- Data loading and modification----------------######
### Data loading
sp_tr <- read.csv("data_traits_MMA_ursus.csv", dec=",", sep=";",
                     header=TRUE,row.names=1,na.strings="NA")

### Trait compilation and ordination
dimorphism <- ordered(sp_tr$dimorphism) 
breeding_site   <- ordered(sp_tr$breeding_site)
social_behavior <- ordered(sp_tr$social_behavior)
weight_max <- log(sp_tr$adult_weight_max)
social_group <- log(sp_tr$social_group_mean)

### Trait Matriw construction
sp_tr_end <- data.frame(main_diet=sp_tr$main_diet,foraging_water_depth=sp_tr$foraging_water_depth, 
                  foraging_location=sp_tr$foraging_location,fasting_strategy=sp_tr$fasting_strategy, 
                  female_sexual_maturity=sp_tr$female_sexual_maturity,weaning=sp_tr$weaning, 
                  gestation=sp_tr$gestation,inter_litter=sp_tr$inter_litter,breeding_site=breeding_site, 
                  social_group=social_group,social_behavior=social_behavior,weight_max=weight_max, 
                  dimorphism=dimorphism) 
rownames(sp_tr_end) <- rownames(sp_tr) 

### Function weigthing vector
v <- c(0.25,0.25,0.25,0.25,0.20,0.20,0.20,0.20,0.20,0.5,0.5,0.5,0.5)

### Gower distance calculation
sp_dist_tr <- daisy(sp_tr_end,metric=c("gower"),type=list(symm=c(4)),weights=v)

### Principal coordinate analyses
Pcoa <- dudi.pco(quasieuclid(sp_dist_tr),scann=F,nf=40)
coords <- Pcoa$li[1:40]

### FUSE calculation
FUSE_res<-FUSE(Mat_dist=sp_dist_tr,Coords=coords,nb_NN=5,GE=sp_tr$IUCN_num,StandGE=T)

FUSE_res2<-FUSE(Mat_dist=sp_dist_tr,Coords=coords,nb_NN=5,GE=sp_tr$IUCN_50,StandGE=T)

FUSE_res3<-FUSE(Mat_dist=sp_dist_tr,Coords=coords,nb_NN=5,GE=sp_tr$IUCN_100,StandGE=T)

########################################################################
###-------------- end of example FUSE calculation ------------------###
#######################################################################


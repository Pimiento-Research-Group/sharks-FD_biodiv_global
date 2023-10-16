#########################################################################                                                                                                               #  
###                                                                   ###
###   FUSE : Functionally Unique, Specialized and Endangered          ### 
###     Authors : Fabien Leprieur (fabien.leprieur@umontpellier.fr)   ###
###         and Camille Albouy (albouycamille@gmail.com)              ###
###                    03/06/2020                                     ###
###                                                                   ###
#########################################################################

#########################################################################
###  Description: The function FUSE provides.....                     ###
###                                                                   ###
### Arguments                                                         ###
### Mat_dist: a dist object provided by the daisy() of the cluster    ###
###     package or dist.ktab() of the ade4 package.                   ###
### Coords: a data.frame with the coordinates of the species on a     ###
###     multidimensional space based on a selected number of          ###
###     axes derived from a Principal Coordinate Analysis (PCOA).     ###
###     The species are in rows and the PCOA axes are in column.      ###
### nb_NN: a numercial value giving the number of nearest neighbour   ###
###     to consider, 5 by default                                     ###
### GE: a numerical vector giving the IUCN status rank (DD = NA,      ###
###     LC=0, NT=1, VU=2, EN=3, CR=4) or the IUCN extinction          ###
###     probability associated with each status, see Mooers et al.    ###
###     (2008) for example with DD=NA, LC=0, NT=0.1, VU=0.4,          ###
###     EN=0.666,CR=0.999)      
### StandGE: option to standerdize the GE values (TRUE or FALSE) 
### Outputs : a dataframe with the species in rows and the different  ###
###      metrics in columns                                           ###
###    FUSE: functionally unique, specialized and endangered          ###
### (see Pimiento et al. 2020 and Griffin et. al 2020)                ###
### FDGE: functionally distinctive and endangered                     ###
### FUn_std: functional uniqueness standardized between 0 and 1       ###
###  (see Mouillot et al. 2013 and Griffin et al. 2020)               ###
### FSp_std:functional specialization standardized between 0 and 1    ###
###  (see Mouillot et al. 2013 and Griffin et al. 2020)               ###
### FDist_std:functional distinctivness standardized between 0 and 1  ###
###  (see Violle et al. 2007 and Griffin et al. 2020)                 ###
#########################################################################

FUSE <- function(Mat_dist,Coords,nb_NN=5,GE,StandGE=F){
  
  require(vegan)
  require(reshape2)
  
  if(!identical(row.names(as.matrix(Mat_dist)),row.names(Coords))){
    stop("Coords lines do not match with the distance matrix")
  } # end of if
  
  nm <- rownames(Coords)
  
  # Specialization calculation
  O <- apply(Coords,2,mean)
  spe <- apply(Coords, 1,function(x){sum((x-O)^2)^0.5})
  
  # Distinctivness calculation 
  dist_sp <- as.matrix(Mat_dist)
  Fdistinct <- apply(dist_sp,1,mean)
  
  # Uniqueness calculation
  uni_res <- get_indicator(Mat_dist=as.matrix(Mat_dist),nb_NN=nb_NN)
  uniqu <- uni_res$Average_uniqueness[,"Mean"]
  
  if(StandGE==TRUE){
    GE <- as.vector(decostand(GE,"range",na.rm=T))
  } #end of if
  
  # FUSE metrics calculation
  FUn_std <- as.vector(decostand(uniqu,"range"))
  FUGE <- log(1+(FUn_std*GE))
  FSp_std <- as.vector(decostand(spe,"range")) 
  FSGE <- log(1+(FSp_std*GE))
  FUS <- FUn_std+FSp_std
  FDist_std <- as.vector(decostand(Fdistinct,"range"))
  FUSE <- setNames(FUGE+FSGE,nm=nm)
  FDGE <- log(1+(FDist_std*GE))
  
  data.frame(cbind(FUSE,FUS,FUGE,FDGE,FSGE,FUn_std,FSp_std,FDist_std))
  
} # end of FUSE function

########################################################################
###------------------------ Internal functions ----------------------###
########################################################################

########################################################################
#####----------------  get_indicator function -------------------- #####
####  Internal function that allow the FUSE function to calculate   ####
####         the nearest neighbors for a all considered species     ####
#### nb_NN: a numercial value giving the number of nearest          ####
####                  neighbour to consider                         ####
####                                                                ####
####  Mat_dist: a matrix object representing the distance in        ####
####  an euclidean space between species based on species trais     ####
####                                                                ####
########################################################################

get_indicator <- function(Mat_dist,nb_NN){
  
  w <- melt(Mat_dist)
  s <- split(w,f=w[,2])
  
  Res <- lapply(s,function(x){get_dist_func(nb_NN=nb_NN,data=x)})
  Res_mean_sd <- do.call(rbind,lapply(1:length(Res),function(i){Res[[i]][[1]]}))
  NN <- lapply(1:length(Res),function(i){Res[[i]][[2]]})
  
  rownames(Res_mean_sd) <- names(NN) <- names(Res)
  list(Average_uniqueness=Res_mean_sd,Nearest_neighbour=NN)
  
} # end of get_indicator

########################################################################
#####----------------  get_dist_func function -------------------- #####
####  Internal function that allow the FUSE function to calculate   ####
####         the nearest neighbors for a given species              ####   
#### nb_NN: a numercial value giving the number of nearest          ####
####                  neighbour to consider                         ####
####                                                                ####
#### data: a dataframe considering all the distance betwwen the     ####
####considered sepcies and all of these neighbours                  ####
########################################################################

get_dist_func <- function(nb_NN,data){
  
  data <- data[order(data[,3], decreasing=F),]
  data <- data[-1,]
  mm <- mean(data[1:nb_NN,3])
  sd <- sd(data[1:nb_NN,3])
  sp <- as.character(data[1:nb_NN,1])
  list(c(Mean=mm,Sd=sd),Species=sp)
  
} # end of get_dist_func

#####################################################################
###-------------------  end of functions ------------------------####
#####################################################################

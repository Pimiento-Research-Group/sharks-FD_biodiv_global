################################################################################################################
#                                                                                                              #
#                                                                                                              # 
#                                 Species and assemblage-level phylogenetic metrics                            #  
#                                                                                                              #  
#                                 Pimiento et al. Functional diversity of sharks and rays                      #
#                                 is highly vulnerable and supported by uniquespecies and                      #  
#                                 locations worldwide Submitted to Nature Communications August 2023           #
#                                                                                                              #  
#                                 Author : fabien leprieur (fabien.leprieur@umontpellier.fr)                   #
#                                                                                                              #  
################################################################################################################



# Main directory 

setwd(".../phylo_analysis") # your directory

# 1- Load the libraries 
library(picante)
library(ape)
library(phyloregion)
library(picante)
library(parallel)
library(PhyloMeasures)
library(phylobase)

# 2-  Load the species occurrence dataset 

load("data/grids_commonspecies.Rdata")
grids <- grids.fd
colnames(grids) <- gsub(" ", "_", colnames(grids))
new_mat<-grids # 1015 species

# 3- Load the file including synonyms between species names in the trait dataset and those in the phylogeny 

Syn_nov2021<-readRDS("data/Syn_nov2021")

# 4-  Load the phylogenetic trees  

tree<-read.nexus("data/10.cal.tree.nex") # phylogenetic trees available in https://vertlife.org/data/sharks/

new_tree_sharks<-sample(tree, 100) # 100 trees take at random 

# 5- Change the tip labels to match species names between the species trait and occurrence datasets

list_tree_sharks_p<-list()

for (i in 1 : length(new_tree_sharks)) {

new_tree_sharks_p<-as.phylo(new_tree_sharks[[i]])  

Vect_pos <- sapply(Syn_nov2021[,2], function(x){grep(as.character(x),new_tree_sharks_p$tip.label)})

Vect_pos <- do.call(c,lapply(1:length(Vect_pos),function(i){
  if(sum(Vect_pos[[i]])==0){
    Vect_pos[[i]]<- NA
  } else {
    Vect_pos[[i]]<- Vect_pos[[i]]
  }
}) )

new_tree_sharks_p$tip.label[Vect_pos[which(Vect_pos!=0)]] <- as.character(Syn_nov2021[which(Vect_pos!=0),1])

list_tree_sharks_p[[i]]<-new_tree_sharks_p

}


# 6- Match between the species included in the species occurrence dataset and those in the phylogenetic tree : species in the phylogeny not found in the species occurrence dataset are removed.

data_sharks_100<-list()
for (i in 1:100) {
data_sharks_100[[i]]<-match.phylo.comm(list_tree_sharks_p[[i]], new_mat)}

# 7- 1 Quantification of assemblage-level phylogenetic diversity metrics : phylogenetic diversity and mean nearest taxon distance (or Phylogenetic Uniqueness PUn)

# 7-1-1 : Phylogenetic Diversity sensu Faith (1992)

matrix_sharks<-as.matrix(data_sharks_100[[1]]$comm)
d<-as(matrix_sharks, "sparseMatrix")
PD<-do.call(rbind, lapply(data_sharks_100, function (x) {PD(d, x$phy)}))
Mean_PD_100T <- apply(PD,2,mean) # mean phylogenetic diversity per cell for the 100 trees taken at random 
saveRDS(Mean_PD_100T, "outputs/Mean_PD_100T")

SR_sharks<-apply(data_sharks_100[[1]]$comm, 1, sum)
residuals_PD_sharks<-residuals(lm(Mean_PD_100T~SR_sharks)) # residuals of the relationship between PD and SR

# 7-2-2: Mean Nearest Taxon Distance (phylogenetic uniqueness)

MNTD<-do.call(rbind, mclapply(data_sharks_100, mc.cores=6, function (x) {mntd.query(x$phy, data_sharks_100[[1]]$comm)}))
Mean_MNTD_100T <- apply(MNTD,2,mean)
saveRDS(Mean_MNTD_100T, "outputs/Mean_MNTD_100T")

# 7-2-3: Compiling the assemblage-level phylogenetic metrics 

latlon<-strsplit(rownames(new_mat), "_")
latlon<-do.call(rbind, latlon)
latlon_final<-data.frame(lat=as.numeric(latlon[,2]), lon=as.numeric(latlon[,1]))

sharks_phylo_results_june23<-data.frame(x=latlon_final$lon, y=latlon_final$lat, SR=SR_sharks, PD=Mean_PD_100T, residuals_PD=residuals_PD_sharks, PUn=Mean_MNTD_100T)
saveRDS(sharks_phylo_results_june23, file="outputs/sharks_phylo_results_june23.rds")

# 8 : Quantifying species-level metrics. : EDGE2 (equivalent to HEDGE) and ED2 (equivalent to HED) per species according to 

# 8-1 : Heightened Evolutionary Distinctiveness HED (called ED2 in Gumbs et al. 2023)

source("functions/HEDGE_function.R") # function provided in Robuchon et al. (2021) Revisiting species and areas of interest for conserving global mammalian phylogenetic diversity. Nature communications 12: 3694 
source("functions/checkphyloarg.R") # function provided in the R package adiv : https://rdrr.io/cran/adiv/
sharks_iucn_final<-readRDS("data/sharks_IUCN_final") # IUCN status for each species and the corresponding extinction probabilities for 50 and 100 years given by Mooers et al. (2008)

proba100<-sharks_iucn_final$iucn_prob100 # we selected the extinction probabilities for 100 years
names(proba100)<-rownames(sharks_iucn_final)

HED<-HED2(data_sharks_100[[1]]$phy, proba100, subtree=FALSE, tol=1e-8) # example with the first tree

HED_sharks<-mclapply(data_sharks_100, mc.cores=6, function(x) {HED2(x$phy, proba100, subtree=FALSE, tol=1e-8)}) # we calculate for the 100 random trees

HED_sharks_final<-as.data.frame(do.call(cbind, lapply(HED_sharks,function(y){y$scores$HED})))
rownames(HED_sharks_final)<-data_sharks_100[[1]]$phy$tip.label
saveRDS(HED_sharks_final, file="outputs/HED_sharks_final.rds")

HED_sharks_final_mean<-apply(HED_sharks_final, 1, mean) # mean of HED for the 100 random trees
HED_sharks_final_sd<-apply(HED_sharks_final, 1, sd) # sd of HED for the 100 random trees 
Res_HED_sharks<-data.frame(Species=rownames(HED_sharks_final), HED_mean=HED_sharks_final_mean, HED_sd=HED_sharks_final_sd)
Res_HED_sharks<-Res_HED_sharks[order(rownames(Res_HED_sharks)),]
saveRDS(Res_HED_sharks, file="outputs/Res_HED_sharks")

# 8-2 : Heightened Evolutionary Distinct and Globally Endangered HEDGE (called EDGE2 in Gumbs et al. 2023)

HEDGE_sharks_final<-as.data.frame(do.call(cbind, lapply(HED_sharks,function(y){y$scores$HEDGE})))
rownames(HEDGE_sharks_final)<-data_sharks_100[[1]]$phy$tip.label
saveRDS(HEDGE_sharks_final, file="outputs/HEDGE_sharks_final.rds")

HEDGE_sharks_final_mean<-apply(HEDGE_sharks_final, 1, mean)
HEDGE_sharks_final_sd<-apply(HEDGE_sharks_final, 1, sd)
Res_HEDGE_sharks<-data.frame(Species=rownames(HEDGE_sharks_final), HEDGE_mean=HEDGE_sharks_final_mean, HEDGE_sd=HEDGE_sharks_final_sd)
Res_HEDGE_sharks<-Res_HEDGE_sharks[order(rownames(Res_HEDGE_sharks)),]
saveRDS(Res_HEDGE_sharks, "outputs/Res_HEDGE_sharks_phy")


# 9 : top 25% ED2 and EDGE2 species 

# 9-2:  Top 25% ED2 #######

Top_ED2_Q25<-Res_HED_sharks[which(Res_HED_sharks$HED_mean>quantile(Res_HED_sharks$HED_mean, 0.75)),]

saveRDS(Top_ED2_Q25, file="outputs/Top_ED2_Q25.rds")

# 9-3 : Top 25% EDGE2 #######

Top_EDGE2_Q25<-Res_HEDGE_sharks[which(Res_HEDGE_sharks$HEDGE_mean>quantile(Res_HEDGE_sharks$HEDGE_mean, 0.75)),]

saveRDS(Top_EDGE2_Q25, file="outputs/Top_EDGE2_Q25.rds")



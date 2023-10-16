#####################################################################################################
#                                                                                                   #                
#                                                                                                   #
#                        Protection coverage of Elasmobranch biodiversity                           #
#                                based on the global network of MPAs                                #
#                                                                                                   #
#                       MS : Pimiento et al. Functional diversity of sharks and rays                #
#                       is highly vulnerable and supported by uniquespecies and locations worldwide #
#                       Submitted to Nature Communications August 2023                              #
#                                                                                                   #
#                 Authors : Théophile Mouton, Fabien Leprieur, Camille Albouy                       #
#                           Catalina Pimiento                                                       #
#####################################################################################################   

# set working directory
setwd(".../MPA/") # your directory 

# Loading the required libraries 

lib_vect <- c("raster","sf","rgeos","maptools","gstat", "reshape", "tidyverse", "ggplot2", "udunits2");
sapply(lib_vect,library,character.only=TRUE)

##################################################################################
#   Part 1 : Filter the global MPA database to include only categories I to III  #
##################################################################################

# 1-1 Global MPAs IUCN status I to III

MPAs_1=sf::st_read(dsn="data/MPA_layers/WDPA_WDOECM_Jan2023_Public_marine_shp_0", 
                   layer = "WDPA_WDOECM_Jan2023_Public_marine_shp-polygons") #This opens only a single layer 

MPAs_1=MPAs_1[(MPAs_1$IUCN_CAT %in% c("Ia", "Ib", "II", "III")),] #Select MPAs with IUCN status I to III

MPAs_2=sf::st_read(dsn="data/MPA_layers/WDPA_WDOECM_Jan2023_Public_marine_shp_1", 
                   layer = "WDPA_WDOECM_Jan2023_Public_marine_shp-polygons") #This opens only a single layer 

MPAs_2=MPAs_2[(MPAs_2$IUCN_CAT %in% c("Ia", "Ib", "II", "III")),] #Select MPAs with IUCN status I to III

MPAs_3=sf::st_read("data/MPA_layers/WDPA_WDOECM_Jan2023_Public_marine_shp_2", 
                   layer = "WDPA_WDOECM_Jan2023_Public_marine_shp-polygons") #This opens only a single layer 

MPAs_3=MPAs_3[(MPAs_3$IUCN_CAT %in% c("Ia", "Ib", "II", "III")),] #Select MPAs with IUCN status I to III

# 1- 2 : Rbind the three layers
MPAs_1st_set=rbind(MPAs_1, MPAs_2)
MPAs_Final_set_I_III=rbind(MPAs_1st_set, MPAs_3)

MPAs_Final_set_I_III=MPAs_Final_set_I_III[!(MPAs_Final_set_I_III$DESIG_ENG
                                            %in% c("Wetland Protected Area", 
                                                   "Forest Nature Reserve",
                                                   "Waterfowl gathering area",
                                                   "Seabird Sanctuary",
                                                   "Managed Flora Reserve",
                                                   "Migratory Bird Sanctuary",
                                                   "Shell Reserve",
                                                   "Archaeological Preserve State Park",
                                                   "Botanical State Park",
                                                   "Bird Reserve", 
                                                   "Special Botanical Reserve",
                                                   "Special Use Forest",
                                                   "Old-growth forest",
                                                   "Wetland Site",
                                                   "Cave",
                                                   "Forest Ecosystem Reserve",
                                                   "Bird Sanctuary",
                                                   "Ramsar Site, Wetland of International Importance",
                                                   "Forest Managed Biological Reserve",
                                                   "Iguana Sanctuary",
                                                   "Bird Sanctuary under the Wild Birds Protection Act 1932",
                                                   "Wetland Park")),] #Remove MPAs not of interest for shark protection

saveRDS(MPAs_Final_set_I_III, file="outputs/MPAs_Final_set_I_III.RDS")

#############################################################
# Part 2 : Protection coverage of elasmobranch biodiversity # 
#############################################################

#2.1  Loading a data frame of all shark biodiversity indices at the global scale

shark_data <- readRDS("data/Data_end_shark_june2023.rds", refhook = NULL)

#2.2 Raster computation for each biodiversity metric (need for for the extraction step 

ProjWGS84 <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"

shark_SR <- rasterFromXYZ(data.frame(x=shark_data$x,y=shark_data$y,z=shark_data$SR),crs=ProjWGS84)
shark_PD <- rasterFromXYZ(data.frame(x=shark_data$x,y=shark_data$y,z=shark_data$PD),crs=ProjWGS84)
shark_FRic <- rasterFromXYZ(data.frame(x=shark_data$x,y=shark_data$y,z=shark_data$FRic),crs=ProjWGS84)
shark_FUn <- rasterFromXYZ(data.frame(x=shark_data$x,y=shark_data$y,z=shark_data$FUn),crs=ProjWGS84)
shark_FSpe <- rasterFromXYZ(data.frame(x=shark_data$x, y=shark_data$y,z=shark_data$FSp),crs=ProjWGS84)
shark_PUn <- rasterFromXYZ(data.frame(x=shark_data$x,y=shark_data$y,z=shark_data$PUn),crs=ProjWGS84)
shark_SR_treat<-rasterFromXYZ(data.frame(x=shark_data$x,y=shark_data$y,z=shark_data$SR_treat),crs=ProjWGS84)
shark_SR_Top_EDGE2<-rasterFromXYZ(data.frame(x=shark_data$x,y=shark_data$y,z=shark_data$SR_Top_EDGE2),crs=ProjWGS84)
shark_SR_Top_FUSE100<-rasterFromXYZ(data.frame(x=shark_data$x,y=shark_data$y,z=shark_data$SR_Top_FUSE100),crs=ProjWGS84) 
shark_SR_Top_ED2<-rasterFromXYZ(data.frame(x=shark_data$x,y=shark_data$y,z=shark_data$SR_Top_ED2_25),crs=ProjWGS84)
shark_SR_Top_FUN25<-rasterFromXYZ(data.frame(x=shark_data$x,y=shark_data$y,z=shark_data$SR_Top_FUN25),crs=ProjWGS84)
shark_SR_Top_FSp25<-rasterFromXYZ(data.frame(x=shark_data$x,y=shark_data$y,z=shark_data$SR_Top_FSp25),crs=ProjWGS84)

#2.3 Loading global MPA network with categories I to III

MPAs_Final_set_I_III  <- readRDS("outputs/MPAs_Final_set_I_III.RDS")
MPA_list <- lapply(MPAs_Final_set_I_III$geometry, function(i) {as(i,"Spatial")})
MPA_area <- sapply(MPA_list, function(x) {gArea(x)})
MPA_final_list <- MPA_list[order(MPA_area,decreasing=T)]
m1 <- do.call(bind,MPA_final_list)

#2.4 Extraction of the biodiversity values for each MPA's polygon

SR <- extract(shark_SR,m1)
SR_O<- na.omit(do.call(c,SR))
SR_O <- SR_O[order(SR_O)]
saveRDS(SR_O, file="outputs/SR_O_I_III")

PD <- extract(shark_PD,m1)
PD_O<- na.omit(do.call(c,PD))
PD_O <- PD_O[order(PD_O)]
saveRDS(PD_O, file="outputs/PD_O_I_III")

Fric<- extract(shark_FRic,m1)
Fric_O<- na.omit(do.call(c,Fric))
Fric_O <- Fric_O[order(Fric_O)]
saveRDS(Fric_O, file="outputs/Fric_O_I_III")

FSpe<- extract(shark_FSpe,m1)
FSpe_O<- na.omit(do.call(c,FSpe))
FSpe_O <- FSpe_O[order(FSpe_O)]
saveRDS(FSpe_O, file="outputs/FSpe_O_I_III")

PUn <- extract(shark_PUn,m1)
PUn_O<- na.omit(do.call(c,PUn))
PUn_O <- PUn_O[order(PUn_O)]
saveRDS(PUn_O, file="outputs/PUn_O_I_III")

FUn<- extract(shark_FUn,m1)
FUn_O<- na.omit(do.call(c,FUn))
FUn_O <- FUn_O[order(FUn_O)]
saveRDS(FUn_O, file="outputs/FUn_O_I_III")

SR_treat <- extract(shark_SR_treat,m1)
SR_treat_O<- na.omit(do.call(c,SR_treat))
SR_treat_O <- SR_treat_O[order(SR_treat_O)]
saveRDS(SR_treat_O, file="outputs/SR_treat_O_I_III")

SR_Top_EDGE2 <- extract(shark_SR_Top_EDGE2,m1)
SR_Top_EDGE2_O<- na.omit(do.call(c,SR_Top_EDGE2))
SR_Top_EDGE2_O <- SR_Top_EDGE2_O[order(SR_Top_EDGE2_O)]
saveRDS(SR_Top_EDGE2_O, file="outputs/SR_Top_EDGE2_O_I_III")

SR_Top_FUSE100 <- extract(shark_SR_Top_FUSE100,m1)
SR_Top_FUSE100_O<- na.omit(do.call(c,SR_Top_FUSE100))
SR_Top_FUSE100_O <- SR_Top_FUSE100_O[order(SR_Top_FUSE100_O)]
saveRDS(SR_Top_FUSE100_O, file="outputs/SR_Top_FUSE100_O_I_III")

SR_Top_ED2<- extract(shark_SR_Top_ED2,m1)
SR_Top_ED2_O<- na.omit(do.call(c,SR_Top_ED2))
SR_Top_ED2_O <- SR_Top_ED2_O[order(SR_Top_ED2_O)]
saveRDS(SR_Top_ED2_O, file="outputs/SR_Top_ED2_O_I_III")

SR_Top_FUN25<- extract(shark_SR_Top_FUN25,m1)
SR_Top_FUN25_O<- na.omit(do.call(c,SR_Top_FUN25))
SR_Top_FUN25_O <- SR_Top_FUN25_O[order(SR_Top_FUN25_O)]
saveRDS(SR_Top_FUN25_O, file="outputs/SR_Top_FUN25_O_I_III")

SR_Top_FSp25<- extract(shark_SR_Top_FSp25,m1)
SR_Top_FSp25_O<- na.omit(do.call(c,SR_Top_FSp25))
SR_Top_FSp25_O<- SR_Top_FSp25_O[order(SR_Top_FSp25_O)]
saveRDS(SR_Top_FSp25_O, file="outputs/SR_Top_FSp25_O_I_III")

#2.5 Data frame compiling biodiversity values inside MPAs

SR_O<-readRDS("outputs/SR_O_I_III")
PD_O<-readRDS("outputs/PD_O_I_III")
Fric_O<-readRDS("outputs/Fric_O_I_III")
FUn_O<-readRDS("outputs/FUn_O_I_III")
FSpe_O<-readRDS("outputs/FSpe_O_I_III")
PUn_O<-readRDS("outputs/PUn_O_I_III")
SR_treat_O<-readRDS("outputs/SR_treat_O_I_III")
SR_Top_EDGE2_O<-readRDS("outputs/SR_Top_EDGE2_O_I_III")
SR_Top_FUSE100_O<-readRDS("outputs/SR_Top_FUSE100_O_I_III")
SR_Top_ED2_O<-readRDS("outputs/SR_Top_ED2_O_I_III")
SR_Top_FUN25_O<-readRDS("outputs/SR_Top_FUN25_O_I_III")
SR_Top_FSp25_O<-readRDS("outputs/SR_Top_FSp25_O_I_III")

cells<-rep(1:length(SR_O))

data_shark_mpa_I_III<-data.frame(cells=cells, SR_O=SR_O, PD_O=PD_O, Fric_O=Fric_O, FSpe_O=FSpe_O, PUn_O=PUn_O, SR_treat_O=SR_treat_O, SR_Top_EDGE2_O, SR_Top_FUSE100_O, SR_Top_ED2_O, 
                   SR_Top_FUN25_O, SR_Top_FSp25_O, FUn_O=FUn_O)


saveRDS(data_shark_mpa_I_III, file="outputs/data_shark_mpa_I_III_june2023")

#2.6 Quantiles for each biodiversity metric

q_PD<-quantile(shark_data$PD, probs = seq(0, 1, 0.025), na.rm = T,
         names = TRUE)
q_SR<-quantile(shark_data$SR, probs = seq(0, 1, 0.025), na.rm = T,
         names = TRUE)
q_FRic<-quantile(shark_data$FRic, probs = seq(0, 1, 0.025), na.rm = T,
         names = TRUE)
q_FSp<-quantile(shark_data$FSp, probs = seq(0, 1, 0.025), na.rm = T,
         names = TRUE)
q_FUn<-quantile(shark_data$FUn, probs = seq(0, 1, 0.025), na.rm = T,
         names = TRUE)
q_PUn<-quantile(shark_data$PUn, probs = seq(0, 1, 0.025), na.rm = T,
         names = TRUE)
q_SR_treat<-quantile(shark_data$SR_treat, probs = seq(0, 1, 0.025), na.rm = T,
         names = TRUE)
q_SR_Top_EDGE2<-quantile(shark_data$SR_Top_EDGE2, probs = seq(0, 1, 0.025), na.rm = T,
         names = TRUE)
q_SR_Top_FUSE100<-quantile(shark_data$SR_Top_FUSE100, probs = seq(0, 1, 0.025), na.rm = T,
         names = TRUE)
q_SR_Top_ED2<-quantile(shark_data$SR_Top_ED2, probs = seq(0, 1, 0.025), na.rm = T,
         names = TRUE)
q_SR_Top_FUN25<-quantile(shark_data$SR_Top_FUN25, probs = seq(0, 1, 0.025), na.rm = T,
         names = TRUE)
q_SR_Top_FSp25<-quantile(shark_data$SR_Top_FSp25, probs = seq(0, 1, 0.025), na.rm = T,
         names = TRUE)

# 2.7 Figure S16 in Supp Mat   

gSR<-ggplot(data_shark_mpa_I_III, aes(x=cells/max(cells)*100, y=SR_O)) + geom_point(size=0.8) +
  geom_hline(yintercept=q_SR["97.5%"], color="red") + geom_hline(yintercept=q_SR["90%"], linetype="dashed", color="red") +
  xlab("Proportion of cells in protected areas (%)") + ylab("Species richness") 
 
gPD<-ggplot(data_shark_mpa_I_III, aes(x=cells/max(cells)*100, y=PD_O)) + geom_point(size=0.8) +
  geom_hline(yintercept=q_PD["97.5%"], color="red") + geom_hline(yintercept=q_PD["90%"], linetype="dashed", color="red") +
  xlab("Proportion of cells in protected areas (%)") + ylab("Phylogenetic diversity")

gFRic<-ggplot(data_shark_mpa_I_III, aes(x=cells/max(cells)*100, y=Fric_O)) + geom_point(size=0.8) +
  geom_hline(yintercept=q_FRic["97.5%"], color="red") + geom_hline(yintercept=q_FRic["90%"], linetype="dashed", color="red") +
  xlab("Proportion of cells in protected areas (%)") + ylab("Functional richness")

gFUn<-ggplot(data_shark_mpa_I_III, aes(x=cells_FUn_O/max(cells_FUn_O)*100, y=FUn_O)) + geom_point(size=0.8) +
  geom_hline(yintercept=q_FUn["97.5%"], color="red") + geom_hline(yintercept=q_FUn["90%"], linetype="dashed", color="red") +
  xlab("Proportion of cells in protected areas (%)") + ylab("Functional uniqueness")

gFSpe<-ggplot(data_shark_mpa_I_III, aes(x=cells/max(cells)*100, y=FSpe_O)) + geom_point(size=0.8) +
  geom_hline(yintercept=q_FSp["97.5%"], color="red") + geom_hline(yintercept=q_FSp["90%"], linetype="dashed", color="red") +
  xlab("Proportion of cells in protected areas (%)") + ylab("Functional specialization")

gPUn<-ggplot(data_shark_mpa_I_III, aes(x=cells/max(cells)*100, y=PUn_O)) + geom_point(size=0.8) +
  geom_hline(yintercept=q_PUn["97.5%"], color="red") + geom_hline(yintercept=q_PUn["90%"], linetype="dashed", color="red") +
  xlab("Proportion of cells in protected areas (%)") + ylab("Phylogenetic Uniqueness")

gSR_treat_O<-ggplot(data_shark_mpa_I_III, aes(x=cells/max(cells)*100, y=SR_treat_O)) + geom_point(size=0.8) +
  geom_hline(yintercept=q_SR_treat["97.5%"], color="red") + geom_hline(yintercept=q_SR_treat["90%"], linetype="dashed", color="red") +
  xlab("Proportion of cells in protected areas (%)") + ylab("Number of threatened species")

gSR_Top_EDGE2_O<-ggplot(data_shark_mpa_I_III, aes(x=cells/max(cells)*100, y=SR_Top_EDGE2_O)) + geom_point(size=0.8) +
  geom_hline(yintercept=q_SR_Top_EDGE2["97.5%"], color="red") + geom_hline(yintercept=q_SR_Top_EDGE2["90%"], linetype="dashed", color="red") +
  xlab("Proportion of cells in protected areas (%)") + ylab("Number of top 25% EDGE2 species")

gSR_Top_FUSE100_O<-ggplot(data_shark_mpa_I_III, aes(x=cells/max(cells)*100, y=SR_Top_FUSE100_O)) + geom_point(size=0.8) +
  geom_hline(yintercept=q_SR_Top_FUSE100["97.5%"], color="red") + geom_hline(yintercept=q_SR_Top_FUSE100["90%"], linetype="dashed", color="red") +
  xlab("Proportion of cells in protected areas (%)") + ylab("Number of top 25% FUSE species")

gSR_Top_ED2_O<-ggplot(data_shark_mpa_I_III, aes(x=cells/max(cells)*100, y=SR_Top_ED2_O)) + geom_point(size=0.8) +
  geom_hline(yintercept=q_SR_Top_ED2["97.5%"], color="red") + geom_hline(yintercept=q_SR_Top_ED2["90%"], linetype="dashed", color="red") +
  xlab("Proportion of cells in protected areas (%)") + ylab("Number of top 25% ED2 species")

gSR_Top_FUN25_O<-ggplot(data_shark_mpa_I_III, aes(x=cells/max(cells)*100, y=SR_Top_FUN25_O)) + geom_point(size=0.8) +
  geom_hline(yintercept=q_SR_Top_FUN25["97.5%"], color="red") + geom_hline(yintercept=q_SR_Top_FUN25["90%"], linetype="dashed", color="red") +
  xlab("Proportion of cells in protected areas (%)") + ylab("Number of top 25% FUn species")

gSR_Top_Fsp25_O<-ggplot(data_shark_mpa_I_III, aes(x=cells/max(cells)*100, y=SR_Top_FSp25_O)) + geom_point(size=0.8) +
  geom_hline(yintercept=q_SR_Top_FSp25["97.5%"], color="red") + geom_hline(yintercept=q_SR_Top_FSp25["90%"], linetype="dashed", color="red") +
  xlab("Proportion of cells in protected areas (%)") + ylab("Number of top 25% FSp species")


require(grid); require(ggpubr)

figure <- ggpubr::ggarrange(gSR + rremove("xlab"), 
                            gPD + rremove("xlab"), 
                            gFRic + rremove("xlab"), 
                            gPUn + rremove("xlab"), 
                            gFUn + rremove("xlab"), 
                            gFSpe + rremove("xlab"), # remove axis labels from plots
                            gSR_treat_O + rremove("xlab"),
                            gSR_Top_EDGE2_O + rremove("xlab"),
                            gSR_Top_FUSE100_O + rremove("xlab"),
                            gSR_Top_ED2_O + rremove("xlab"),
                            gSR_Top_FUN25_O + rremove("xlab"),
                            gSR_Top_Fsp25_O + rremove("xlab"),
                            labels = c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L"),
                            ncol = 3, nrow = 4,
                            common.legend = TRUE, legend = "right",
                            align = "hv", 
                            font.label = list(size = 12, color = "black", face = "bold", family = NULL, position = "top"))

figure = annotate_figure(figure, bottom = textGrob("Proportion of cells in protected areas (%)", gp = gpar(cex = 1)))

ggsave(file = "outputs/Fig_ S16_june2023.pdf", 
       figure, 
       height = 360, width = 360, units="mm")

#2.8 Identification of the % values in Fig S17 = % of hotspots vs. non hotspots inside the MPAs

data_shark_mpa_I_IV<-readRDS("outputs/data_shark_mpa_I_III_june2023")

value_inter_SR<-round((length(which((data_shark_mpa_I_IV$SR_O>=q_SR["97.5%"])))/length(data_shark_mpa_I_IV$SR_O))*100, 2)
value_inter_PD<-round((length(which((data_shark_mpa_I_IV$PD_O>=q_PD["97.5%"])))/length(data_shark_mpa_I_IV$PD_O))*100, 2) 
value_inter_FRiC<-round((length(which((data_shark_mpa_I_IV$Fric_O>=q_FRic["97.5%"])))/length(data_shark_mpa_I_IV$Fric_O))*100, 2) 
value_inter_FSp<-round((length(which((data_shark_mpa_I_IV$FSpe_O>=q_FSp["97.5%"])))/length(data_shark_mpa_I_IV$FSpe_O))*100, 2)  
value_inter_FUn<-round((length(which((data_shark_mpa_I_IV$FUn_O>=q_FUn["97.5%"])))/length(data_shark_mpa_I_IV_FUn$FUn_O))*100,2) 
value_inter_PUn<-round((length(which((data_shark_mpa_I_IV$PUn_O>=q_PUn["97.5%"])))/length(data_shark_mpa_I_IV$PUn_O))*100,2)
value_inter_SR_treat<-round((length(which((data_shark_mpa_I_IV$SR_treat_O)>=q_SR_treat["97.5%"]))/length(data_shark_mpa_I_IV$SR_treat_O))*100,2)
value_inter_SR_Top_EDGE2<-round((length(which((data_shark_mpa_I_IV$SR_Top_EDGE2_O>=q_SR_Top_EDGE2["97.5%"])))/length(data_shark_mpa_I_IV$SR_Top_EDGE2_O))*100,2)
value_inter_SR_Top_FUSE<-round((length(which((data_shark_mpa_I_IV$SR_Top_FUSE100_O>=q_SR_Top_FUSE100["97.5%"])))/length(data_shark_mpa_I_IV$SR_Top_FUSE100_O))*100,2)
value_inter_SR_Top_FUn<-round((length(which((data_shark_mpa_I_IV$SR_Top_FUN25_O>=q_SR_Top_FUN25["97.5%"])))/length(data_shark_mpa_I_IV$SR_Top_FUN25_O))*100,2)
value_inter_SR_Top_FSp<-round((length(which((data_shark_mpa_I_IV$SR_Top_FSp25_O>=q_SR_Top_FSp25["97.5%"])))/length(data_shark_mpa_I_IV$SR_Top_FSp25_O))*100,2)
value_inter_SR_Top_ED2<-round((length(which((data_shark_mpa_I_IV$SR_Top_ED2_O>=q_SR_Top_ED2["97.5%"])))/length(data_shark_mpa_I_IV$SR_Top_ED2_O))*100,2)

Value<-c(value_inter_SR, value_inter_PD, value_inter_FRiC, value_inter_FSp, value_inter_FUn,value_inter_PUn,
         value_inter_SR_treat, value_inter_SR_Top_EDGE2, value_inter_SR_Top_FUSE,value_inter_SR_Top_ED2, value_inter_SR_Top_FUn,
         value_inter_SR_Top_FSp)

Index <-c("SR", "PD", "FRiC", "FSp", "FUn", "PUn", "SR_threat", "SR_Top_EDGE2", "SR_Top_FUSE",
          "SR_Top_ED2", "SR_Top_FUn", "SR_Top_FSp")

value_curve_MPA_I_III_june2023<-as.data.frame(cbind(Index, Value))

saveRDS(value_curve_MPA_I_III_june2023, file="value_curve_MPA_I_III_june2023")


# 2.9 Percentage of hotspots protected (Fig. S18) according to different threshold values to define biodiversity hotspots (2.5%, 5%, 10% and 20%)

# 2.9.1  Loading the required data
shark_data <- readRDS("data/Data_end_shark_june2023.RDS", refhook = NULL) # biodiversity values for all the cells 
data_k<-readRDS("outputs/data_shark_mpa_I_III_june2023") # biodiversity values in MPAs 

# 2.9.2 Quantiles for each biodiversity metric

q_PD<-quantile(shark_data$PD, probs = seq(0, 1, 0.025), na.rm = T,
               names = TRUE)
q_SR<-quantile(shark_data$SR, probs = seq(0, 1, 0.025), na.rm = T,
               names = TRUE)
q_FRic<-quantile(shark_data$FRic, probs = seq(0, 1, 0.025), na.rm = T,
                 names = TRUE)
q_FSp<-quantile(shark_data$FSp, probs = seq(0, 1, 0.025), na.rm = T,
                names = TRUE)
q_FUn<-quantile(shark_data$FUn, probs = seq(0, 1, 0.025), na.rm = T,
                names = TRUE)
q_PUn<-quantile(shark_data$PUn, probs = seq(0, 1, 0.025), na.rm = T,
                names = TRUE)
q_SR_treat<-quantile(shark_data$SR_treat, probs = seq(0, 1, 0.025), na.rm = T,
                     names = TRUE)
q_SR_Top_EDGE2<-quantile(shark_data$SR_Top_EDGE2, probs = seq(0, 1, 0.025), na.rm = T,
                         names = TRUE)
q_SR_Top_FUSE100<-quantile(shark_data$SR_Top_FUSE100, probs = seq(0, 1, 0.025), na.rm = T,
                           names = TRUE)
q_SR_Top_ED2<-quantile(shark_data$SR_Top_ED2, probs = seq(0, 1, 0.025), na.rm = T,
                       names = TRUE)
q_SR_Top_FUN25<-quantile(shark_data$SR_Top_FUN25, probs = seq(0, 1, 0.025), na.rm = T,
                         names = TRUE)
q_SR_Top_FSp25<-quantile(shark_data$SR_Top_FSp25, probs = seq(0, 1, 0.025), na.rm = T,
                         names = TRUE)

# 2.9.3 Calculate percentage of cells within hotspot / Fig S18

df = rbind.data.frame(
  c(
    length(data_k[data_k$SR_O > q_SR["97.5%"],2])/length(shark_data[shark_data$SR > q_SR["97.5%"],3])*100, # quantile 2.5% for species richness
    length(data_k[data_k$SR_O > q_SR["95%"],2])/length(shark_data[shark_data$SR > q_SR["95%"],3])*100, # quantile 5% for species richness
    length(data_k[data_k$SR_O > q_SR["90%"],2])/length(shark_data[shark_data$SR > q_SR["90%"],3])*100, #  quantile 10% for species richness
    length(data_k[data_k$SR_O > q_SR["80%"],2])/length(shark_data[shark_data$SR > q_SR["80%"],3])*100) # #  quantile 10% for species richness
  ,
  PD = c(
    length(data_k[data_k$PD_O > q_PD["97.5%"],3])/length(shark_data[shark_data$PD > q_PD["97.5%"],4])*100,
    length(data_k[data_k$PD_O > q_PD["95%"],3])/length(shark_data[shark_data$PD > q_PD["95%"],4])*100,
    length(data_k[data_k$PD_O > q_PD["90%"],3])/length(shark_data[shark_data$PD > q_PD["90%"],4])*100,
    length(data_k[data_k$PD_O > q_PD["80%"],3])/length(shark_data[shark_data$PD > q_PD["80%"],4])*100)
  ,
  FRic = c(
    length(data_k[data_k$Fric_O > q_FRic["97.5%"],4])/length(shark_data[shark_data$FRic > q_FRic["97.5%"],6])*100,
    length(data_k[data_k$Fric_O > q_FRic["95%"],4])/length(shark_data[shark_data$FRic > q_FRic["95%"],6])*100,
    length(data_k[data_k$Fric_O > q_FRic["90%"],4])/length(shark_data[shark_data$FRic > q_FRic["90%"],6])*100,
    length(data_k[data_k$Fric_O > q_FRic["80%"],4])/length(shark_data[shark_data$FRic > q_FRic["80%"],6])*100)
  ,
  PUn = c(
    length(data_k[data_k$PUn_O > q_PUn["97.5%"],6])/length(shark_data[shark_data$PUn > q_PUn["97.5%"],5])*100,
    length(data_k[data_k$PUn_O > q_PUn["95%"],6])/length(shark_data[shark_data$PUn > q_PUn["95%"],5])*100,
    length(data_k[data_k$PUn_O > q_PUn["90%"],6])/length(shark_data[shark_data$PUn > q_PUn["90%"],5])*100,
    length(data_k[data_k$PUn_O > q_PUn["80%"],6])/length(shark_data[shark_data$PUn > q_PUn["80%"],5])*100)
  ,
  FUn = c(
    length(data_k_FUn[data_k$FUn_O > q_FUn["97.5%"],2])/length(shark_data[shark_data$FUn > q_FUn["97.5%"],8])*100,
    length(data_k_FUn[data_k$FUn_O > q_FUn["95%"],2])/length(shark_data[shark_data$FUn > q_FUn["95%"],8])*100,
    length(data_k_FUn[data_k$FUn_O > q_FUn["90%"],2])/length(shark_data[shark_data$FUn > q_FUn["90%"],8])*100,
    length(data_k_FUn[data_k$FUn_O > q_FUn["80%"],2])/length(shark_data[shark_data$FUn > q_FUn["80%"],8])*100)
  ,
  FSpe = c(
    length(data_k[data_k$FSpe_O > q_FSp["97.5%"],5])/length(shark_data[shark_data$FSp > q_FSp["97.5%"],7])*100,
    length(data_k[data_k$FSpe_O > q_FSp["95%"],5])/length(shark_data[shark_data$FSp > q_FSp["95%"],7])*100,
    length(data_k[data_k$FSpe_O > q_FSp["90%"],5])/length(shark_data[shark_data$FSp > q_FSp["90%"],7])*100,
    length(data_k[data_k$FSpe_O > q_FSp["80%"],5])/length(shark_data[shark_data$FSp > q_FSp["80%"],7])*100)
  ,
  SR_treat = c(
    length(data_k[data_k$SR_treat_O > q_SR_treat["97.5%"],7])/length(shark_data[shark_data$SR_treat > q_SR_treat["97.5%"],9])*100,
    length(data_k[data_k$SR_treat_O > q_SR_treat["95%"],7])/length(shark_data[shark_data$SR_treat > q_SR_treat["95%"],9])*100,
    length(data_k[data_k$SR_treat_O > q_SR_treat["90%"],7])/length(shark_data[shark_data$SR_treat > q_SR_treat["90%"],9])*100,
    length(data_k[data_k$SR_treat_O > q_SR_treat["80%"],7])/length(shark_data[shark_data$SR_treat > q_SR_treat["80%"],9])*100)
  ,
  Top_EDGE2 = c(
    length(data_k[data_k$SR_Top_EDGE2_O > q_SR_Top_EDGE2["97.5%"],8])/length(shark_data[shark_data$SR_Top_EDGE2 > q_SR_Top_EDGE2["97.5%"],10])*100,
    length(data_k[data_k$SR_Top_EDGE2_O > q_SR_Top_EDGE2["95%"],8])/length(shark_data[shark_data$SR_Top_EDGE2 > q_SR_Top_EDGE2["95%"],10])*100,
    length(data_k[data_k$SR_Top_EDGE2_O > q_SR_Top_EDGE2["90%"],8])/length(shark_data[shark_data$SR_Top_EDGE2 > q_SR_Top_EDGE2["90%"],10])*100,
    length(data_k[data_k$SR_Top_EDGE2_O > q_SR_Top_EDGE2["80%"],8])/length(shark_data[shark_data$SR_Top_EDGE2 > q_SR_Top_EDGE2["80%"],10])*100)
  ,
  Top_FUSE = c(
    length(data_k[data_k$SR_Top_FUSE100_O > q_SR_Top_FUSE100["97.5%"],9])/length(shark_data[shark_data$SR_Top_FUSE100 > q_SR_Top_FUSE100["97.5%"],11])*100,
    length(data_k[data_k$SR_Top_FUSE100_O > q_SR_Top_FUSE100["95%"],9])/length(shark_data[shark_data$SR_Top_FUSE100 > q_SR_Top_FUSE100["95%"],11])*100,
    length(data_k[data_k$SR_Top_FUSE100_O > q_SR_Top_FUSE100["90%"],9])/length(shark_data[shark_data$SR_Top_FUSE100 > q_SR_Top_FUSE100["90%"],11])*100,
    length(data_k[data_k$SR_Top_FUSE100_O > q_SR_Top_FUSE100["80%"],9])/length(shark_data[shark_data$SR_Top_FUSE100 > q_SR_Top_FUSE100["80%"],11])*100)
  ,
  Top_ED2 = c(
    length(data_k[data_k$SR_Top_ED2_O > q_SR_Top_ED2["97.5%"],10])/length(shark_data[shark_data$SR_Top_ED2 > q_SR_Top_ED2["97.5%"],12])*100,
    length(data_k[data_k$SR_Top_ED2_O > q_SR_Top_ED2["95%"],10])/length(shark_data[shark_data$SR_Top_ED2 > q_SR_Top_ED2["95%"],12])*100,
    length(data_k[data_k$SR_Top_ED2_O > q_SR_Top_ED2["90%"],10])/length(shark_data[shark_data$SR_Top_ED2 > q_SR_Top_ED2["90%"],12])*100,
    length(data_k[data_k$SR_Top_ED2_O > q_SR_Top_ED2["80%"],10])/length(shark_data[shark_data$SR_Top_ED2 > q_SR_Top_ED2["80%"],12])*100)
  ,
  Top_FUn = c(
    length(data_k[data_k$SR_Top_FUN25_O > q_SR_Top_FUN25["97.5%"],11])/length(shark_data[shark_data$SR_Top_FUN25 > q_SR_Top_FUN25["97.5%"],13])*100,
    length(data_k[data_k$SR_Top_FUN25_O > q_SR_Top_FUN25["95%"],11])/length(shark_data[shark_data$SR_Top_FUN25 > q_SR_Top_FUN25["95%"],13])*100,
    length(data_k[data_k$SR_Top_FUN25_O > q_SR_Top_FUN25["90%"],11])/length(shark_data[shark_data$SR_Top_FUN25 > q_SR_Top_FUN25["90%"],13])*100,
    length(data_k[data_k$SR_Top_FUN25_O > q_SR_Top_FUN25["80%"],11])/length(shark_data[shark_data$SR_Top_FUN25 > q_SR_Top_FUN25["80%"],13])*100)
  ,
  Top_FSp = c(
    length(data_k[data_k$SR_Top_FSp25_O > q_SR_Top_FSp25["97.5%"],12])/length(shark_data[shark_data$SR_Top_FSp25 > q_SR_Top_FSp25["97.5%"],14])*100,
    length(data_k[data_k$SR_Top_FSp25_O > q_SR_Top_FSp25["95%"],12])/length(shark_data[shark_data$SR_Top_FSp25 > q_SR_Top_FSp25["95%"],14])*100,
    length(data_k[data_k$SR_Top_FSp25_O > q_SR_Top_FSp25["90%"],12])/length(shark_data[shark_data$SR_Top_FSp25 > q_SR_Top_FSp25["90%"],14])*100,
    length(data_k[data_k$SR_Top_FSp25_O > q_SR_Top_FSp25["80%"],12])/length(shark_data[shark_data$SR_Top_FSp25 > q_SR_Top_FSp25["80%"],14])*100)
  
)

colnames(df) = c("2.5%", "5%", "10%", "20%")

df$BioInd = c("Species richness", "Phylogenetic diversity", "Functional richness", "Phylogenetic uniqueness", 
              "Functional uniqueness", "Functional specialisation",
              "Threatened species", "Top 25% EDGE2", "Top 25% FUSE", "Top 25% ED2", "Top 25% FUn", "Top 25% FSp")

df = reshape::melt.data.frame(df, id = "BioInd")

df$BioInd = factor(df$BioInd, levels=c("Species richness", "Phylogenetic diversity", "Functional richness", "Phylogenetic uniqueness", 
                                       "Functional uniqueness", "Functional specialisation",
                                       "Threatened species", "Top 25% EDGE2", "Top 25% FUSE", "Top 25% ED2", "Top 25% FUn", "Top 25% FSp"))

saveRDS(df, file="protected_biodiv_mpaI_III_june2023")

g_Hot = ggplot(df, aes(x=variable, y=value)) +
  geom_bar(position = position_dodge(), stat="identity") +
  scale_y_continuous(expand=c(0,0), limits = c(0,100)) + 
  scale_x_discrete(expand = c(0,0)) +
  facet_wrap(~BioInd, ncol = 3) +
  theme_bw() +
  theme(legend.position = "none",
        legend.title.align=0.5 ,
        plot.margin = unit(c(1,1,1,1), "lines"),
        strip.background =element_rect(fill="black"),
        strip.text = element_text(colour = 'white')) +
  xlab("Top cells\n     ") + ylab("% of hotspot protected")

ggsave(file = "Fig_S18_june2023.jpeg", 
       g_Hot, 
       height = 360*0.8, width = 360*0.7, units="mm")


# 2.10 Barplot of the Figure 4 in the main text 

library(reshape)
library(tidyverse)

MPA_one_tree <- readRDS("outputs/protected_biodiv_mpaI_III_june2023")

bars <- MPA_one_tree %>%
  filter(variable=="2.5%") %>%
  mutate(metric=BioInd) %>%
  mutate(metric = recode(BioInd,"Species richness"="SR","Phylogenetic diversity"="PD",
                         "Functional richness"="FRic","Phylogenetic uniqueness"="PUn",
                         "Functional uniqueness"="FUn","Functional specialisation"="FSp")) %>%
  mutate(type = if_else(str_detect(metric, "Top"),
                        "species", "community"))%>%
  mutate(type = replace(type, metric=="Threatened species", "species")) %>%
  mutate(inside_MPAs = value ) %>%
  
  mutate(outside_MPAs=100-inside_MPAs) %>%
  pivot_longer(cols = inside_MPAs:outside_MPAs, names_to = "cells", values_to = "Hotspot cells") %>%
  mutate(cells = factor(cells))

bars$metric <- fct_reorder(bars$metric, bars$value)

cols <- c("#177e89","#ffc857", "#084c61","#db3a34")

pdf("Fig_4_barplots.pdf", 
    width = 7, height = 8) 

ggplot(bars, 
       aes(y=`Hotspot cells`, x=metric))+
  geom_bar(position="stack", stat="identity", aes(fill= cells))+
  scale_color_manual(values= cols[4:3])+
  scale_fill_manual(values= cols[4:3])+
  coord_flip()+
  theme_bw()+
  theme_classic()+
  labs(x = NULL, y= "Hotspot cells (%)")

dev.off()


#2.11 Fig_S17

intersects <- readRDS("outputs/value_curve_MPA_I_III_june2023")  %>%
  mutate(metric = recode(as_factor(Index),"SR_threat"="Threatened species","SR_Top_EDGE2"="Top 25% EDGE2",
                         "SR_Top_FUSE"="Top 25% FUSE","SR_Top_ED2"="Top 25% ED2",
                         "SR_Top_FUn"="Top 25% FUn","SR_Top_FSp"="Top 25% FSp")) %>%
  mutate(type = if_else(str_detect(metric, "Top"), 
                        "species", "community"))%>%
  mutate(type = replace(type, metric=="Threatened species", "species")) %>%
  mutate(hotspots=as.numeric(Value))%>%
  mutate(not_hotspots=100-hotspots) %>%
  pivot_longer(cols = hotspots:not_hotspots, names_to = "cells", values_to = "MPA cells")

intersects$metric <- fct_reorder(intersects$metric, as.numeric(intersects$Value))

intersects %>% 
  filter(cells=="not_hotspots")%>%
  summarise(mean=mean(`MPA cells`))

pdf("Fig_S17_june2023.pdf", 
    width = 6, height = 4)

ggplot(intersects, 
       aes(y=`MPA cells`, x=metric))+
  geom_bar(position="stack", stat="identity", aes(fill= `cells`))+
  scale_color_manual(values= cols[4:3])+
  scale_fill_manual(values= cols[4:3])+
  coord_flip()+
  theme_bw()+
  theme_classic()+
  labs(x = NULL, y= "MPA cells (%)")

dev.off()




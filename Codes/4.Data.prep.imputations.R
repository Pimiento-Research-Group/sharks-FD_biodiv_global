# Data prep for imputations
# Date: 26 Mar 2020
# Author: Catalina Pimiento


library(readxl)
library(stringr)
library(tidyverse)
library(dplyr)
library(janitor)

# Load traits dataset
  # traits.raw<-read_xlsx("~/Dropbox (Smithsonian)/SHARKS/Traits/shark_traits.xlsx")
  # traits.new <- read_xlsx("~/Dropbox (Smithsonian)/SHARKS/Traits/shark_traits_corrected.xlsx") 
    traits.curated<-read_xlsx("~/Dropbox (Smithsonian)/SHARKS/Sahrks FD/Analyses/R_imputs/shark_traits_curated.xlsx")

# Create data frame for each trait. All treated the same way as in FD analyes, plus phylogenetic info.

  # orders<-traits.raw %>%
    orders<-traits.curated %>%
    dplyr::select(Species, Order) %>%
    mutate(Order.Carcharhiniformes=0, Order.Heterodontiformes=0, Order.Hexanchiformes=0, Order.Lamniformes=0, Order.Myliobatiformes=0,
           Order.Orectolobiformes=0, Order.Pristiophoriformes=0, Order.Rajiformes=0, Order.Rhinopristiformes=0,
           Order.Squaliformes=0, Order.Squatiniformes=0, Order.Torpediniformes=0) %>%
    mutate(Order.Carcharhiniformes=replace(Order.Carcharhiniformes, Order=="Carcharhiniformes",1)) %>%
    mutate(Order.Heterodontiformes=replace(Order.Heterodontiformes, Order=="Heterodontiformes",1)) %>%
    mutate(Order.Hexanchiformes=replace(Order.Hexanchiformes, Order=="Hexanchiformes",1)) %>%
    mutate(Order.Lamniformes=replace(Order.Lamniformes, Order=="Lamniformes",1)) %>%
    mutate(Order.Myliobatiformes=replace(Order.Myliobatiformes, Order=="Myliobatiformes",1)) %>%
    mutate(Order.Orectolobiformes=replace(Order.Orectolobiformes, Order=="Orectolobiformes",1)) %>%
    mutate(Order.Pristiophoriformes=replace(Order.Pristiophoriformes, Order=="Pristiophoriformes",1)) %>%
    mutate(Order.Rajiformes=replace(Order.Rajiformes, Order=="Rajiformes",1)) %>%
    mutate(Order.Rhinopristiformes=replace(Order.Rhinopristiformes, Order=="Rhinopristiformes",1)) %>%
    mutate(Order.Squaliformes=replace(Order.Squaliformes, Order=="Squaliformes",1)) %>%
    mutate(Order.Squatiniformes=replace(Order.Squatiniformes, Order=="Squatiniformes",1)) %>%
    mutate(Order.Torpediniformes=replace(Order.Torpediniformes, Order=="Torpediniformes",1)) %>%
    dplyr::select(-Order)
    # orders$Order <- as.numeric(as.factor(orders$Order))

  
  # families<-traits.raw %>%
    families<-traits.curated%>%
    dplyr::select(Species, Family) %>%
    mutate(Family.Alopiidae=0, Family.Anacanthobatidae=0, Family.Brachaeluridae=0, Family.Carcharhinidae=0,
           Family.Centrophoridae=0, Family.Cetorhinidae=0, Family.Chlamydoselachidae=0, Family.Dalatiidae=0,
           Family.Dasyatidae=0, Family.Echinorhinidae=0, Family.Ginglymostomatidae=0, Family.Gurgesiellidae=0,
           Family.Gymnuridae=0, Family.Hemigaleidae=0, Family.Hemiscylliidae=0, Family.Heterodontidae=0,
           Family.Hexanchidae=0, Family.Hexatrygonidae=0, Family.Hexatrygonidae=0, Family.Hypnidae=0,
           Family.Lamnidae=0, Family.Leptochariidae=0, Family.Megachasmidae=0, Family.Mitsukurinidae=0,
           Family.Mobulidae=0, Family.Myliobatidae=0, Family.Narcinidae=0, Family.Narkidae=0,
           Family.Odontaspididae=0, Family.Orectolobidae=0, Family.Oxynotidae=0, Family.Parascylliidae=0,
           Family.Platyrhinidae=0, Family.Plesiobatidae=0, Family.Potamotrygonidae=0, Family.Pristidae=0,
           Family.Pristiophoridae=0, Family.Proscylliidae=0, Family.Pseudocarchariidae=0, Family.Rajidae=0,
           Family.Rhincodontidae=0, Family.Rhinidae=0, Family.Rhinobatidae=0, Family.Rhinopteridae=0,
           Family.Rhynchobatidae=0, Family.Scyliorhinidae=0, Family.Somniosidae=0, Family.Sphyrnidae=0,
           Family.Squalidae=0, Family.Squatinidae=0, Family.Stegostomatidae=0, Family.Torpedinidae=0,
           Family.Triakidae=0, Family.Urolophidae=0, Family.Urolophidae=0, Family.Zanobatidae=0,
           Family.Urotrygonidae=0, Family.Arhynchobatidae=0, Family.Pseudotriakidae=0, Family.Etmopteridae=0) %>%
    mutate(Family.Heterodontidae=replace(Family.Heterodontidae, Family=="Heterodontidae",1)) %>%
    mutate(Family.Chlamydoselachidae=replace(Family.Chlamydoselachidae, Family=="Chlamydoselachidae",1)) %>%
    mutate(Family.Hexanchidae=replace(Family.Hexanchidae, Family=="Hexanchidae",1)) %>%
    mutate(Family.Alopiidae=replace(Family.Alopiidae, Family=="Alopiidae",1)) %>%
    mutate(Family.Odontaspididae=replace(Family.Odontaspididae, Family=="Odontaspididae",1)) %>%
    mutate(Family.Lamnidae=replace(Family.Lamnidae, Family=="Lamnidae",1)) %>%
    mutate(Family.Cetorhinidae=replace(Family.Cetorhinidae, Family=="Cetorhinidae",1)) %>%
    mutate(Family.Megachasmidae=replace(Family.Megachasmidae, Family=="Megachasmidae",1)) %>%
    mutate(Family.Mitsukurinidae=replace(Family.Mitsukurinidae, Family=="Mitsukurinidae",1)) %>%
    mutate(Family.Pseudocarchariidae=replace(Family.Pseudocarchariidae, Family=="Pseudocarchariidae",1)) %>%
    mutate(Family.Myliobatidae=replace(Family.Myliobatidae, Family=="Myliobatidae",1)) %>%
    mutate(Family.Scyliorhinidae=replace(Family.Scyliorhinidae, Family=="Scyliorhinidae",1)) %>%
    mutate(Family.Dasyatidae=replace(Family.Dasyatidae, Family=="Dasyatidae",1)) %>%
    mutate(Family.Gymnuridae=replace(Family.Gymnuridae, Family=="Gymnuridae",1)) %>%
    mutate(Family.Potamotrygonidae=replace(Family.Potamotrygonidae, Family=="Potamotrygonidae",1)) %>%
    mutate(Family.Hexatrygonidae=replace(Family.Hexatrygonidae, Family=="Hexatrygonidae",1)) %>%
    mutate(Family.Carcharhinidae=replace(Family.Carcharhinidae, Family=="Carcharhinidae",1)) %>%
    mutate(Family.Mobulidae=replace(Family.Mobulidae, Family=="Mobulidae",1)) %>%
    mutate(Family.Hemigaleidae=replace(Family.Hemigaleidae, Family=="Hemigaleidae",1)) %>%
    mutate(Family.Plesiobatidae=replace(Family.Plesiobatidae, Family=="Plesiobatidae",1)) %>%
    mutate(Family.Proscylliidae=replace(Family.Proscylliidae, Family=="Proscylliidae",1)) %>%
    mutate(Family.Rhinopteridae=replace(Family.Rhinopteridae, Family=="Rhinopteridae",1)) %>%
    mutate(Family.Urolophidae=replace(Family.Urolophidae, Family=="Urolophidae",1)) %>%
    mutate(Family.Urotrygonidae=replace(Family.Urotrygonidae, Family=="Urotrygonidae",1)) %>%
    mutate(Family.Brachaeluridae=replace(Family.Brachaeluridae, Family=="Brachaeluridae",1)) %>%
    mutate(Family.Hemiscylliidae=replace(Family.Hemiscylliidae, Family=="Hemiscylliidae",1)) %>%
    mutate(Family.Parascylliidae=replace(Family.Parascylliidae, Family=="Parascylliidae",1)) %>%
    mutate(Family.Orectolobidae=replace(Family.Orectolobidae, Family=="Orectolobidae",1)) %>%
    mutate(Family.Ginglymostomatidae=replace(Family.Ginglymostomatidae, Family=="Ginglymostomatidae",1)) %>%
    mutate(Family.Rhincodontidae=replace(Family.Rhincodontidae, Family=="Rhincodontidae",1)) %>%
    mutate(Family.Stegostomatidae=replace(Family.Stegostomatidae, Family=="Stegostomatidae",1)) %>%
    mutate(Family.Pristiophoridae=replace(Family.Pristiophoridae, Family=="Pristiophoridae",1)) %>%
    mutate(Family.Rajidae=replace(Family.Rajidae, Family=="Rajidae",1)) %>%
    mutate(Family.Anacanthobatidae=replace(Family.Anacanthobatidae, Family=="Anacanthobatidae",1)) %>%
    mutate(Family.Arhynchobatidae=replace(Family.Arhynchobatidae, Family=="Arhynchobatidae",1)) %>%
    mutate(Family.Sphyrnidae=replace(Family.Sphyrnidae, Family=="Sphyrnidae",1)) %>%
    mutate(Family.Triakidae=replace(Family.Triakidae, Family=="Triakidae",1)) %>%
    mutate(Family.Pseudotriakidae=replace(Family.Pseudotriakidae, Family=="Pseudotriakidae",1)) %>%
    mutate(Family.Leptochariidae=replace(Family.Leptochariidae, Family=="Leptochariidae",1)) %>%
    mutate(Family.Gurgesiellidae=replace(Family.Gurgesiellidae, Family=="Gurgesiellidae",1)) %>%
    mutate(Family.Rhinobatidae=replace(Family.Rhinobatidae, Family=="Rhinobatidae",1)) %>%
    mutate(Family.Pristidae=replace(Family.Pristidae, Family=="Pristidae",1)) %>%
    mutate(Family.Rhinidae=replace(Family.Rhinidae, Family=="Rhinidae",1)) %>%
    mutate(Family.Rhynchobatidae=replace(Family.Rhynchobatidae, Family=="Rhynchobatidae",1)) %>%
    mutate(Family.Zanobatidae=replace(Family.Zanobatidae, Family=="Zanobatidae",1)) %>%
    mutate(Family.Etmopteridae=replace(Family.Etmopteridae, Family=="Etmopteridae",1)) %>%
    mutate(Family.Centrophoridae=replace(Family.Centrophoridae, Family=="Centrophoridae",1)) %>%
    mutate(Family.Somniosidae=replace(Family.Somniosidae, Family=="Somniosidae",1)) %>%
    mutate(Family.Squalidae=replace(Family.Squalidae, Family=="Squalidae",1)) %>%
    mutate(Family.Dalatiidae=replace(Family.Dalatiidae, Family=="Dalatiidae",1)) %>%
    mutate(Family.Echinorhinidae=replace(Family.Echinorhinidae, Family=="Echinorhinidae",1)) %>%
    mutate(Family.Oxynotidae=replace(Family.Oxynotidae, Family=="Oxynotidae",1)) %>%
    mutate(Family.Squatinidae=replace(Family.Squatinidae, Family=="Squatinidae",1)) %>%
    mutate(Family.Narcinidae=replace(Family.Narcinidae, Family=="Narcinidae",1)) %>%
    mutate(Family.Narkidae=replace(Family.Narkidae, Family=="Narkidae",1)) %>%
    mutate(Family.Hypnidae=replace(Family.Hypnidae, Family=="Hypnidae",1)) %>%
    mutate(Family.Platyrhinidae=replace(Family.Platyrhinidae, Family=="Platyrhinidae",1)) %>%
    mutate(Family.Torpedinidae=replace(Family.Torpedinidae, Family=="Torpedinidae",1))%>%
    dplyr::select(-Family)
    # families$Family <- as.numeric(as.factor(families$Family))
          
# binary traits   
    
  # habitat <- traits.new %>%
    habitat <- traits.curated %>%
    dplyr::select(Species, habitat) %>%
    mutate(Habitat.coastal=0, Habitat.oceanic=0) %>%
    
    mutate(Habitat.coastal=replace(Habitat.coastal, habitat=="all"|habitat=="shelf"
                                   |habitat=="shelf/slope", 1)) %>%
    mutate(Habitat.coastal=replace(Habitat.coastal, habitat=="NA", NA)) %>%
    
    mutate(Habitat.oceanic=replace(Habitat.oceanic, habitat=="all"|habitat=="shelf/slope"
                                   |habitat=="offshore"|habitat=="slope/offshore"|habitat=="slope",1)) %>%
    mutate(Habitat.oceanic=replace(Habitat.oceanic, habitat=="NA", NA)) %>%
    
    dplyr::select(-habitat)
  

  # vertical<- traits.new %>%
    vertical<- traits.curated %>%
    dplyr::select(Species, vertical) %>%
    mutate(Vertical.benthic=0, Vertical.pelagic=0) %>%
    
    mutate(Vertical.benthic=replace(Vertical.benthic, vertical=="benthic"|vertical=="benthopelagic",1)) %>%
    mutate(Vertical.benthic=replace(Vertical.benthic, vertical=="NA",NA)) %>%
    
    mutate(Vertical.pelagic=replace(Vertical.pelagic, vertical=="benthopelagic"|vertical=="pelagic", 1)) %>%
    mutate(Vertical.pelagic=replace(Vertical.pelagic, vertical=="NA", NA)) %>%
    
    dplyr::select(-vertical)
  
  # diet<- traits.new %>%
    diet<- traits.curated %>%
    dplyr::select(Species, diet) %>%
    mutate(Diet.plankton=0,Diet.inverts=0, Diet.fish=0, Diet.highverts=0) %>%
    
    mutate(Diet.plankton=replace(Diet.plankton, diet=="plankton"|diet=="plankton/fish", 1)) %>%
    mutate(Diet.plankton=replace(Diet.plankton, diet=="NA", NA)) %>%
    
    mutate(Diet.inverts=replace(Diet.inverts, diet=="inverts"|diet=="inverts/fish"
                                |diet=="inverts/fish/high verts", 1)) %>%
    mutate(Diet.inverts=replace(Diet.inverts, diet=="NA", NA)) %>%
    
    mutate(Diet.fish=replace(Diet.fish, diet=="fish"|diet=="inverts/fish"|diet=="plankton/fish"
                             |diet=="inverts/fish/high verts", 1)) %>%
    mutate(Diet.fish=replace(Diet.fish, diet=="NA", NA)) %>%
    
    mutate(Diet.highverts=replace(Diet.highverts, diet=="inverts/fish/high verts", 1)) %>%
    mutate(Diet.highverts=replace(Diet.highverts, diet=="NA", NA))%>%
    
    dplyr::select(-diet)
  
  # diet.to.infer <- traits.raw %>%
  #     dplyr::select(Species, diet) 
  #     diet.to.infer$diet <- as.factor(diet.to.infer$diet)
  #     diet.to.infer$diet<-fct_relevel(diet.to.infer$diet, c("plankton", "inverts", "fish", "inverts/fish", 
  #                                                           "inverts/fish/high verts", "NA"))
  #     diet.to.infer $diet <- as.numeric(diet.to.infer$diet)
  #     diet.to.infer $diet[diet.to.infer$diet==6] <- NA
  #     
  
# ordinal traits
  
  # terrestriality<- traits.new %>%
    terrestriality<- traits.curated %>%
    dplyr::select(Species, terrestriality) 
    terrestriality$terrestriality <- as.factor(terrestriality$terrestriality)
    terrestriality$terrestriality<-fct_relevel(terrestriality$terrestriality, c("no","brackish","freshwater", "NA"))
    terrestriality$terrestriality<-as.numeric(terrestriality$terrestriality)
    terrestriality$terrestriality[terrestriality$terrestriality==4]<-NA
  
  # thermo<- traits.raw %>%
    thermo<- traits.curated %>%
    dplyr::select(Species, thermo)
    thermo$thermo <- as.factor(thermo$thermo)
    thermo$thermo<-fct_relevel(thermo$thermo, c("ecto","meso"))
    thermo$thermo<-as.numeric(thermo$thermo)
  
  # feeding<- traits.raw %>%
    feeding<- traits.curated %>%
    dplyr::select(Species, feeding) 
    feeding$feeding <- as.factor(feeding$feeding)
    feeding$feeding<-fct_relevel(feeding$feeding, c("filter feeder", "macropredator"))
    feeding$feeding<-as.numeric(feeding$feeding)
    
# IUCN
  iucn <- as_tibble(read.csv( "~/Dropbox (Smithsonian)/SHARKS/IUCN/species_updated.iucn_no.syn.csv")) %>%
    rename(iucn = category) %>%
    dplyr::select(Species, iucn) %>%
    right_join(traits.raw[,2], by="Species") %>%
    mutate(iucn=replace(iucn, is.na(iucn), "NE")) %>%
    mutate(iucn=replace(iucn, iucn=="CR",4))%>%
    mutate(iucn=replace(iucn, iucn=="EN",3))%>%  
    mutate(iucn=replace(iucn, iucn=="VU",2))%>%
    mutate(iucn=replace(iucn, iucn=="NT",1))%>%
    mutate(iucn=replace(iucn, iucn=="LC",0))%>%
    mutate(iucn=replace(iucn, iucn=="DD",NA))%>%
    mutate(iucn=replace(iucn, iucn=="NE",NA))
  iucn$iucn <- as.numeric(iucn$iucn)
  
  setdiff(iucn$Species, traits.raw$Species); setdiff(traits.raw$Species, iucn$Species)
  which(duplicated(traits.raw$Species)==T)
  
 #no.status <- traits.raw[setdiff(traits.raw$Species, iucn$Species),6]

# add biomes
  regions_as.traits <- as_tibble(read.csv("~/Dropbox (Smithsonian)/SHARKS/SHARKS-TR/Sharks FD/Analyses/R_outputs/all_species_regions.csv")) %>%
    dplyr::select(-X) %>%
    relocate(Species)
  
  # add species without sptial data but inferred from text in Fishbase (not using it as it would add a different treatment)
  # sp.to.add <- traits.raw[,c(2, 21:27)] %>%
  #   relocate(Species, Arctic, `North Atlantic`, `North Pacific`, Indian, `South Pacific`, `South Atlantic`, Southern) %>%
  #   setNames(names(regions_as.traits)) %>%
  #   filter(Species %in% setdiff(traits.raw$Species, regions_as.traits$Species))
  
  regions_as.traits <- regions_as.traits %>%
    right_join(traits.raw[,2], by="Species")
    # inner_join(traits.raw[,2], by="Species") %>%
    # bind_rows(sp.to.add)

  setdiff(regions_as.traits$Species, traits.raw$Species); setdiff(traits.raw$Species, regions_as.traits$Species)
  setdiff(habitat$Species, regions_as.traits$Species)
  
# add geographical ranges
  range_as.traits <- as_tibble(read.csv("~/Dropbox (Smithsonian)/SHARKS/SHARKS-TR/Sharks FD/Analyses/R_outputs/all_species_range.csv")) %>%
    dplyr::select(-X)%>%
    right_join(traits.raw[,2], by="Species")
    
  # setdiff(range_as.traits$Species, regions_as.traits$Species)
    
# combine all
shark.traits <- traits.new  %>%
  dplyr::select(Species, max.length) %>%
  left_join(orders, by="Species") %>%
  left_join(families, by="Species") %>%
  left_join(habitat, by="Species") %>%
  left_join(vertical, by="Species") %>%
  left_join(diet, by="Species") %>%
  left_join(terrestriality, by="Species") %>%
  left_join(thermo, by="Species") %>%
  left_join(feeding, by="Species") %>%
  left_join(iucn, by="Species") %>%
  left_join(regions_as.traits, by = "Species") %>%
  left_join(range_as.traits, by= "Species") %>%
  dplyr::rename(size = max.length)

setdiff(traits.raw$Species, shark.traits$Species); setdiff(shark.traits$Species,traits.raw$Species)

# save
# write.csv(shark.traits, file="~/Dropbox (Smithsonian)/SHARKS/SHARKS-TR/Sharks FD/Analyses/Imputations/Oct_2021/shark.traits.imputations_all.csv")
# save using one taxonomy as one-hot encoding 
 write.csv(shark.traits, file="~/Dropbox (Smithsonian)/SHARKS/SHARKS-TR/Sharks FD/Analyses/Imputations/Oct_2021/shark.traits.imputations_all2.csv")
 
  
# # For the purposes of imputing diet:
#   shark.traits.base <- shark.traits %>%
#      dplyr::select(-Diet.plankton, -Diet.inverts, -Diet.fish, -Diet.highverts) %>%
#      left_join(diet.to.infer, by="Species")
#   # base table with complete cases only 
#     shark.traits.to.infer <-  na.omit(shark.traits.base)
#     write.csv(shark.traits.to.infer, file="~/Dropbox (Smithsonian)/SHARKS/SHARKS-TR/Sharks FD/Analyses/Diet inferences/shark.traits.base.csv") 
#  # table with traits needed to be imputed, with NAs
#    shark.traits.imp <-  shark.traits.base %>%
#      filter(is.na(diet))
#     write.csv(shark.traits.imp, file="~/Dropbox (Smithsonian)/SHARKS/SHARKS-TR/Sharks FD/Analyses/Diet inferences/shark.traits.imp.csv")
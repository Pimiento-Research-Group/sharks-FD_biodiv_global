# This script does the following:
# 1. assigns biomes/regions to shark species from IUCN data
# 2. Transform grids into occurrences to then calculaye per-species geographic range
# Date: 16 Jul 2021
# Author: Catalina Pimiento

library(readxl)
library(stringr)
library(tidyverse)
library(dplyr)
library(janitor)
library(reshape2)
library (raster)

#########################################################################################################
# REGIONS
#########################################################################################################
# read look up table with all the shark names and synonyms
  names.raw<-read_xlsx("~/Dropbox (Smithsonian)/SHARKS/Taxonomy/Lookup_Taxonomy.xlsx", sheet = "Species")

  names.simple<- names.raw %>%
    dplyr::select(-Family, -Order, -Superorder)

# load spatial data
  distr <- as_tibble(
    readRDS("~/Dropbox (Smithsonian)/SHARKS/SHARKS-TR/Sharks FD/Analyses/R_inputs/Map_PA_biomes_0_5.RDS"))
  
# load("~/Dropbox (Smithsonian)/SHARKS/SHARKS-TR/Sharks FD/Analyses/Mat_Pa_CHONDRICHTHYES.rdata")
# distr2<-as_tibble(Mat_Pa_CHONDRICHTHYES)

# load data frame with fixed synonyms
  load('~/Dropbox (Smithsonian)/SHARKS/SHARKS-TR/Sharks FD/Analyses/R_outputs/iucn.names_traits.names.RData') # distr.names

# replace names in distr with valid names (there will be duplicates)
  new.names <- c("lon","lat", distr.names$new.names)
  length(new.names)
  which(duplicated(new.names)==TRUE)

  distr.no.syn<- distr %>%
    rename(lon = X_ENTIER, lat = Y_ENTIER) %>%
    dplyr:: select(-Region, -Biome, -Biome1) %>%
    setNames(new.names)

# merge occurrence values from duplicated names (merge junior synonyms)
  distr.no.syn <- 
    as_tibble(do.call(cbind,lapply(split(seq_len(ncol(distr.no.syn)),names(distr.no.syn)),
                                   function(x) rowSums(distr.no.syn[x]))))
  
  which(duplicated(names(distr.no.syn))==TRUE)

# re-oragnise and re-convert to 0s and 1s and bind biomes
  grids <- as_tibble(distr.no.syn) %>%
    dplyr::select(-lon,-lat) %>%
    replace(. > 1, 1) %>%
    bind_cols(lon= distr.no.syn$lon, lat=distr.no.syn$lat, biomes=distr$Biome1) %>%
    dplyr::relocate(lon, lat, biomes) 
  
  save(grids,file=
         '~/Dropbox (Smithsonian)/SHARKS/SHARKS-TR/Sharks FD/Analyses/R_outputs/grids_biomes_no.syn.RData') 
  
  grids <- na.omit(grids) %>%
    dplyr::select(-lon, -lat)

# for loop to extract biomes of each species
  species_names <- names(grids)[2:dim(grids)[2]]
  biome_names <- unique(grids$biomes)
  biomes_matrix <- matrix(nrow = length(biome_names), ncol = length(species_names))
  
  for (i in 1:length(species_names)) {
    biomes_matrix[, i] = 0
    biomes_sp_i = unique(grids[which(grids[,i+1]==1),1])$biomes
    biomes_matrix[match(biomes_sp_i, biome_names), i] = 1
  }
  
  regions <- as.data.frame(t(biomes_matrix))
  colnames(regions) <- biome_names
  regions$Species <- species_names
  
  setdiff(names.simple$Species, regions$Species)
  
  write.csv(regions,
            "~/Dropbox (Smithsonian)/SHARKS/SHARKS-TR/Sharks FD/Analyses/R_outputs/all_species_regions.csv")
   
   no_spatial_data <- setdiff(names.raw$Species, regions$Species)[c(1, 11, 39, 41:42,61)] 
   
  #######################################################################################################
  # GEOGRAPHICAL RANGE
  #######################################################################################################
   # load grids and biomes without synonyms
     load('~/Dropbox (Smithsonian)/SHARKS/SHARKS-TR/Sharks FD/Analyses/R_outputs/grids_biomes_no.syn.RData')

   # convert grids into occurrences
   #   occ <- grids %>%
   #     unite("grid", lon:lat, remove = FALSE) %>%
   #     dplyr::select(-lat, -lon, -biomes)
   # 
   #   occ_melted <- melt(occ) # only run this with a clean environment, otherwise memory will be exhausted
   #   occ_melted <- occ_melted[occ_melted$value != 0, ]
   # 
   #   occ_melted$decimallongitude <-  sapply(strsplit(occ_melted$grid, "_"),"[",1)
   #   occ_melted$decimallatitude <- sapply(strsplit(occ_melted$grid, "_"),"[",2)
   # 
   # save(occ_melted,file='~/Dropbox (Smithsonian)/SHARKS/IUCN/occurrences_iucn_0.5.RData')
   
   load("~/Dropbox (Smithsonian)/SHARKS/IUCN/occurrences_iucn_0.5.RData")
   
   iucn.occ <- as_tibble(occ_melted) %>%
     dplyr::select(-grid, -value) %>%
     rename(Species=variable) %>%
     mutate(decimallongitude = as.numeric(decimallongitude),
            decimallatitude = as.numeric(decimallatitude))
   
   range <-  iucn.occ %>%
     group_by(Species) %>%
     summarise(min_long = min(decimallongitude), 
               max_long = max(decimallongitude),
               mid_long = median(decimallongitude),
               
               min_lat = min(decimallatitude),
               max_lat = max(decimallatitude), 
               mid_lat = median(decimallatitude))
   
write.csv(range, "~/Dropbox (Smithsonian)/SHARKS/SHARKS-TR/Sharks FD/Analyses/R_outputs/all_species_range.csv")
   
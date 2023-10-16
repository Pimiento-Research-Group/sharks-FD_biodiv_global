# This script:
# 1. Creates a data frame with species names from IUCN and new names given according with trait dataset
# 2. Uses grid data from IUCN at a 0.5 resolution
# Date: 16 Jul 2021
# Author: Catalina Pimiento


library(readxl)
library(stringr)
library(tidyverse)
library(dplyr)
library(janitor)
library(reshape2)

#  read look up table with all the shark names and synonyms
   names.raw<-read_xlsx("~/Dropbox (Smithsonian)/SHARKS/Taxonomy/Lookup_Taxonomy.xlsx", sheet = "Species")

   which(duplicated(names.raw$Species)==T)

   names.simple<- names.raw %>%
      dplyr::select(-Family, -Order, -Superorder)

   length(names.simple$Species)

#  load spatial data
   distr <- as_tibble(
      readRDS("~/Dropbox (Smithsonian)/SHARKS/SHARKS-TR/Sharks FD/Analyses/R_inputs/Map_PA_biomes_0_5.RDS"))

# delete chimeras    
   chim1 <- str_subset(names(distr), "Callorhinchus")
   chim2 <- str_subset(names(distr), "Chimaera")
   chim3 <- str_subset(names(distr), "Hydrolagus")
   chim4 <- str_subset(names(distr), "Harriotta")

   distr <- distr %>%
      dplyr::select(-chim1, -chim2, -chim3, -chim4)

# FIX SYNONYMS
   # data frame with names in distr data and their corresponding valid name
      distr.names <- tibble(names(distr)) %>%
         setNames("Species") %>%
         filter(Species!= "X_ENTIER" & Species!="Y_ENTIER" & 
                   Species!="Region" & Species!="Biome"& Species!="Biome1") %>%
         left_join(names.simple, by = c("Species" = "Synonyms.1")) %>%
         left_join(names.simple, by = c("Species" = "Synonyms.2")) %>%
         left_join(names.simple, by = c("Species" = "Synonyms.3")) %>%
         left_join(names.simple, by = c("Species" = "Synonyms.4")) %>%
         left_join(names.simple, by = c("Species" = "Synonyms.5")) %>%
         left_join(names.simple, by = c("Species" = "Synonyms.6")) %>%
         left_join(names.simple, by = c("Species" = "Synonyms.7")) %>%
         left_join(names.simple, by = c("Species" = "Synonyms.8")) %>%
         left_join(names.simple, by = c("Species" = "Synonyms.9")) %>%
         left_join(names.simple, by = c("Species" = "Synonyms.10")) %>%
         remove_empty("cols") %>%
         dplyr::select(-contains("Synonyms")) %>%
         mutate(new.names = dplyr::coalesce(Species.y,Species.y.y, Species.y.y.y, 
                                            Species.y.y.y.y, Species.y.y.y.y.y,
                                            Species.y.y.y.y.y.y, Species.y.y.y.y.y.y.y,
                                            Species.y.y.y.y.y.y.y.y,Species.y.y.y.y.y.y.y.y.y, 
                                            Species.y.y.y.y.y.y.y.y.y.y, Species)) %>%
         dplyr::select(-contains("y"))

   length(names(distr))-2 # species from original distr data: 1114
   length(unique(distr.names$new.names)) # corresponding valid names: 1091

   save(distr.names,file=
           '~/Dropbox (Smithsonian)/SHARKS/SHARKS-TR/Sharks FD/Analyses/R_outputs/iucn.names_traits.names.RData')

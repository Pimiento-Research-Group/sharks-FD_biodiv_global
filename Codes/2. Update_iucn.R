# This script updates IUCN updated_iucn of sharks and rays
# Date: 15 Jul 2021
# Author: Catalina Pimiento

# load libraries
  library(dplyr)
  library(tidyverse)
  library(readxl)
  library(janitor)
  library(rredlist)
  library(purrr)
  
# # using rredlist requires an API token
#   Sys.setenv(IUCN_KEY = "d13b7c9bbc6d1c2e2099771acf4ab79c1b1a45080bb1e5b2e4b9b57b1706665a")
#   Sys.getenv("IUCN_KEY")
#   apikey <- Sys.getenv("IUCN_KEY") # needed to access rredlist package
# 
# # species names downloaded from iucn for spatial analyses
#   load("~/Dropbox (Smithsonian)/SHARKS/SHARKS-TR/Sharks FD/Analyses/R_outputs/iucn.names_traits.names.RData")
# 
# # get iucn info on those species names
#   species_iucn<-distr.names %>%
#     dplyr::select(Species) %>% # now apply the rl_search function to each species using purrr::map()
#     mutate(iucn_pull = map(Species, rl_search, key = apikey))
# 
# # clean dataset and gather only iucn statuses
#   species_iucn_clean <-  species_iucn %>%
#     mutate(category = map_chr(iucn_pull, pluck, "result", "category", .default = NA)) %>%
#     dplyr::select(Species, category) %>%
#     bind_cols(distr.names$new.names) %>%
#     rename(trait.names=...3)
# 
# write.csv(species_iucn_clean, "~/Dropbox (Smithsonian)/SHARKS/IUCN/iucn_cat_syn.names.csv")

species_iucn_clean <- read.csv("~/Dropbox (Smithsonian)/SHARKS/IUCN/iucn_cat_syn.names.csv")
  
# address synonyms and use valid names, which are also used in trait analyses
  # identify duplicated names
    synonyms <- species_iucn_clean %>%
      group_by(trait.names) %>%
      filter(n()>1)
    
  # identify synonyms that did not result in duplicated names
    which(is.na(match(species_iucn_clean$Species, species_iucn_clean$trait.names))==T)
    
  # eliminate duplicated and leave status of matching species names
    synonyms <-  synonyms[row.names(synonyms) %in% 
                            unique(match(synonyms$trait.names, synonyms$Species)),]
   
   species_iucn_updated <- species_iucn_clean %>%
     filter(!trait.names %in% synonyms$trait.names) %>%
     bind_rows(synonyms) %>%
     dplyr::select(-Species)
   
   which(duplicated(species_iucn_updated$trait.names)==T) #no duplicates
   
  # species which names changed in iucn website when we downloaded the data
   wrong.names <- species_iucn_updated %>%
     dplyr:: filter(is.na(category))
   
   species_iucn_updated_fixed <- species_iucn_updated %>%
     mutate(category=replace(category, trait.names=="Dipturus australis", "NT")) %>%
     mutate(category=replace(category, trait.names=="Dipturus cerva", "NT")) %>%
     mutate(category=replace(category, trait.names=="Dipturus confusus", "CR")) %>%
     mutate(category=replace(category, trait.names=="Dipturus endeavouri", "NT")) %>%
     mutate(category=replace(category, trait.names=="Dipturus grahami", "LC")) %>%
     mutate(category=replace(category, trait.names=="Dipturus healdi", "LC")) %>%
     mutate(category=replace(category, trait.names=="Zearaja maugeana", "EN")) %>%
     mutate(category=replace(category, trait.names=="Dipturus oculus", "LC")) %>%
     mutate(category=replace(category, trait.names=="Dipturus polyommata", "LC")) %>%
     mutate(category=replace(category, trait.names=="Himantura oxyrhyncha", "EN")) %>%
     mutate(category=replace(category, trait.names=="Myliobatis freminvillei", "VU")) %>%
     mutate(category=replace(category, trait.names=="Narcine lasti", "LC")) %>%
     mutate(category=replace(category, trait.names=="Narcine nelsoni", "LC")) %>%
     mutate(category=replace(category, trait.names=="Narcine ornata", "LC")) %>%
     mutate(category=replace(category, trait.names=="Narcine tasmaniensis", "LC")) %>%
     mutate(category=replace(category, trait.names=="Narcine westraliensis", "LC")) %>%
     mutate(category=replace(category, trait.names=="Psammobatis parvacauda", "LC")) %>%
     mutate(category=replace(category, trait.names=="Raja inornata", "LC")) %>%
     mutate(category=replace(category, trait.names=="Raja stellulata", "LC")) %>%
     mutate(category=replace(category, trait.names=="Taeniura grabata", "NT"))
   
# add missing species from trait dataset
  # load look up table with synonyms
    all.names <-read_xlsx("~/Dropbox (Smithsonian)/SHARKS/Taxonomy/Lookup_Taxonomy.xlsx", 
                          sheet = "Species") %>%
      dplyr::select(-Family, -Order, -Superorder) 
    
  # what species are not already in distr data?
    missing.sp <- all.names %>%
      filter(!Species %in% species_iucn_updated_fixed$trait.names) 
  
        syn.1 =  missing.sp %>%
          dplyr::select(Synonyms.1) %>%
          dplyr::filter(!Synonyms.1=="NA") %>%
          dplyr::rename(Species = Synonyms.1)
         syn.2 =  missing.sp %>%
           dplyr::select(Synonyms.2) %>%
           dplyr::filter(!Synonyms.2=="NA")%>%
           dplyr::rename(Species = Synonyms.2)
         # syn.3 =  missing.sp %>%
         #   dplyr::select(Synonyms.3) %>%
         #   dplyr::filter(!Synonyms.3=="NA")%>%
         #   dplyr::rename(Species = Synonyms.3)
         # syn.4 =  missing.sp %>%
         #   dplyr::select(Synonyms.4) %>%
         #   dplyr::filter(!Synonyms.4=="NA")%>%
         #   dplyr::rename(Species = Synonyms.4)

# check the iucn category of these species that are missing and their synonyms
  missing.sp.syn <- missing.sp %>%
    dplyr::select(Species) %>%
    bind_rows(syn.1,syn.2) %>%
    dplyr::select(Species) %>% # now apply the rl_search function to each species using purrr::map()
    mutate(iucn_pull = map(Species, rl_search, key = apikey)) 
  
  missing.sp.syn_clean <- missing.sp.syn %>%
    mutate(category = map_chr(iucn_pull, pluck, "result", "category", .default = NA)) %>% 
    dplyr::select(Species, category) %>% 
    dplyr::filter(!is.na(category)) %>%
    dplyr::select(category, Species)
  
 intersect(missing.sp.syn_clean$Species, missing.sp$Species)
 
 species_iucn_updated_no.syn <- species_iucn_updated_fixed %>%
   rename(Species = trait.names) %>%
   bind_rows(missing.sp.syn_clean)
 
# get still missing species
  setdiff(all.names$Species, species_iucn_updated_no.syn$Species) 
  # these are all not evaluated and are considered junior synonyms of other species but they are not according to Weigmann
  
# add not evaluated species
  species_iucn_updated_no.syn <- species_iucn_updated_no.syn %>%
    bind_rows(
      tibble(setdiff(all.names$Species,species_iucn_updated_no.syn$Species)) %>%
    setNames("Species") %>%
    mutate(category="NE")
    )
  
  write.csv(species_iucn_updated_no.syn, 
            "~/Dropbox (Smithsonian)/SHARKS/IUCN/species_updated.iucn_no.syn.csv")
  
  
# get their synonyms from iucn
 missing.syns <- as_tibble(setdiff(all.names$Species, species_iucn_updated$Species)) %>%
   setNames("Species")%>%
   mutate(iucn_pull = map(Species, rl_synonyms, key = apikey))
        
 missing.syns_clean <- missing.syns %>% 
   mutate(accepted_name = map_chr(iucn_pull, pluck, "result", "accepted_name", .default = NA)) %>% 
   select(Species, accepted_name) %>%
   filter(!is.na(accepted_name))
   
 
################################################################################################################################
# checking based on new publication by Dulvy
 
 pre <- read_csv("~/Dropbox (Smithsonian)/SHARKS/IUCN/species_updated.iucn_no.syn_pre.csv")   
 species_iucn_updated_no.syn<- read_csv("~/Dropbox (Smithsonian)/SHARKS/IUCN/species_updated.iucn_no.syn.csv")   
 
 identical(pre$category, species_iucn_updated_no.syn$category)
 which(ifelse(pre$category==species_iucn_updated_no.syn$category,"Yes","No")=="No")
 
 updated_names <- read_csv("~/Dropbox (Smithsonian)/SHARKS/IUCN/Dulvy_2021_species.csv")   
 
 setdiff(updated_names$Species, species_iucn_updated_no.syn$Species)
 
################################################################################################################################
# old way
#
# load datasets
  setwd("~/Dropbox (Smithsonian)/SHARKS/Traits")
  traits<-read_xlsx("shark_traits.xlsx")
  
  setwd("~/Dropbox (Smithsonian)/SHARKS/SHARKS-TR/Sharks FD/Analyses/redlist_species_data_eed15753-4bac-4841-957c-e853652020d2/")
  
  iucn<-as_tibble(read.csv("assessments.csv")) %>%
    mutate(Species=as.character(scientificName)) %>%
    mutate(updated_iucn=as.character(status))%>%
    dplyr::select(Species, updated_iucn)
  
   # delete chimeras
     chim1 <- str_subset(iucn$Species, "Callorhinchus")
     chim2 <- str_subset(iucn$Species, "Chimaera")
     chim3 <- str_subset(iucn$Species, "Hydrolagus")
     chim4 <- str_subset(iucn$Species, "Harriotta")
     chim5 <- str_subset(iucn$Species, "Neoharriotta")
     chim6 <- str_subset(iucn$Species, "Rhinochimaera")

    iucn <- iucn %>%
      filter(!Species %in% chim1 & !Species%in%chim2 & !Species%in%chim3 &
             !Species%in%chim4 & !Species%in%chim5 & !Species%in%chim6)

  # check for inconsistentcies  
    setdiff(iucn$Species, traits$Species)
    setdiff(traits$Species, iucn$Species)

  # check if due to synonyms
    setwd("~/Dropbox (Smithsonian)/SHARKS/Taxonomy/")
    names<-read_xlsx("Lookup_Taxonomy.xlsx", sheet = "Species")%>%
      dplyr::select(-Family, -Order, -Superorder)

    intersect(iucn$Species, names$Synonyms.1)# yes
    intersect(iucn$Species, names$Synonyms.2)# yes
    intersect(iucn$Species, names$Synonyms.3)# yes
    intersect(iucn$Species, names$Synonyms.4)# yes
    intersect(iucn$Species, names$Synonyms.5)# yes
    intersect(iucn$Species, names$Synonyms.6)# yes
    intersect(iucn$Species, names$Synonyms.7)# no  
    
iucn.names <- 
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
    mutate(mycol = dplyr::coalesce(Species.y,Species.y.y, Species.y.y.y, Species.y.y.y.y, Species.y.y.y.y.y,
                                   Species.y.y.y.y.y.y, Species.y.y.y.y.y.y.y,Species.y.y.y.y.y.y.y.y,
                                   Species.y.y.y.y.y.y.y.y.y, Species.y.y.y.y.y.y.y.y.y.y, Species)) %>%
    mutate(Species=mycol) %>%
    dplyr::select(-contains("y"), -mycol) %>%
    distinct(Species,.keep_all = TRUE) # eliminate duplicates created because of synonyms
    
  # check for inconsistentcies  
    setdiff(traits$Species, iucn.names$Species)
    setdiff(iucn.names$Species, traits$Species)
    
  # save
    write.csv(iucn.names, file="~/Dropbox (Smithsonian)/SHARKS/SHARKS-TR/Sharks FD/Analyses/iucn.updated.csv")

    
  

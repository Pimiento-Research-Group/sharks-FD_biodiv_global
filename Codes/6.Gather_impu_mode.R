# This script does the following:
# 1. Checks proportion of missing data per trait
# 2. Computes modal value across imputations
# 3. Computes correlations between inferred values across imputations
# Date: 7 Dec 2021
# Authors: Catalina Pimiento and Daniele Silvestro


library(readxl)
library(stringr)
library(tidyverse)
library(dplyr)
library(ggcorrplot)
library(cowplot)
library(writexl)

# WHICH TRAITS WERE IMPUTED?
  shark.traits <- read_csv("~/Dropbox (Smithsonian)/SHARKS/SHARKS-TR/Sharks FD/Analyses/Imputations/Oct_2021/shark.traits.imputations_all2.csv") 
  traits.curated<-read_xlsx("~/Dropbox (Smithsonian)/SHARKS/SHARKS-TR/Sharks FD/Analyses/R_inputs/shark_traits_curated.xlsx") 
  
# which had NAs?
  imp.data <- 
    data.frame(nas=rbind(
      length(which(traits.curated$habitat=="NA")), length(which(traits.curated$vertical=="NA")),
      length(which(traits.curated$terrestriality=="NA")), length(which(traits.curated$thermo=="NA")),
      length(which(traits.curated$feeding=="NA")), length(which(traits.curated$diet=="NA")),
      length(which(traits.curated$max.length=="NA")), length(which(is.na(shark.traits$iucn))),
      length(which(is.na(shark.traits$Arctic_Ocean)))))
  
  rownames(imp.data) <- c("habitat", "vertical", "terrestriality", "thermo","feeding","diet","size",
                          "iucn", "oceans")

  imp.data
  which(imp.data>0)
  
# GET THE MODE ACROSS IMPUTATIONS
# using imputations that used traits, IUCN and geographic info together
# get imputed data: joint diet and IUCN  
  load("~/Dropbox (Smithsonian)/SHARKS/SHARKS-TR/Sharks FD/Analyses/Imputations/Oct_2021/shark_imputations_1.rda")
  
  # gather all data imputed from each trait
    habitat.c <- res_all %>% 
      map(., function(x) as.numeric(as.character(x[[1]][["Habitat.coastal"]])))
    habitat.o <- res_all %>% 
      map(., function(x) as.numeric(as.character(x[[1]][["Habitat.oceanic"]])))
    terrestriality <- res_all %>% 
      map(., function(x) as.numeric(as.character(x[[1]][["terrestriality"]])))
    diet.p <- res_all %>% 
      map(., function(x) as.numeric(as.character(x[[1]][["Diet.plankton"]])))
    diet.i<- res_all %>% 
      map(., function(x) as.numeric(as.character(x[[1]][["Diet.inverts"]]))) 
    diet.f <- res_all %>% 
      map(., function(x) as.numeric(as.character(x[[1]][["Diet.fish"]]))) 
    diet.h <- res_all %>% 
      map(., function(x) as.numeric(as.character(x[[1]][["Diet.highverts"]]))) 
    size <- res_all %>% 
      map(., function(x) as.numeric(as.character(x[[1]][["size"]])))
    iucn <- res_all %>% 
      map(., function(x) as.numeric(as.character(x[[1]][["iucn"]])))
    ocean.a <- res_all %>% 
      map(., function(x) as.numeric(as.character(x[[1]][["Arctic_Ocean"]])))
    ocean.na <- res_all %>% 
      map(., function(x) as.numeric(as.character(x[[1]][["North_Atlantic_Ocean"]])))
    ocean.np <- res_all %>% 
      map(., function(x) as.numeric(as.character(x[[1]][["North_Pacific_Ocean"]])))
    ocean.i <- res_all %>% 
      map(., function(x) as.numeric(as.character(x[[1]][["Indian_Ocean"]])))
    ocean.sp <- res_all %>% 
      map(., function(x) as.numeric(as.character(x[[1]][["South_Pacific_Ocean"]])))
    ocean.sa <- res_all %>% 
      map(., function(x) as.numeric(as.character(x[[1]][["South_Atlantic_Ocean"]])))
    ocean.s <- res_all %>% 
      map(., function(x) as.numeric(as.character(x[[1]][["Southern.Ocean"]])))

# mode function
   mymode <- function(x) {
    t <- table(x)
    names(t)[ which.max(t) ]
  }

# compute the mode
  habitat.c.mode <- as_tibble(do.call(cbind,habitat.c)) %>%
    rowwise()%>%
    mutate(mode=mymode(c(V1, V2, V3, V4, V5, V6, V7, V8, V9, V10)))
  
  habitat.o.mode <- as_tibble(do.call(cbind,habitat.o)) %>%
    rowwise()%>%
    mutate(mode=mymode(c(V1, V2, V3, V4, V5, V6, V7, V8, V9, V10)))
  
  terrestriality.mode <- as_tibble(do.call(cbind,terrestriality)) %>%
    rowwise()%>%
    mutate(mode=mymode(c(V1, V2, V3, V4, V5, V6, V7, V8, V9, V10)))
  
  diet.p.mode <- as_tibble(do.call(cbind,diet.p)) %>%
    rowwise()%>%
    mutate(mode=mymode(c(V1, V2, V3, V4, V5, V6, V7, V8, V9, V10)))
  
  diet.i.mode <- as_tibble(do.call(cbind,diet.i)) %>%
    rowwise()%>%
    mutate(mode=mymode(c(V1, V2, V3, V4, V5, V6, V7, V8, V9, V10)))
  
  diet.f.mode <- as_tibble(do.call(cbind,diet.f)) %>%
    rowwise()%>%
    mutate(mode=mymode(c(V1, V2, V3, V4, V5, V6, V7, V8, V9, V10)))
  
  diet.h.mode <- as_tibble(do.call(cbind,diet.h)) %>%
    rowwise()%>%
    mutate(mode=mymode(c(V1, V2, V3, V4, V5, V6, V7, V8, V9, V10)))
  
  size.mode <- as_tibble(do.call(cbind,size)) %>%
    rowwise()%>%
    mutate(mode=mymode(c(V1, V2, V3, V4, V5, V6, V7, V8, V9, V10)))
  
  iucn.mode <- as_tibble(do.call(cbind,iucn)) %>%
    rowwise()%>%
    mutate(mode=mymode(c(V1, V2, V3, V4, V5, V6, V7, V8, V9, V10)))
  
  ocean.a.mode <- as_tibble(do.call(cbind,ocean.a)) %>%
    rowwise()%>%
    mutate(mode=mymode(c(V1, V2, V3, V4, V5, V6, V7, V8, V9, V10)))
  
  ocean.na.mode <- as_tibble(do.call(cbind,ocean.na)) %>%
    rowwise()%>%
    mutate(mode=mymode(c(V1, V2, V3, V4, V5, V6, V7, V8, V9, V10)))
  
  ocean.np.mode <- as_tibble(do.call(cbind,ocean.np)) %>%
    rowwise()%>%
    mutate(mode=mymode(c(V1, V2, V3, V4, V5, V6, V7, V8, V9, V10)))
  
  ocean.i.mode <- as_tibble(do.call(cbind,ocean.i)) %>%
    rowwise()%>%
    mutate(mode=mymode(c(V1, V2, V3, V4, V5, V6, V7, V8, V9, V10)))
  
  ocean.sp.mode <- as_tibble(do.call(cbind,ocean.sp)) %>%
    rowwise()%>%
    mutate(mode=mymode(c(V1, V2, V3, V4, V5, V6, V7, V8, V9, V10)))
  
  ocean.sa.mode <- as_tibble(do.call(cbind,ocean.sa)) %>%
    rowwise()%>%
    mutate(mode=mymode(c(V1, V2, V3, V4, V5, V6, V7, V8, V9, V10)))
  
  ocean.s.mode <- as_tibble(do.call(cbind,ocean.s)) %>%
    rowwise()%>%
    mutate(mode=mymode(c(V1, V2, V3, V4, V5, V6, V7, V8, V9, V10)))
  
# put everything together in a tibble with only the mode
  
  imp.mode <- shark.traits %>%
    dplyr::select(Species) %>%
      bind_cols(habitat.c.mode$mode, habitat.o.mode$mode, terrestriality.mode$mode,
              diet.p.mode$mode, diet.i.mode$mode, diet.f.mode$mode, diet.h.mode$mode,
              size.mode$mode, iucn.mode$mode,ocean.a.mode$mode, ocean.na.mode$mode,
              ocean.np.mode$mode, ocean.i.mode$mode,ocean.sp.mode$mode, ocean.sa.mode$mode,
              ocean.s.mode$mode) %>%
      setNames(c("Species", "habitat.coastal", "habitat.oceanic", "terrestriality",
               "diet.plankton", "diet.inverts", "diet.fish", "diet.highverts",
               "size", "iucn", "arctic", "north.atlantic","north.pacific", "indian",
               "south.pacific", "south.atlantic", "southern"))
          
  
  save(imp.mode,file='~/Dropbox (Smithsonian)/SHARKS/SHARKS-TR/Sharks FD/Analyses/R_outputs/imp.mode.RData')
  write_xlsx(imp.mode, "~/Dropbox (Smithsonian)/SHARKS/SHARKS-TR/Sharks FD/Analyses/R_outputs/imp.mode.xlsx")
  

# put everything together to calculate fuse across all imputations and obtain uncertaintity
  imp.all <- shark.traits %>%
    dplyr::select(Species) %>%
    bind_cols(habitat.c.mode %>%
                dplyr::select(V1:V10), 
              habitat.o.mode%>%
                dplyr::select(V1:V10), 
              terrestriality.mode%>%
                dplyr::select(V1:V10),
              diet.p.mode%>%
                dplyr::select(V1:V10), 
              diet.i.mode%>%
                dplyr::select(V1:V10), 
              diet.f.mode%>%
                dplyr::select(V1:V10), 
              diet.h.mode%>%
                dplyr::select(V1:V10),
              size.mode%>%
                dplyr::select(V1:V10), 
              iucn.mode%>%
                dplyr::select(V1:V10)) %>%
    setNames(c("Species", 
               seq(1, 10) %>%
                 paste0("habitat.coastal", "_", .),
               seq(1, 10) %>%
                 paste0("habitat.oceanic", "_", .), 
               seq(1, 10) %>%
                 paste0("terrestriality", "_", .), 
               seq(1, 10) %>%
                 paste0("diet.plankton", "_", .), 
               seq(1, 10) %>%
                 paste0("diet.inverts", "_", .), 
               seq(1, 10) %>%
                 paste0("diet.fish", "_", .), 
               seq(1, 10) %>%
                 paste0("diet.highverts", "_", .), 
               seq(1, 10) %>%
                 paste0("size", "_", .), 
               seq(1, 10) %>%
                 paste0("iucn", "_", .)))
    
  save(imp.all,file='~/Dropbox (Smithsonian)/SHARKS/SHARKS-TR/Sharks FD/Analyses/R_outputs/imp.all.RData')
  
  plot_grid(
    
 imp.all %>%
  dplyr::select(matches("habitat.coastal")) %>%
  cor(use="pairwise.complete.obs", method = "spearman") %>% 
  ggcorrplot(show.diag=FALSE, type="lower", lab=TRUE, lab_size=2),

 imp.all %>%
   dplyr::select(matches("habitat.oceanic")) %>%
   cor(use="pairwise.complete.obs", method = "spearman") %>% 
   ggcorrplot(show.diag=FALSE, type="lower", lab=TRUE, lab_size=2),
 
 imp.all %>%
   dplyr::select(matches("terrestriality")) %>%
   cor(use="pairwise.complete.obs", method = "spearman") %>% 
   ggcorrplot(show.diag=FALSE, type="lower", lab=TRUE, lab_size=2),
 
 imp.all %>%
   dplyr::select(matches("diet.plankton")) %>%
   cor(use="pairwise.complete.obs", method = "spearman") %>% 
   ggcorrplot(show.diag=FALSE, type="lower", lab=TRUE, lab_size=2),
 
 imp.all %>%
   dplyr::select(matches("diet.inverts")) %>%
   cor(use="pairwise.complete.obs", method = "spearman") %>% 
   ggcorrplot(show.diag=FALSE, type="lower", lab=TRUE, lab_size=2),
 
 imp.all %>%
   dplyr::select(matches("diet.fish")) %>%
   cor(use="pairwise.complete.obs", method = "spearman") %>% 
   ggcorrplot(show.diag=FALSE, type="lower", lab=TRUE, lab_size=2),
 
 imp.all %>%
   dplyr::select(matches("diet.highverts")) %>%
   cor(use="pairwise.complete.obs", method = "spearman") %>% 
   ggcorrplot(show.diag=FALSE, type="lower", lab=TRUE, lab_size=2),
 
 imp.all %>%
   dplyr::select(matches("size")) %>%
   cor(use="pairwise.complete.obs") %>% 
   ggcorrplot(show.diag=FALSE, type="lower", lab=TRUE, lab_size=2),
 
 imp.all %>%
   dplyr::select(matches("iucn")) %>%
   cor(use="pairwise.complete.obs", method = "spearman") %>% 
   ggcorrplot(show.diag=FALSE, type="lower", lab=TRUE, lab_size=2),
 
 ncol=3)
 
####################################################################################################################################

# uisng imputations that used traits and IUCN separately
  load("~/Dropbox (Smithsonian)/SHARKS/SHARKS-TR/Sharks FD/Analyses/Imputations/Oct_2021/shark_diet_imputations_1.rda")
  load("~/Dropbox (Smithsonian)/SHARKS/SHARKS-TR/Sharks FD/Analyses/Imputations/Oct_2021/shark_iucn_imputations_1.rda")
  
  
# gather all data imputed from each trait
  habitat.c <- res_diet %>% 
    map(., function(x) as.numeric(as.character(x[[1]][["Habitat.coastal"]])))
  habitat.o <- res_diet %>% 
    map(., function(x) as.numeric(as.character(x[[1]][["Habitat.oceanic"]])))
  terrestriality <- res_diet %>% 
    map(., function(x) as.numeric(as.character(x[[1]][["terrestriality"]])))
  diet.p <- res_diet %>% 
    map(., function(x) as.numeric(as.character(x[[1]][["Diet.plankton"]])))
  diet.i<- res_diet %>% 
    map(., function(x) as.numeric(as.character(x[[1]][["Diet.inverts"]]))) 
  diet.f <- res_diet %>% 
    map(., function(x) as.numeric(as.character(x[[1]][["Diet.fish"]]))) 
  diet.h <- res_diet %>% 
    map(., function(x) as.numeric(as.character(x[[1]][["Diet.highverts"]]))) 
  size <- res_diet %>% 
    map(., function(x) as.numeric(as.character(x[[1]][["size"]])))
  iucn <- res_iucn %>% 
    map(., function(x) as.numeric(as.character(x[[1]][["iucn"]])))
  ocean.a <- res_iucn %>% 
    map(., function(x) as.numeric(as.character(x[[1]][["Arctic_Ocean"]])))
  ocean.na <- res_iucn %>% 
    map(., function(x) as.numeric(as.character(x[[1]][["North_Atlantic_Ocean"]])))
  ocean.np <- res_iucn %>% 
    map(., function(x) as.numeric(as.character(x[[1]][["North_Pacific_Ocean"]])))
  ocean.i <- res_iucn %>% 
    map(., function(x) as.numeric(as.character(x[[1]][["Indian_Ocean"]])))
  ocean.sp <- res_iucn %>% 
    map(., function(x) as.numeric(as.character(x[[1]][["South_Pacific_Ocean"]])))
  ocean.sa <- res_iucn %>% 
    map(., function(x) as.numeric(as.character(x[[1]][["South_Atlantic_Ocean"]])))
  ocean.s <- res_iucn %>% 
    map(., function(x) as.numeric(as.character(x[[1]][["Southern.Ocean"]])))
  
 # mode function
  mymode <- function(x) {
    t <- table(x)
    names(t)[ which.max(t) ]
  }
  
# compute the mode
  habitat.c.mode <- as_tibble(do.call(cbind,habitat.c)) %>%
    rowwise()%>%
    mutate(mode=mymode(c(V1, V2, V3, V4, V5, V6, V7, V8, V9, V10)))
  
  habitat.o.mode <- as_tibble(do.call(cbind,habitat.o)) %>%
    rowwise()%>%
    mutate(mode=mymode(c(V1, V2, V3, V4, V5, V6, V7, V8, V9, V10)))
  
  terrestriality.mode <- as_tibble(do.call(cbind,terrestriality)) %>%
    rowwise()%>%
    mutate(mode=mymode(c(V1, V2, V3, V4, V5, V6, V7, V8, V9, V10)))
  
  diet.p.mode <- as_tibble(do.call(cbind,diet.p)) %>%
    rowwise()%>%
    mutate(mode=mymode(c(V1, V2, V3, V4, V5, V6, V7, V8, V9, V10)))
  
  diet.i.mode <- as_tibble(do.call(cbind,diet.i)) %>%
    rowwise()%>%
    mutate(mode=mymode(c(V1, V2, V3, V4, V5, V6, V7, V8, V9, V10)))
  
  diet.f.mode <- as_tibble(do.call(cbind,diet.f)) %>%
    rowwise()%>%
    mutate(mode=mymode(c(V1, V2, V3, V4, V5, V6, V7, V8, V9, V10)))
  
  diet.h.mode <- as_tibble(do.call(cbind,diet.h)) %>%
    rowwise()%>%
    mutate(mode=mymode(c(V1, V2, V3, V4, V5, V6, V7, V8, V9, V10)))
  
  size.mode <- as_tibble(do.call(cbind,size)) %>%
    rowwise()%>%
    mutate(mode=mymode(c(V1, V2, V3, V4, V5, V6, V7, V8, V9, V10)))
  
  iucn.mode <- as_tibble(do.call(cbind,iucn)) %>%
    rowwise()%>%
    mutate(mode=mymode(c(V1, V2, V3, V4, V5, V6, V7, V8, V9, V10)))
  
  ocean.a.mode <- as_tibble(do.call(cbind,ocean.a)) %>%
    rowwise()%>%
    mutate(mode=mymode(c(V1, V2, V3, V4, V5, V6, V7, V8, V9, V10)))
  
  ocean.na.mode <- as_tibble(do.call(cbind,ocean.na)) %>%
    rowwise()%>%
    mutate(mode=mymode(c(V1, V2, V3, V4, V5, V6, V7, V8, V9, V10)))
  
  ocean.np.mode <- as_tibble(do.call(cbind,ocean.np)) %>%
    rowwise()%>%
    mutate(mode=mymode(c(V1, V2, V3, V4, V5, V6, V7, V8, V9, V10)))
  
  ocean.i.mode <- as_tibble(do.call(cbind,ocean.i)) %>%
    rowwise()%>%
    mutate(mode=mymode(c(V1, V2, V3, V4, V5, V6, V7, V8, V9, V10)))
  
  ocean.sp.mode <- as_tibble(do.call(cbind,ocean.sp)) %>%
    rowwise()%>%
    mutate(mode=mymode(c(V1, V2, V3, V4, V5, V6, V7, V8, V9, V10)))
  
  ocean.sa.mode <- as_tibble(do.call(cbind,ocean.sa)) %>%
    rowwise()%>%
    mutate(mode=mymode(c(V1, V2, V3, V4, V5, V6, V7, V8, V9, V10)))
  
  ocean.s.mode <- as_tibble(do.call(cbind,ocean.s)) %>%
    rowwise()%>%
    mutate(mode=mymode(c(V1, V2, V3, V4, V5, V6, V7, V8, V9, V10)))
  
# put everything together in a tibble
  imp.mode_2 <- shark.traits %>%
    dplyr::select(Species) %>%
    bind_cols(habitat.c.mode$mode, habitat.o.mode$mode, terrestriality.mode$mode,
              diet.p.mode$mode, diet.i.mode$mode, diet.f.mode$mode, diet.h.mode$mode,
              size.mode$mode, iucn.mode$mode,ocean.a.mode$mode, ocean.na.mode$mode,
              ocean.np.mode$mode, ocean.i.mode$mode,ocean.sp.mode$mode, ocean.sa.mode$mode,
              ocean.s.mode$mode) %>%
    setNames(c("Species", "habitat.coastal", "habitat.oceanic", "terrestriality",
               "diet.plankton", "diet.inverts", "diet.fish", "diet.highverts",
               "size", "iucn", "arctic", "north.atlantic","north.pacific", "indian",
               "south.pacific", "south.atlantic", "southern"))
  
  
  save(imp.mode_2,file='~/Dropbox (Smithsonian)/SHARKS/SHARKS-TR/Sharks FD/Analyses/R_outputs/imp.mode_2.RData')
  
####################################################################################################################################
####################################################################################################################################
# check differences in trait values
# visual inspection

diet_imp <- imp.mode_2 %>%
    dplyr::select(Species, diet.plankton, diet.inverts, diet.fish, diet.highverts) %>%
    left_join(
      imp.mode %>%
                dplyr::select(Species, diet.plankton, diet.inverts, diet.fish, diet.highverts), by="Species" ) %>%
    pivot_longer(diet.plankton.x:diet.highverts.y, names_to = "trait", values_to = "diet") %>%
    filter(diet==1)%>% 
    mutate(imp = ifelse(trait %in% c("diet.plankton.x","diet.inverts.x","diet.fish.x","diet.highverts.x"),"imp1","imp2")) %>%
    select(-diet)
  
  ggplot(diet_imp,
         aes(y=Species,x=trait, col=imp))+
    geom_point()
  
  ggsave("~/Dropbox (Smithsonian)/SHARKS/SHARKS-TR/Sharks FD/Analyses/R_outputs/Plots/imputations_diet_diff.jpg",
         plot = last_plot(),width = 30,height = 60, units = c("cm"), dpi = 300)

  iucn_imp <- imp.mode_2 %>%
    dplyr::select(Species, iucn) %>%
    left_join(
      imp.mode %>%
        dplyr::select(Species, iucn), by="Species") %>%
    pivot_longer(iucn.x:iucn.y, names_to = "trait", values_to = "iucn")%>%
    mutate(imp = ifelse(trait %in% "iucn.x","imp1","imp2")) %>%
    select(-trait)
    
  ggplot(iucn_imp,
         aes(y=Species,x=iucn, col=imp))+
    geom_point()+
    facet_wrap(vars(imp))
  
  ggsave("~/Dropbox (Smithsonian)/SHARKS/SHARKS-TR/Sharks FD/Analyses/R_outputs/Plots/imputations_iucn_diff.jpg",
         plot = last_plot(),width = 30,height = 60, units = c("cm"), dpi = 300)
  

# differences
 
  iucn_diff <- imp.mode_2 %>%
    dplyr::select(Species, iucn) %>%
    left_join(
      imp.mode %>%
        dplyr::select(Species, iucn), by="Species") %>%
    mutate(dif=as.numeric(iucn.x)-as.numeric(iucn.y)) %>%
    filter(dif!=0)
  
  (dim(iucn_diff)[1]*100)/dim(imp.mode_2)[1]
  
  setdiff(iucn_diff$iucn.x,iucn_diff$iucn.y)
  
  diet_diff <- imp.mode_2 %>%
    dplyr::select(Species, diet.plankton, diet.inverts, diet.fish, diet.highverts) %>%
    left_join(
      imp.mode %>%
        dplyr::select(Species, diet.plankton, diet.inverts, diet.fish, diet.highverts), by="Species")%>%
    mutate(plankton=as.numeric(diet.plankton.x)-as.numeric(diet.plankton.y),
           inverts=as.numeric(diet.inverts.x)-as.numeric(diet.inverts.y),
           fish=as.numeric(diet.fish.x)-as.numeric(diet.fish.y),
           highverts=as.numeric(diet.highverts.x)-as.numeric(diet.highverts.y))
  
  length(which(diet_diff$plankton!=0))*100/dim(imp.mode_2)[1]
  length(which(diet_diff$inverts!=0))*100/dim(imp.mode_2)[1]
  length(which(diet_diff$fish!=0))*100/dim(imp.mode_2)[1]
  length(which(diet_diff$highverts!=0))*100/dim(imp.mode_2)[1]
  
           
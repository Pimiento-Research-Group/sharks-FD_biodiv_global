## R code for the analysis of: Elasmobranch functional diversity is highly vulnerable and supported by unique species and locations worldwide


### Authors

Catalina Pimiento, Camille Albouy, Daniele Silvestro, Theophile Mouton, Laure Velez, David Mouillot, Aaron B. Judah, John N. Griffin, Fabien Leprieur

### Intro

This folder contains scripts for the FUNCTIONAL DIVERSITY analysis of data in R. Required packages: `tidyverse`, `dplyr`, `tidyr`, `readxl`,`stringr`, `purrr`, `readxl`, `tidyselect`, `data.table`, `reshape2` and `janitor` for data handling and manipulation; `ade4`,`geometry`,`raster`,`missForest`,`mFD` and `vegan` for analyses; and `cowplot`, `grid`, `visreg`, `ggcorrplot`,`wesanderson`,`gghighligh` and `RColorBrewer` for plotting.

These can be downloaded and installed using the following commands:

```
install.packages("data.table")
install.packages("dplyr")
install.packages("tidyr")
install.packages("reshape2")
install.packages("ade4")
install.packages("geometry")
install.packages("readxl")
install.packages("stringr")
install.packages("tidyverse")
install.packages("janitor")
install.packages("rredlist")
install.packages("purrr")
install.package("raster")
install.package("missForest")
install.package("ggcorrplot")
install.package("cowplot")
install.package("beepr")
install.package("rgdal")
install.package("tidyselect") 
install.package("mFD")
install.package("vegan")
install.package("cowplot")
install.package("grid")
install.package("visreg")
install.package("ggcorrplot")
install.package("wesanderson") 
install.package("gghighlight")
install.package("RColorBrewer)"
install.package("wesanderson")

```

### Codes
The scripts are as follows
1.Fix_Synonyms.R: Identifies discrepancies in species names due to synonyms between distribution data (downloaded from IUCN) and trait data. Produces:
  - iucn.names_traits.names.RData
2. Update_iucn.R: Connects to IUCN API and gathers IUCN status. Produces: 
  - iucn_cat_syn.names.csv (dataframe with trait and IUCN names and their updated status)
  - species_updated.iucn_no.syn.csv  (trait names and updated IUCN)
3. Regions_Ouccurrences.R: Transform gridded data (from 0 above) into occurrences to then calculate per-species geographic range. Produces: 
  - grids_biomes_no.syn.RData
  - all_species_regions.csv
  - all_species_range.csv
4. Data.prep.imputations.R: It prepares the data for multiple imputations and transform taxonomy into one-hot-ecoding. Produces:
  - shark.traits.imputations_all2.csv
5. Imputations.R: Performs multiple imputations *using one-hot-encoding* and calculates error. Produces:
  - shark_diet_imputations_1.rda
  - shark_iucn_imputations_1.rda
  - shark_imputations_1.rda
6. Gather_impu_mode.R: Checks proportion of missing data per trait, computes modal value across imputations and correlations between inferred values across imputations. Produces:
  - imp.mode.RData
  - imp.all.RData
7. FD-Analyses.Rmd: Performs all functional diversity analyses and produces all figures and tables cited in the paper. It requires the following functions:
  - quality_funct_space_fromdist.R
  - fonction_FRIC_Global_full.R
  - get_indicator_function 2.R
  
All codes and data for the phylogenetic and spatial analyses are provided in separate folders as follows:

Phylogenetic analyses
  - The phylo_analysis folder is related to all the analyses and data for quantifying species and assemblage-level phylogenetic metrics. 
    It produces 2 objects used in FD-Analyses.Rmd
    - Res_HED_sharks.rds
    - Res_HEDGE_sharks.rds
MPA_analysis
  - The MPA_analysis folder includes all the information to quantify and map MPAs
Spatial analyses
 - The congruence_mapping folder includes all the analyses made to identify and map the hotspots and to quantify the congruence among metrics (including fishing pressure).




    
    
    
    
    
    
    
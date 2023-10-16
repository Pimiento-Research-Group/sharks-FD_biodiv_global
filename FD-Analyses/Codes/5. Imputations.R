# Trait imputations
# Author: Daniele Silvestro
# Date: Oct 2021

setwd("~/Dropbox (Smithsonian)/SHARKS/SHARKS-TR/Sharks FD/Analyses/Imputations/Oct_2021")

one_hot =1

if (one_hot == 1){
    tbl_file = "shark.traits.imputations_all2.csv"    
} else{
    tbl_file = "shark.traits.imputations_all.csv"
}


library(missForest)
tbl <- read.csv(tbl_file, row.names=1)
head(tbl)

# prep data types
red_data <- tbl
dtypes = sapply(red_data, class)
red_data$size <- log10(red_data$size)

red_data[,grep("_long",colnames(tbl))] <- red_data[,grep("_long",colnames(tbl))]/100
red_data[,grep("_lat",colnames(tbl))] <- red_data[,grep("_lat",colnames(tbl))]/100
# turn categorical to factors
ind = which(sapply(red_data, typeof) == "integer")
for ( i in ind){
    red_data[,i] = factor(red_data[,i])
}

if (one_hot == 0){
    # remove family data
    red_data <- red_data[, -grep("Family", colnames(red_data))]
}

### IMPUTE *DIET* USING MISSFOREST (based on all traits excluding IUCN, ocean, long/lat)

# drop columns
red_data_diet <- subset(red_data, select=Species:feeding)

res = NULL
n_imputations = 10
for (i in 1:n_imputations){
    set.seed(i)
    print(c("Running imputation n.", i))
    #drop first column (with species names) 
    #run missForest
    missForest_imputation <- missForest(xmis = red_data_diet[- c(1)], 
                                        maxiter = 100, 
                                        ntree = 100, 
                                        variablewise = TRUE)
    head(missForest_imputation$ximp) #missForest printed results
    res[[i]] <- missForest_imputation
}

for (i in 1:n_imputations){
    r <- res[[i]]$ximp
    # rescale body size
    r$size <- 10^(r$size)
    res[[i]]$ximp <- r
}
res_diet <- res
save(res_diet, file=paste0("shark_diet_imputations_", one_hot, ".rda"))


### IMPUTE *IUCN* USING MISSFOREST (based on IUCN, oceean, long/lat, and phylo info)

# drop columns
red_data_iucn <- subset(red_data, select=c(Species, iucn:mid_lat))
red_data_iucn <- cbind(red_data_iucn, red_data[,grep('Order', colnames(red_data))])
red_data_iucn <- cbind(red_data_iucn, red_data[,grep('Family', colnames(red_data))])
res = NULL
n_imputations = 10
for (i in 1:n_imputations){
    set.seed(i)
    print(c("Running imputation n.", i))
    #drop first column (with species names) 
    #run missForest
    missForest_imputation <- missForest(xmis = red_data_iucn[- c(1)], 
                                        maxiter = 100, 
                                        ntree = 100, 
                                        variablewise = TRUE)
    head(missForest_imputation$ximp) #missForest printed results
    res[[i]] <- missForest_imputation
}

for (i in 1:n_imputations){
    r <- res[[i]]$ximp
    # rescale back longitude
    r[,grep("_long",colnames(r))] <- r[,grep("_long",colnames(r))] * 100
    # rescale back longitude
    r[,grep("_lat",colnames(r))] <- r[,grep("_lat",colnames(r))] * 180
    # fix long/lat column names
    colnames(r) <- gsub("_long", "LAT", colnames(r))
    colnames(r) <- gsub("_lat", "_longitude", colnames(r))
    colnames(r) <- gsub("LAT", "_latitude", colnames(r))
    res[[i]]$ximp <- r
}
res_iucn <- res
save(res_iucn, file=paste0("shark_iucn_imputations_", one_hot, ".rda"))


### IMPUTE USING MISSFOREST (joint imputation of IUCN status and traits)
res = NULL
n_imputations = 10
for (i in 1:n_imputations){
    set.seed(i)
    print(c("Running imputation n.", i))
    #drop first column (with species names) 
    #run missForest
    missForest_imputation <- missForest(xmis = red_data[- c(1)], 
                                        maxiter = 100, 
                                        ntree = 100, 
                                        variablewise = TRUE)
    head(missForest_imputation$ximp) #missForest printed results
    res[[i]] <- missForest_imputation
}

for (i in 1:n_imputations){
    r <- res[[i]]$ximp
    # rescale back longitude
    r[,grep("_long",colnames(r))] <- r[,grep("_long",colnames(r))] * 180
    # rescale back longitude
    r[,grep("_lat",colnames(r))] <- r[,grep("_lat",colnames(r))] * 100
    # rescale body size
    r$size <- 10^(r$size)
    res[[i]]$ximp <- r
}

res_all <- res
save(res_all, file=paste0("shark_imputations_", one_hot, ".rda"))


# CHECK ERROR

one_hot = 1 

load(paste0("shark_diet_imputations_", one_hot, ".rda"))
load(paste0("shark_iucn_imputations_", one_hot, ".rda"))
load(paste0("shark_imputations_", one_hot, ".rda"))


# quantify OOB error
getOOBerror <- function(res, colname, n_imputations = 10){
    OOBs = c()
    for (i in 1:n_imputations){
        indx = grep(colname, colnames(res[[i]]$ximp))
        OOBs <- c(OOBs, as.numeric(res[[i]]$OOBerror[indx]))
    }
    return(OOBs)
}

mean(getOOBerror(res_all, "Diet"))*100 
mean(getOOBerror(res_diet, "Diet"))*100 # DIET* USING MISSFOREST (based on all traits excluding IUCN, ocean, long/lat)

mean(getOOBerror(res_all, "iucn"))*100
mean(getOOBerror(res_iucn, "iucn"))*100 # IUCN* USING MISSFOREST (based on IUCN, oceean, long/lat, and phylo info)




as.numeric(as.data.frame(res[[1]]$ximp)) - as.numeric(as.data.frame(res[[2]]$ximp))


f = 1 - (as.data.frame(res[[3]]$ximp) == as.data.frame(res[[1]]$ximp))

apply(f, FUN=sum, 1)

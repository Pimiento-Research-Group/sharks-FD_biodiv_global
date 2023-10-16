####################################################################

get_dist_func <- function(nb_NN,data=s[[1]]){
  
  data <- data[order(data[,3], decreasing=F),]
  data <- data[-1,]
  mm <- mean(data[1:nb_NN,3])
  sd <- sd(data[1:nb_NN,3])
  sp <- as.character(data[1:nb_NN,1])
  
  return(list(c(Mean=mm,Sd=sd),Species=sp))
}

get_indicator <- function(Mat_dist = traits.mat.list[[1]],nb_NN=nb_NN){
  w <- melt(Mat_dist)
  s<- split(w,f=w[,2])
  
  Res <- lapply(s,function(x) {get_dist_func(nb_NN=nb_NN,data=x)})
  
  Res_mean_sd <- do.call(rbind,lapply(1:length(Res),function(i){Res[[i]][[1]]}))
  rownames(Res_mean_sd) <- names(Res)
  
  NN <- lapply(1:length(Res),function(i){Res[[i]][[2]]})
  names(NN) <- names(Res)
  
  return(list(Average_uniqueness=Res_mean_sd,Nearest_neighbour=NN))
}
  
###################################################################################
  
get_FU_region <- function(Mat_dist = traits.mat.list[[1]], nb_NN, Selected_sp){

Region_Mat_dist<-Mat_dist[Selected_sp,Selected_sp]

FU_sp<-get_indicator(Region_Mat_dist, nb_NN)

FU_Region_Mean<-mean(FU_sp$Average_uniqueness[,"Mean"])
FU_Region_sd<-sd(FU_sp$Average_uniqueness[,"Sd"])

data.frame(FU_Region_Mean=FU_Region_Mean,FU_Region_sd=FU_Region_sd)}

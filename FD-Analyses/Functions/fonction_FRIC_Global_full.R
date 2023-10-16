gb_fric_function <- function(pcoa=pcoa, ax=ax){
  coord_d<-pcoa$li[,ax]
  rs <- nrow(coord_d)
  
  # functional Richness
  FRic <- round(convhulln(coord_d,"FA")$vol,10)
  
  # identity of vertices
  vert0<-convhulln(coord_d,"Fx TO 'vert.txt'")
  vert1<-scan("vert.txt",quiet=T)
  vert2<-(vert1+1)[-1]
  
  FE_vert_k<- row.names(coord_d)[vert2] 
  FE_vertices_cells <- FE_vert_k[!is.na(FE_vert_k)] # arendre
  
  # number of vertices
  nbFE_vertices<-length(FE_vert_k)
  
  ### Calcul du FDiv
  # coord values of vertices
  trvertices<-coord_d[vert2,]
  
  # coordinates of the center of gravity of the vertices (B)
  B<-apply(trvertices,2,function(x) {mean(x, na.rm=TRUE)})
  
  # computing euclidian dstances to B (dB)
  dB<-apply(coord_d, 1, function(x) { (sum((x-B)^2) )^0.5} )
  
  # mean of dB values, deviations to mean 
  meandB<-mean(dB)
  devdB<-dB-meandB
  
  # computation of FDiv
  FDiv<-round((sum(devdB)+meandB) / (sum(abs(devdB))+meandB) ,6)# arendre
  
  #computation of Originality
  dist_sp<-as.matrix(dist(coord_d,method="euclidean")) ; dist_sp[which(dist_sp==0)]<-NA
  orig_sp<-apply(dist_sp, 1, min, na.rm=T )
  Orig_vect <- c(mean=mean(orig_sp,na.rm=T), min=min(orig_sp, na.rm=T), max=max(orig_sp, na.rm=T), se=sd(orig_sp)/sqrt(length(orig_sp)))
  
  #computing functional specialization 
  O<-apply(coord_d, 2, mean)
  speS<-apply(coord_d, 1, function(x) { (sum((x-O)^2) )^0.5} )
  Specialization <- c(mean=mean(speS,na.rm=T), min=min(speS, na.rm=T), max=max(speS, na.rm=T), se=sd(speS)/sqrt(length(speS)))
  
  #computing FEs 
  #FEs <-unique(FE[,2])
  
  #computing FR
  #FR<- rs/length(FEs)
  
  #computing FV
  #r <- table(FE[,2])
  #FV <-(length(FEs)-sum(sapply(r,function(x){min(x-1,1)})))/length(FEs)
  
  list(data=data.frame(RS = rs, FRic = FRic, FDiv = FDiv,nb_vertices = nbFE_vertices), specialization= Specialization, originality=Orig_vect, Sp_vertice=FE_vertices_cells)
  
}





get_FV_Sp <- function(ax=ax,Fric_tot=richness_funct_global_ST[[1]][2],pcoa=pcoa,Selected_sp){

  coord_d<-pcoa$li[,ax]
  coord_d_k <- na.omit(coord_d[as.character(Selected_sp),])
  rs <- dim(coord_d_k)[1]

  # functional Richness
  FRic <- round(convhulln(coord_d_k,"FA")$vol,10)

  # identity of vertices
  vert0<-convhulln(coord_d_k,"Fx TO 'vert.txt'")
  vert1<-scan("vert.txt",quiet=T)
  vert2<-(vert1+1)[-1]

  FE_vert_k<- row.names(coord_d_k)[vert2]
  FE_vertices_cells <- FE_vert_k[!is.na(FE_vert_k)] # arendre

  # number of vertices
  nbFE_vertices<-length(FE_vert_k)

  ### Calcul du FDiv
  # coord values of vertices
  trvertices<-coord_d_k[vert2,]

  # coordinates of the center of gravity of the vertices (B)
  B<-apply(trvertices,2,function(x) {mean(x, na.rm=TRUE)})

  # computing euclidian dstances to B (dB)
  dB<-apply(coord_d_k, 1, function(x) { (sum((x-B)^2) )^0.5} )

  # mean of dB values, deviations to mean
  meandB<-mean(dB)
  devdB<-dB-meandB

  # computation of FDiv
  FDiv<-round((sum(devdB)+meandB) / (sum(abs(devdB))+meandB) ,6)# arendre

  # computation of FRic
  Fric_r <- round(FRic/convhulln(coord_d,"FA")$vol, 6)

  # computation of originality of each species: distance to nearest neighbour among the global pool of species
  dist_sp<-as.matrix(dist(coord_d_k,method="euclidean")) ; dist_sp[which(dist_sp==0)]<-NA
  orig_sp<-apply(dist_sp, 1, min, na.rm=T )
  Orig_vect <- c(mean=mean(orig_sp,na.rm=T), min=min(orig_sp, na.rm=T), max=max(orig_sp, na.rm=T), se=sd(orig_sp)/sqrt(length(orig_sp)))

  #computing functional specialization
  O<-apply(coord_d_k, 2, mean)
  speS<-apply(coord_d_k, 1, function(x) { (sum((x-O)^2) )^0.5} )
  Specialization <- c(mean=mean(speS,na.rm=T), min=min(speS, na.rm=T), max=max(speS, na.rm=T), se=sd(speS)/sqrt(length(speS)))

  #computing FEs
  #FEs <-unique(FE[Selected_sp,2])

  #computing FR
  #FR<- (length(Selected_sp))/length(FEs)

  #computing FV
  #r <- table(FE[Selected_sp,2])
  #FV <-(length(FEs)-sum(sapply(r,function(x){min(x-1,1)}))/length(FEs))

  list(data=data.frame(RS = rs, FRic = FRic, FDiv = FDiv,nb_vertices = nbFE_vertices, Fric_r=Fric_r), specialization= Specialization, originality=Orig_vect, Sp_vertice=FE_vertices_cells)
}  # end of function get_FV_Sp 





# mixted colors function 
biColorInter <- function(col1,col2,frac){
  
  #colors needs to be in RGB [0;1] space (will add for [0;255]
  #decomposing colors
  r1 <- col1[1]
  g1 <- col1[2]
  b1 <- col1[3]
  
  r2 <- col2[1]
  g2 <- col2[2]
  b2 <- col2[3]
  
  #deltas calculation
  Dr <- r2 - r1
  Dg <- g2 - g1
  Db <- b2 - b1
  
  #new r, g and b
  red = r1 + (Dr * frac)
  gre = g1 + (Dg * frac)
  blu = b1 + (Db * frac)
  
  unlist(c(red,gre,blu))
  
}#eo biColorInter

################### 4 colors interpolation
#c1---------c2
#
#
#c3---------c4

quadriColorInter <- function(colx1=roug,colx2=bleu,coly1=jaun,coly2=vert,fracx,fracy){
  
  finColX1 <- biColorInter(colx1,colx2,fracx)
  finColX2 <- biColorInter(coly1,coly2,fracx)
  fin <- biColorInter(finColX1,finColX2,fracy)
  rgb(fin[1],fin[2],fin[3])
  
}#eo quadriColorInter

#fonction pour borner une variable entre 0 et 1
getfunborn = function(x){
  #mini <- min(x)
  maxi <- max(x)
  #x <- (x - mini) / (maxi - mini)
  x <- x/maxi
  #x<-sqrt(x)
   return (x)
}

# fonction pour obtenir un vecteur de couleur à partir d'une condition d'un delta 
getdatacolor <- function (data=data_2MYA[,3:4],col1=jaunefabli, col2=roug, col3=bleu, col4=vert){
  
  data[,1] <- getfunborn (data[,1] )
  data[,2] <- getfunborn (data[,2] )
  
  col <- vector()
  
  for (i in 1:dim(data)[1]){
    col[i] <- quadriColorInter (col1,col2,col3,col4,fracx=data[i,1],fracy=data[i,2])  
  } # end of for 
  
  return(cbind(data,col))
  
} # end of getdatacolor()



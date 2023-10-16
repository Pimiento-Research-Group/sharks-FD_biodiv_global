#####################################################################################################
#                                                                                                   #
#                       Congruence analyses among biodiversity metrics & Mapping                    #
#                                                                                                   #
#                       Pimiento et al. Functional diversity of sharks and rays                     #
#                       is highly vulnerable and supported by uniquespecies and                     #  
#                       locations worldwide Submitted to Nature Communications August 2023          #
#                                                                                                   #             
#                         Authors : Camille Albouy (albouycamille@gmail.com)                        #
#                        fabien leprieur (fabien.leprieur@umontpellier.fr)                          #
#                                                                                                   #
#                                                                                                   #
##################################################################################################### 

# library loading 
lib_vect <- c("raster","gstat","rgeos","maptools","rgdal","maptools","parallel","RColorBrewer", "corrplot");
sapply(lib_vect,library,character.only=TRUE)

# set working directory

setwd(".../congruence_mapping") # your directory 

# Projection
ProjWGS84 <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"

# Loading world shape file
coast <- readShapePoly("data/GSHHS_l_L1_disolve.shp")

# Projection for mapping 
Proj <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs +towgs84=0,0,0"

# Functions loading 
source('functions/Conservation et services2_Hotspots_Congruence.R',chdir=T)
source("functions/Map_the_world.R")
source("functions/Graphic_functions.R")
source("functions/Function_carto_grad.R")
source("functions/Graphic_functions.R")
source("functions/Map_the_world_new.R")

###############################################################################
#                        Part 1 : Data preparation                            #
###############################################################################


# 1-1  Data loading
DataPD <- readRDS("data/sharks_phylo_results_june23.rds") # phylogenetic diversity metrics 
DataFD <-  readRDS("data/fd.metrics.grid.rds") # functional diversity metrics
Mat_PA <- readRDS("data/Mat_Pa.RDS") # species occurrence matrix (grid cells)
colnames(Mat_PA) <- gsub(" ","_",colnames(Mat_PA))

Coord <- apply(do.call(rbind,strsplit(x=rownames(Mat_PA),split="_")),2,as.numeric)
colnames(Coord) <- c("x","y")

Data_sp <- readRDS("data/sharks_iucn_final.rds") 
Treat_sp <- Data_sp[which(Data_sp$iucn_GE==2| Data_sp$iucn_GE==3 |Data_sp$iucn_GE==4),] # select species with CR, EN, VU status

# 1-2 Species richness (SR)  of threatened species

Mat_pa_treatsp <- Mat_PA[,Treat_sp$Species]
SR_treat <- apply(Mat_pa_treatsp,1,sum)


# 1-3 SR of the top 25 EDGE2 species
Data_TopEdge2 <- readRDS("data/Top_EDGE2_Q25_june2023.rds") 
Top_EDGE2_sp <- Mat_PA[,Data_TopEdge2$Species]
SR_Top_EDGE2 <- apply(Top_EDGE2_sp,1,sum)

# 1- 4 SR of the top 25 FUSE species
Data_TopFuse_100 <- readRDS("data/Top_FUSEP100_Q25_final_june2023.rds") 
Top_FUSE100_sp <- Mat_PA[,Data_TopFuse_100$Species]
SR_Top_FUSE100 <- apply(Top_FUSE100_sp,1,sum)

# 1- 5 SR of the top 25 ED2 species
Top_ED2_25 <- readRDS("data/Top_ED2_Q25_june2023.rds")
PA_Top_ED2_25 <- Mat_PA[,Top_ED2_25$Species]
SR_Top_ED2_25 <- apply(PA_Top_ED2_25,1,sum)

# 1-6 SR of the top 25 FUN species
Top_FUN25 <- readRDS("data/Top_FUn_Q25_final_june2023.rds")
PA_Top_FUN25 <- Mat_PA[,Top_FUN25$Species]
SR_Top_FUN25 <- apply(PA_Top_FUN25,1,sum)

# 1-7 SR of the top 25 FSp species
Top_FSp25 <- readRDS("data/Top_FSp_Q25_final_june2023.rds")
PA_Top_FSp25 <- Mat_PA[,Top_FSp25$Species]
SR_Top_FSp25 <- apply(PA_Top_FSp25,1,sum)

# 1-8 Residuals of the linear relationship between FRic and PD 

FD_PDres <- residuals(lm(DataFD$Fric~DataPD$PD))

# 1- 9 Residuals of the linear relationship between FRic and SR (species richness) 

FD_SRres <- residuals(lm(DataFD$Fric~DataFD$s))


# 1-10  Compiling biodiversity metrics 

Data_end_shark <- data.frame(Coord,SR=DataPD$RS,PD=DataPD$PD, PUn=DataPD$PUn,FRic=DataFD$Fric, FD_PDres=FD_PDres,FD_SRres,
                             FSp=DataFD$FSp,FUn=DataFD$FUn,ED=DataPD$ED_mean, ED2=DataPD$ED2_mean, SR_treat=SR_treat,SR_Top_EDGE2=SR_Top_EDGE2,
                             SR_Top_FUSE100=SR_Top_FUSE100,SR_Top_ED2_25=SR_Top_ED2_25, SR_Top_FUN25=SR_Top_FUN25, SR_Top_FSp25=SR_Top_FSp25)

saveRDS(Data_end_shark,"outputs/Data_end_shark_june2023.RDS")

### 1- 11 Fishing data loading

Fishing_data <- readRDS("data/Fishing_data.rds") # Fishing impact indices derived from the global fishing wath database : https://globalfishingwatch.org/datasets-and-code/
colnames(Fishing_data)[1:2] <- c("x","y")

## 1- 12 Krigging fishing data to match the spatial resolution of the 
Fishing_impact_2020 <- krige(Fishing_impact_2020~1,location=~x+y, data=data.frame(Fishing_data), newdata=data.frame(Coord),maxdist=3)
CummulFishing_impact_2012_2020 <- krige(CummulFishing_impact_2012_2020~1,location=~x+y, data=data.frame(Fishing_data), newdata=data.frame(Coord),maxdist=3)
Fishing_data_05<- data.frame(Coord,Fishing_impact_2020=Fishing_impact_2020[,3],CummulFishing_impact_2012_2020=CummulFishing_impact_2012_2020[,3])

#######################################################################
#                         Part 2 : Congruence analyses                #
#######################################################################

# 2-1 Congruence data 
Cong_data <- data.frame(x=DataPD$x,y=DataPD$y,RS=DataPD$RS,FRIC=DataFD$Fric, PD=DataPD$PD, FUn=DataFD$FUn1NN, 
                        PUn=DataPD$PUn,FSp=DataFD$FSp,SR_treat,SR_Top_EDGE2,SR_Top_FUSE100,SR_Top_ED_25, SR_Top_ED2_25,SR_Top_FUN25,
                        SR_Top_FSp25, ED2=DataPD$ED2_mean, ED=DataPD$ED_mean, Fishing_data_05[,-c(1:2)])


Cong_data <-  na.omit(Cong_data) 

Data_end <- Cong_data[,-c(1:2)] # Final data frame for the analyses

# 2-2 Congruence test using the the Cons_Serv_Permut() function

res_hot_spot <- Cons_Serv_Permut(Data_end,999,seuil=2.5/100) # using the 2.5% greatest values (2.5% hotspots)
res_hot_spot_5 <- Cons_Serv_Permut(Data_end,999,seuil=5/100) # using the 5% greatest values (5% hotspots)
res_hot_spot_10 <- Cons_Serv_Permut(Data_end,999,seuil=10/100) # using the 10% greatest values (10% hotspots)

saveRDS(res_hot_spot,file="outputs/res_hot_spot.RDS")
saveRDS(res_hot_spot_5,file="outputs/res_hot_spot_5.RDS")
saveRDS(res_hot_spot_10,file="outputs/res_hot_spot_10.RDS")

Hotspot <- cbind(res_hot_spot$Mat.Hotspot,X=Cong_data$x,Y=Cong_data$y)
Hotspot5 <- cbind(res_hot_spot_5$Mat.Hotspot,X=Cong_data$x,Y=Cong_data$y)
Hotspot10 <- cbind(res_hot_spot_10$Mat.Hotspot,X=Cong_data$x,Y=Cong_data$y)

###############################################################################                                                                           #
#                             Part 3  Mapping                                 #                                                                                                  #
###############################################################################

# 3- 1 Library loading
lib_vect <- c("raster","rgdal","maptools","parallel","RColorBrewer", "corrplot")
sapply(lib_vect,library,character.only=TRUE)


# 3-2 Colors definition 

col <- rev(c("#d73027","#fdae61","#fee090","#e0f3f8","#abd9e9","#6497E2"))
col1 <- rev(c("#d73027","#f46d43","#fdae61","#fee090","#ffffbf","#e0f3f8","#abd9e9","#74add1","#4575b4"))
colfunc <- colorRampPalette(c("#4575b4", "white", "#d73027"))
col2<-colfunc(16)

#############################
# 3- 3 : Figure 2 of the MS #
#############################

pdf("outputs/Figure_2_june2023.pdf",width=7.2,height=7.5)
d<-layout(mat=rbind(c(1,1,1,1,1,3,3,3,3,3),c(1,1,1,1,1,3,3,3,3,3),c(1,1,1,1,1,3,3,3,3,3),c(2,2,2,2,2,4,4,4,4,4),
                    c(5,5,5,5,5,7,7,7,7,7),c(5,5,5,5,5,7,7,7,7,7),c(5,5,5,5,5,7,7,7,7,7),c(6,6,6,6,6,8,8,8,8,8),
                    c(9,9,9,9,9,11,11,11,11,11),c(9,9,9,9,9,11,11,11,11,11),c(9,9,9,9,9,11,11,11,11,11),c(10,10,10,10,10,12,12,12,12,12)))

ylim=c(-53,73); xlim=c(-167,167)

### Plot FRic (A)

par(mar=c(0,2,2,0.5)); 

Map_the_world_3 (Data_PA=Data_end_shark$FRic,coord_X=Data_end_shark$x,coord_Y=Data_end_shark$y,
                 col=col1,breaks= seq(0,0.7,0.05),include.lowest=TRUE,xlim=xlim,ylim=ylim,coast=coast,names_fig="A",
                 cex_names_fig=0.7,x_names=0.5,y_names=87,cex.axis.leg=0.4,cex_point=0.25,
                 bg_col="gray80",bathy=FALSE,legend=FALSE,ncol=2,col_coast="gray20")
par(mar=c(0,2,0.2,0.5))
get_leg(data=Data_end_shark$FRic,breaks=seq(0,0.7,0.05),cex=0.6,Colour=col1,title="Functional Richness (FRic)",
        ncol=5,bty="n",pt.cex=1)


####### Plot Residuals of the linear relationship between FRic and SR (D)
par(mar=c(0.5,2,1.5,0.5))
Map_the_world_3 (Data_PA=FD_SRres,coord_X=Coord[,1],coord_Y=Coord[,2],
                 col=c("#173869","#4575B4","#87a3cc", "#fcfafa", "#EC9E9A", "#a1312b", "#821611"),breaks=c(seq(-0.5,0,0.1),seq(0.1, 0.3, 0.1)),include.lowest=T,xlim=xlim,ylim=ylim,coast=coast,names_fig="D",
                 cex_names_fig=0.7,x_names=0.5,y_names=87,cex.axis.leg=0.4,cex_point=0.25,legend=F,bg_col="gray80",col_coast="gray20")
par(mar=c(0,2,0.1,0.5))
get_leg(data=FD_SRres,breaks=c(seq(-0.5,0,0.1),seq(0.1,0.3,0.1)),cex=0.6,Colour=c("#173869","#4575B4","#87a3cc", "#fcfafa", "#EC9E9A", "#a1312b", "#821611"),
        title="Residuals of the linear relationship between FRic and SR",
        ncol=5,bty="n",pt.cex=1)


#### Plot FUn: mean nearest taxon distance based on functional distances = functional uniqueness (B)

par(mar=c(0.5,2,1.5,0.5))
Map_the_world_3 (Data_PA=Data_end_shark$FUn,coord_X=Coord[,1],coord_Y=Coord[,2],
                 col=col1,breaks=c(seq(0,0.23,0.01)),include.lowest=T,xlim=xlim,ylim=ylim,coast=coast,names_fig="B",
                 cex_names_fig=0.7,x_names=0.5,y_names=87,cex.axis.leg=0.4,cex_point=0.25,legend=F,bg_col="gray80",col_coast="gray20")
par(mar=c(0,2,0.1,0.5))
get_leg(data=DataFD$FUn,breaks=c(seq(0,0.23,0.01)),cex=0.6,Colour=col1,
        title="Functional Uniqueness (FUn)", ncol=5,bty="n",pt.cex=1)

#### Plot PUn : mean nearest taxon distance based on phylogenetic distances = phylogenetic uniqueness (E)

par(mar=c(0.5,2,1.5,0.5))
Map_the_world_3 (Data_PA=Data_end_shark$PUn,coord_X=Coord[,1],coord_Y=Coord[,2],
                 col=col1,breaks=c(seq(50,500,50)),include.lowest=T,xlim=xlim,ylim=ylim,coast=coast,names_fig="E",
                 cex_names_fig=0.7,x_names=0.5,y_names=87,cex.axis.leg=0.4,cex_point=0.25,legend=F,bg_col="gray80",col_coast="gray20")
par(mar=c(0,2,0.1,0.5))
get_leg(data=DataFD$FUn,breaks=c(seq(50,500,50)),cex=0.6,Colour=col1,
        title="Phylogenetic Uniqueness (PUn)", ncol=5,bty="n",pt.cex=1)

#### Functional Specialization (C)

par(mar=c(0.5,2,1.5,0.5))
Map_the_world_3 (Data_PA=Data_end_shark$FSp,coord_X=Coord[,1],coord_Y=Coord[,2],
                 col=col1,breaks=c(seq(0,0.23,0.02)),include.lowest=T,xlim=xlim,ylim=ylim,coast=coast,names_fig="C",
                 cex_names_fig=0.7,x_names=0.5,y_names=87,cex.axis.leg=0.4,cex_point=0.25,legend=F,bg_col="gray80",col_coast="gray20")
par(mar=c(0,2,0.1,0.5))
get_leg(data=DataFD$FUn,breaks=c(seq(0,0.23,0.02)),cex=0.6,Colour=col1,
        title="Functional Specialization (FSp)", ncol=5,bty="n",pt.cex=1)

#### Number of top 25% FUSE (F)

par(mar=c(0.5,2,1.5,0.5))
Map_the_world_3 (Data_PA=SR_Top_FUSE100,coord_X=Coord[,1],coord_Y=Coord[,2],
                 col=col1,breaks=c(seq(0,50,5)),include.lowest=T,xlim=xlim,ylim=ylim,coast=coast,names_fig="F",
                 cex_names_fig=0.7,x_names=0.5,y_names=87,cex.axis.leg=0.4,cex_point=0.25,legend=F,bg_col="gray80",col_coast="gray20")
par(mar=c(0,2,0.1,0.5))
get_leg(data=DataFD$FUn,breaks=c(seq(0,50, 5)),cex=0.6,Colour=col1,
        title="Number of top 25% FUSE species", ncol=5,bty="n",pt.cex=1)

dev.off()

###########################
# 3-4 Figure 3 of the MS. #
###########################

pdf("outputs/Figure_3_june2023.pdf",width=7.5,height=7.8)
layout(mat=rbind(c(1,1,1,1,1,2,2,2,2,2),c(1,1,1,1,1,2,2,2,2,2),c(1,1,1,1,1,2,2,2,2,2),c(1,1,1,1,1,2,2,2,2,2),
                 c(3,3,3,3,3,4,4,4,4,4),c(3,3,3,3,3,4,4,4,4,4),c(3,3,3,3,3,4,4,4,4,4),c(3,3,3,3,3,4,4,4,4,4),
                 c(5,5,5,5,5,6,6,6,6,6),c(5,5,5,5,5,6,6,6,6,6),c(5,5,5,5,5,6,6,6,6,6),c(5,5,5,5,5,6,6,6,6,6)))

# Hotspots of FRic (2.5% greatest value - hotspots) - A

par(mar=c(3,2,1.5,0.5))
carto_congruence_single (nbvar=1,data="Hotspot",x="FRIC",
                         names_fig="A",col_pt=c("#B70F0D"),cex_point=0.2,cex_text_fig=1,
                         cex_legend=0.6,main_legend="FRic hotspots",ncol=1,ylim=c(-53,73), xlim=c(-167,167),
                         cex.axis=0.4,axisx=T,axisy=T,legend=FALSE)

mtext("A",side=2,line=0.5,at=87,cex=0.7,bg="white",las=2)

# Congruence between EDGE2 and FUSE - D

par(mar=c(3,2,1.5,1))
carto_function(y=Hotspot,x="SR_Top_FUSE100", index2="SR_Top_EDGE2",names_fig="D",cex_point=0.2,cex_text_fig=0.7,cex_legend=0.6,
               text=c("Non-hotspot values","Values in top 2.5% for FUSE","Values in top 2.5% for EDGE","Congruence zone"),
               xlim=c(-167,167),ylim=c(-53,73),pos_legY=87,pos_legX=0.5,cex.axis=0.4,ncol=1,pos_leg="bottomleft",pt.cexleg=1,legend="")

# Congruence between FRic and SR - B

par(mar=c(3,2,1.5,0.5))
carto_function(y=Hotspot,x="RS",index2="FRIC",names_fig="B",cex_point=0.2,cex_text_fig=0.7,cex_legend=0.6,
               text=c("Non-hotspot values","Values in top 2.5% for FUSE","Values in top 2.5% for EDGE","Congruence zone"),
               xlim=c(-167,167),ylim=c(-53,73),pos_legY=87,pos_legX=0.5,cex.axis=0.4,ncol=1,pos_leg="bottomleft",pt.cexleg=1,legend="")

# Fishing hotspots - E

par(mar=c(3,2,1.5,1))
carto_congruence_single (nbvar=1,data="Hotspot",x="CummulFishing_impact_2012_2020",
                         names_fig="",col_pt=c("#B70F0D"),cex_point=0.28,cex_text_fig=1,cex_legend=0.6,
                         main_legend="FRic hotspots",ncol=1,ylim=c(-53,73), xlim=c(-167,167),
                         cex.axis=0.4,axisx=T,axisy=T,legend=F)

mtext("E",side=2,line=0.5,at=87,cex=0.7,bg="white",las=2)


# Congruence between FUn and PUn - C

par(mar=c(3,2,1.5,0.5))
carto_function(y=Hotspot,x="PUn",index2="FUn",names_fig="C",cex_point=0.2,cex_text_fig=0.7,cex_legend=0.6,
               text=c("Non-hotspot values","Values in top 2.5% for FUSE","Values in top 2.5% for EDGE","Congruence zone"),
               xlim=c(-167,167),ylim=c(-53,73),pos_legY=87,pos_legX=0.5,cex.axis=0.4,ncol=1,pos_leg="bottomleft",pt.cexleg=1,legend="C")


# Congruence between EDGE2,  FUSE and fishing impact - F

par(mar=c(3,2,1.5,1))
carto_congruence_t2.5(nbvar=3, data="Hotspot",x="SR_Top_EDGE2",index2="SR_Top_FUSE100",index3="CummulFishing_impact_2012_2020",
                      names_fig="",ylim=c(-53,73), xlim=c(-167,167), 
                      cex_point=0.27,cex_text_fig=0.4,cex_legend=0.45,
                      col_pt=c("#B70F0D"),pos_legY=86.5,pos_legX=-17.5,cex.axis=0.4,ncol=1, legend=FALSE)

mtext("F",side=2,line=0.5,at=87,cex=0.7,bg="white",las=2)


#### add the legend at the end
par(new=T, fig=c(0,1,0,1),mar=c(0,2,0.1,0.5))
plot(NULL,xaxt='n',yaxt='n',bty='n',ylab='',xlab='',xlim=0:1,ylim=0:1)

#A
legend(x=0 ,y=0.72,legend=c("No values","Non hotspots values","2.5% of FRic values"),
       pch=19,bty="n",col=c("#B6E0EE","#6BAACC","#B70F0D"),
       cex=0.65,horiz=T)
#D
legend(x=0.54,y=0.72,legend=c( text=c("No values","Non-hotspot values","Values in top 2.5% for EDGE2",
                                      "Values in top 2.5% for FUSE","Congruence zone")),
       pch=19,bty="n",col=c("#B6E0EE","#6BAACC","#B70F0D","#EDDF19","#fdae61"),
       cex=0.65,ncol=3)
#B

legend(x=0 ,y=0.355,legend=c( text=c("No values","Non-hotspot values","Values in top 2.5% for FRic",
                                     "Values in top 2.5% for SR","Congruence zone")),
       pch=19,bty="n",col=c("#B6E0EE","#6BAACC","#B70F0D","#EDDF19","#fdae61"),
       cex=0.65,ncol=3)
# E

legend(x=0.54,y=0.355,legend=c("No values","Non hotspots values","2.5% of fishing impact"),
       pch=19,bty="n",col=c("#B6E0EE","#6BAACC","red2"),
       cex=0.65,horiz=T)
# C

legend(x=0,y=0.001,legend=c( text=c("No values","Non-hotspot values","Values in top 2.5% for FUn",
                                    "Values in top 2.5% for PUn","Congruence zone")),
       pch=19,bty="n",col=c("#B6E0EE","#6BAACC","#B70F0D","#EDDF19","#fdae61"),
       cex=0.65,ncol=3)

# F

legend(x=0.52,y=0.001,legend=c("No values","Non hotspots values","2.5% of EDGE2, FUSE, Fishing"),
       pch=19,bty="n",col=c("#B6E0EE","#6BAACC","red2"),
       cex=0.65,horiz=T)


dev.off()


##############################
# 3-5 Figure S9 in Supp Mat  #        
##############################

pdf(file="outputs/FigureS9_june2023.pdf",width = 12, height = 7)

par(mfrow = c(2,3),mar=c(3,3,2,1))

smoothScatter(x=DataPD$RS,y=DataPD$PD,xlab="SR",ylab="PD",cex=0.1,pch=1,nrpoints=1,
              bg="darkblue",col="black",mgp=c(2,0.5,0))
mod1<-lm(DataPD$PD~DataPD$RS)
abline(mod1,lwd=2,col="red")
mtext("(a)",side=2,line=1.5,at=12400,cex=1,bg="white",las=2)    

smoothScatter(x=DataPD$RS,y=DataFD$Fric,xlab="SR",ylab="FRic",nrpoints=1,cex=0.1,pch=1,bg="darkblue",col="black",mgp=c(2,0.5,0))
mod2<-lm(DataFD$Fric~DataPD$RS)
abline(mod2,lwd=2,col="red")
mtext("(b)",side=2,line=1.2,at=0.72,cex=1,bg="white",las=2)  

smoothScatter(y=DataFD$Fric,x=DataPD$PD,ylab="FRic",xlab="PD",nrpoints=1,cex=0.15,pch=1,bg="darkblue",col="black",mgp=c(2,0.5,0))
mod3<-lm(DataFD$Fric~DataPD$PD)
abline(mod3,lwd=2,col="red")
mtext("(c)",side=2,line=1.2,at=0.72,cex=1,bg="white",las=2)

smoothScatter(y=DataFD$FUn1NN,x=DataFD$s,ylab="FUn",xlab="SR",nrpoints=1,cex=0.15,pch=1,bg="darkblue",col="black",mgp=c(2,0.5,0))
lines(loess.smooth(y=DataFD$FUn1NN,x=DataFD$s),lwd=2,col="red")
mtext("(d)",side=2,line=1.2,at=0.235,cex=1,bg="white",las=2)

smoothScatter(y=DataPD$PUn,x=DataPD$RS,ylab="PUn",xlab="SR",nrpoints=1,cex=0.15,pch=1,bg="darkblue",col="black",mgp=c(2,0.5,0))
lines(loess.smooth(y=DataPD$PUn,x=DataPD$RS),lwd=2,col="red")
mtext("(e)",side=2,line=1.2,at=520,cex=1,bg="white",las=2)

smoothScatter(y=DataFD$FSp,x=DataFD$s,ylab="FSp",xlab="SR",nrpoints=1,cex=0.15,pch=1,bg="darkblue",col="black",mgp=c(2,0.5,0))
lines(loess.smooth(y=DataFD$FSp,x=DataFD$s),lwd=2,col="red")
mtext("(f)",side=2,line=1.2,at=0.2375,cex=1,bg="white",las=2)

dev.off()


##############################
# 3-6 Figure S10 in Supp Mat #                             
##############################

pdf("outputs/Figure_S10_june2023.pdf",width=7.2,height=7.5)
#d1<-layout(mat=rbind(c(1,1,1,1,1,3,3,3,3,3),c(1,1,1,1,1,3,3,3,3,3),c(1,1,1,1,1,3,3,3,3,3), c(2,2,2,2,2,4,4,4,4,4),
#                     c(5,5,5,5,5,7,7,7,7,7),c(5,5,5,5,5,7,7,7,7,7),c(5,5,5,5,5,7,7,7,7,7), c(6,6,6,6,6,8,8,8,8,8)))

layout(mat=rbind(c(1,1,1,1,1,3,3,3,3,3),c(1,1,1,1,1,3,3,3,3,3),c(1,1,1,1,1,3,3,3,3,3),c(2,2,2,2,2,4,4,4,4,4),
                    c(5,5,5,5,5,7,7,7,7,7),c(5,5,5,5,5,7,7,7,7,7),c(5,5,5,5,5,7,7,7,7,7),c(6,6,6,6,6,8,8,8,8,8),
                    c(9,9,9,9,9,11,11,11,11,11),c(9,9,9,9,9,11,11,11,11,11),c(9,9,9,9,9,11,11,11,11,11),c(10,10,10,10,10,12,12,12,12,12)))

ylim=c(-53,73) ; xlim=c(-167,167)

### Plot SR (A)

par(mar=c(0,2,2,0.5)); 

Map_the_world_3 (Data_PA=Data_end_shark$SR,coord_X=Data_end_shark$x,coord_Y=Data_end_shark$y,
                 col=col1,breaks= seq(0,150,10),include.lowest=TRUE,xlim=xlim,ylim=ylim,coast=coast,names_fig="A",
                 cex_names_fig=0.7,x_names=0.5,y_names=87,cex.axis.leg=0.4,cex_point=0.25,
                 bg_col="gray80",bathy=FALSE,legend=FALSE,ncol=2,col_coast="gray20")
par(mar=c(0,2,0.2,0.5))
get_leg(data=Data_end_shark$SR,breaks=seq(0,150,10),cex=0.6,Colour=col1,title="Species Richness (SR))",
        ncol=5,bty="n",pt.cex=1)


####### Plot Residuals of the linear relationship between FRic and PD (C)

par(mar=c(0,2,2,0.5))

Map_the_world_3 (Data_PA=FD_PDres,coord_X=Coord[,1],coord_Y=Coord[,2],
                 col=c("#173869","#4575B4","#87a3cc", "#fcfafa", "#EC9E9A", "#a1312b", "#821611"),breaks=c(seq(-0.4,0,0.1),seq(0.1,0.3,0.1)),include.lowest=T,xlim=xlim,ylim=ylim,coast=coast,names_fig="C",
                 cex_names_fig=0.7,x_names=0.5,y_names=87,cex.axis.leg=0.4,cex_point=0.25,legend=F,bg_col="gray80",col_coast="gray20")
par(mar=c(0,2,0.1,0.5))
get_leg(data=FD_PDres,breaks=c(seq(-0.4,0,0.1),seq(0.1,0.3,0.1)),cex=0.6,Colour=c("#173869","#4575B4","#87a3cc", "#fcfafa", "#EC9E9A", "#a1312b", "#821611"),
        title="Residuals of the linear relationship between FRic and PD",
        ncol=5,bty="n",pt.cex=1)


#### Plot PD: Phylogenetic Diversity (B)

par(mar=c(0.5,2,1.5,0.5))
Map_the_world_3 (Data_PA=Data_end_shark$PD,coord_X=Coord[,1],coord_Y=Coord[,2],
                 col=col1,breaks=c(seq(400,12000,600)),include.lowest=T,xlim=xlim,ylim=ylim,coast=coast,names_fig="B",
                 cex_names_fig=0.7,x_names=0.5,y_names=87,cex.axis.leg=0.4,cex_point=0.25,legend=F,bg_col="gray80",col_coast="gray20")
par(mar=c(0,2,0.1,0.5))
get_leg(data=Data_end_shark$PD,breaks=c(seq(400,12000,600)),cex=0.6,Colour=col1,
        title="Phylogenetic Diversity (PD)", ncol=5,bty="n",pt.cex=1)

####### Plot Residuals of the linear relationship between PD and SR (D)

par(mar=c(0.5,2,1.5,0.5))

Map_the_world_3 (Data_PA=DataPD$residuals_PD,coord_X=Coord[,1],coord_Y=Coord[,2],
                 col=c("#173869","#4575B4","#87a3cc", "#fcfafa", "#EC9E9A", "#821611"),breaks=c(seq(-2500,0,500),seq(500, 1000, 500)),include.lowest=T,xlim=xlim,ylim=ylim,coast=coast,names_fig="D",
                 cex_names_fig=0.7,x_names=0.5,y_names=87,cex.axis.leg=0.4,cex_point=0.25,legend=F,bg_col="gray80",col_coast="gray20")
par(mar=c(0,2,0.1,0.5))
get_leg(data=DataPD$residuals_PD,breaks=c(seq(-2500,0,500),seq(500,1000,500)),cex=0.6,Colour=c("#173869","#4575B4","#87a3cc", "#fcfafa", "#EC9E9A", "#821611"),
        title="Residuals of the linear relationship between PD and SR",
        ncol=5,bty="n",pt.cex=1)

### Plot FUn considering the global pool of species (E)

par(mar=c(0,2,2,0.5)); 

Map_the_world_3(Data_PA=DataFD$mean.FUn_std.grids,coord_X=Data_end_shark$x,coord_Y=Data_end_shark$y,
                col=col1,breaks= seq(0,0.12,0.01),include.lowest=TRUE,xlim=xlim,ylim=ylim,coast=coast,names_fig="E",
                cex_names_fig=0.7,x_names=0.5,y_names=87,cex.axis.leg=0.4,cex_point=0.25,
                bg_col="gray80",bathy=FALSE,legend=FALSE,ncol=2,col_coast="gray20")
par(mar=c(0,2,0.2,0.5))
get_leg(data=DataFD$mean.FUn_std.grids,breaks=seq(0,0.12,0.01),cex=0.6,Colour=col1,title="FUn",
        ncol=5,bty="n",pt.cex=1)

dev.off()


##############################
# 3-7 Figure S11 in Supp Mat #                             
##############################                 

pdf("outputs/Figure_S11_june2023.pdf",width=7.2,height=7.5)
d<-layout(mat=rbind(c(1,1,1,1,1,3,3,3,3,3),c(1,1,1,1,1,3,3,3,3,3),c(1,1,1,1,1,3,3,3,3,3),c(2,2,2,2,2,4,4,4,4,4),
                    c(5,5,5,5,5,7,7,7,7,7),c(5,5,5,5,5,7,7,7,7,7),c(5,5,5,5,5,7,7,7,7,7),c(6,6,6,6,6,8,8,8,8,8),
                    c(9,9,9,9,9,11,11,11,11,11),c(9,9,9,9,9,11,11,11,11,11),c(9,9,9,9,9,11,11,11,11,11),c(10,10,10,10,10,12,12,12,12,12)))

ylim=c(-53,73); xlim=c(-167,167)

### Number of Threatened species (A)

par(mar=c(0,2,2,0.5)); 

Map_the_world_3 (Data_PA=Data_end_shark$SR_treat,coord_X=Data_end_shark$x,coord_Y=Data_end_shark$y,
                 col=col1,breaks= seq(0,70,5),include.lowest=TRUE,xlim=xlim,ylim=ylim,coast=coast,names_fig="A",
                 cex_names_fig=0.7,x_names=0.5,y_names=87,cex.axis.leg=0.4,cex_point=0.25,
                 bg_col="gray80",bathy=FALSE,legend=FALSE,ncol=2,col_coast="gray20")
par(mar=c(0,2,0.2,0.5))
get_leg(data=Data_end_shark$FRic,breaks=seq(0,70,5),cex=0.6,Colour=col1,title="Number of threatened species",
        ncol=5,bty="n",pt.cex=1)

### Number of top 25% ED2 species (C)

par(mar=c(0,2,2,0.5)); 

Map_the_world_3 (Data_PA=Data_end_shark$SR_Top_ED2_25,coord_X=Data_end_shark$x,coord_Y=Data_end_shark$y,
                 col=col1,breaks= seq(0,60,5),include.lowest=TRUE,xlim=xlim,ylim=ylim,coast=coast,names_fig="D",
                 cex_names_fig=0.7,x_names=0.5,y_names=87,cex.axis.leg=0.4,cex_point=0.25,
                 bg_col="gray80",bathy=FALSE,legend=FALSE,ncol=2,col_coast="gray20")
par(mar=c(0,2,0.2,0.5))
get_leg(data=Data_end_shark$SR_Top_ED2_25,breaks=seq(0,60,5),cex=0.6,Colour=col1,title="Number of top 25% ED2 species",
        ncol=5,bty="n",pt.cex=1)

### Number of top 25% FUn species (B)

par(mar=c(0,2,2,0.5)); 

Map_the_world_3 (Data_PA=Data_end_shark$SR_Top_FUN25,coord_X=Data_end_shark$x,coord_Y=Data_end_shark$y,
                 col=col1,breaks= seq(0,50,5),include.lowest=TRUE,xlim=xlim,ylim=ylim,coast=coast,names_fig="B",
                 cex_names_fig=0.7,x_names=0.5,y_names=87,cex.axis.leg=0.4,cex_point=0.25,
                 bg_col="gray80",bathy=FALSE,legend=FALSE,ncol=2,col_coast="gray20")
par(mar=c(0,2,0.2,0.5))
get_leg(data=Data_end_shark$SR_Top_FUN25,breaks=seq(0,50,5),cex=0.6,Colour=col1,title="Number of top 25% FUn species",
        ncol=5,bty="n",pt.cex=1)

#### Number of top 25% EDGE2 species (E)

par(mar=c(0,2,2,0.5));
Map_the_world_3 (Data_PA=Data_end_shark$SR_Top_EDGE2,coord_X=Coord[,1],coord_Y=Coord[,2],
                 col=col1,breaks=c(seq(0,50,5)),include.lowest=T,xlim=xlim,ylim=ylim,coast=coast,names_fig="E",
                 cex_names_fig=0.7,x_names=0.5,y_names=87,cex.axis.leg=0.4,cex_point=0.25,legend=F,bg_col="gray80",col_coast="gray20")
par(mar=c(0,2,0.1,0.5))
get_leg(data=Data_end_shark$SR_Top_EDGE2,breaks=c(seq(0,50,5)),cex=0.6,Colour=col1,
        title="Number of top 25% EDGE2 species", ncol=5,bty="n",pt.cex=1)

#### Number of top 25% FSp species (C)

par(mar=c(0,2,2,0.5));
Map_the_world_3 (Data_PA=Data_end_shark$SR_Top_FSp25,coord_X=Coord[,1],coord_Y=Coord[,2],
                 col=col1,breaks=c(seq(0,50,5)),include.lowest=T,xlim=xlim,ylim=ylim,coast=coast,names_fig="C",
                 cex_names_fig=0.7,x_names=0.5,y_names=87,cex.axis.leg=0.4,cex_point=0.25,legend=F,bg_col="gray80",col_coast="gray20")
par(mar=c(0,2,0.1,0.5))
get_leg(data=Data_end_shark$SR_Top_FSp25,breaks=c(seq(0,50,5)),cex=0.6,Colour=col1,
        title="Number of top 25% FSp species", ncol=5,bty="n",pt.cex=1)

#### Number of top 25% FUSE species (F)

par(mar=c(0,2,2,0.5));
Map_the_world_3 (Data_PA=Data_end_shark$SR_Top_FUSE100,coord_X=Coord[,1],coord_Y=Coord[,2],
                 col=col1,breaks=c(seq(0,45,5)),include.lowest=T,xlim=xlim,ylim=ylim,coast=coast,names_fig="F",
                 cex_names_fig=0.7,x_names=0.5,y_names=87,cex.axis.leg=0.4,cex_point=0.25,legend=F,bg_col="gray80",col_coast="gray20")
par(mar=c(0,2,0.1,0.5))
get_leg(data=Data_end_shark$SR_Top_FUSE100,breaks=c(seq(0,45,5)),cex=0.6,Colour=col1,
        title="Number of top 25% FUSE species", ncol=5,bty="n",pt.cex=1)

dev.off()


##############################
# 3-8 Figure S12 in Supp Mat #
##############################

pdf("outputs/Figure_S12_june2023.pdf",width=8,height=6)
layout(mat=rbind(c(1,1,1,1,1,2,2,2,2,2),c(1,1,1,1,1,2,2,2,2,2),c(1,1,1,1,1,2,2,2,2,2),c(1,1,1,1,1,2,2,2,2,2),
                 c(3,3,3,3,3,4,4,4,4,4),c(3,3,3,3,3,4,4,4,4,4),c(3,3,3,3,3,4,4,4,4,4),c(3,3,3,3,3,4,4,4,4,4)))

# Hotspots of FSp (2.5% greatest value - hotspots) - A

par(mar=c(3,2,1.5,0.5))
carto_congruence_single (nbvar=1,data="Hotspot",x="FSp",
                         names_fig="A",col_pt=c("#B70F0D"),cex_point=0.2,cex_text_fig=0.75,
                         cex_legend=0.6,main_legend="FRic hotspots",ncol=1,ylim=c(-53,73), xlim=c(-167,167),
                         cex.axis=0.4,axisx=T,axisy=T,legend=FALSE)

mtext("A",side=2,line=0.5,at=87,cex=0.75,bg="white",las=2)

# Congruence between FRic and PD - B

par(mar=c(3,2,1.5,1))
carto_function(y=Hotspot,x="PD", index2="FRIC",names_fig="B",cex_point=0.2,cex_text_fig=0.75,cex_legend=0.6,
               text=c("Non-hotspot values","Values in top 2.5% for FUSE","Values in top 2.5% for EDGE","Congruence zone"),
               xlim=xlim,ylim=ylim,pos_legY=87,pos_legX=0.5,cex.axis=0.4,ncol=1,pos_leg="bottomleft",pt.cexleg=1,legend="")

# Congruence between SR, PD, FRic and fishing impact - C

par(mar=c(3,2,1.5,1))

carto_congruence_t(nbvar=4, data="Hotspot",data2="Hotspot5" ,data3="Hotspot10",x="RS",index2="FRIC",index3= "PD", index4="CummulFishing_impact_2012_2020",
                   names_fig="",ylim=c(-53,73), xlim=c(-167,167),
                   cex_point=0.2,cex_text_fig=0.45,cex_legend=0.45,
                   text=c(""), col_pt=c("red","#fdae61","purple"),pos_legY=86.5,pos_legX=-14.75,cex.axis=0.4,ncol=1, legend=FALSE)

mtext("C",side=2,line=0.5,at=87,cex=0.75,bg="white",las=2)

# Congruence between EDGE2,  FUSE and fishing impact - D

carto_congruence_t(nvar=3, data="Hotspot",data2="Hotspot5" ,data3="Hotspot10",x="SR_Top_EDGE2",index2="SR_Top_FUSE100",index3="CummulFishing_impact_2012_2020",
                   names_fig="",ylim=c(-53,73), xlim=c(-167,167),
                   cex_point=0.2,cex_text_fig=0.45,cex_legend=0.45,
                   text=c(""),
                   col_pt=c("red","#fdae61","purple"),pos_legY=86.5,pos_legX=-14.75,cex.axis=0.4,ncol=1, legend=FALSE)

mtext("D",side=2,line=0.5,at=87,cex=0.75,bg="white",las=2)


#### add the legend at the end
par(new=T, fig=c(0,1,0,1),mar=c(0,2,0.1,0.5))
plot(NULL,xaxt='n',yaxt='n',bty='n',ylab='',xlab='',xlim=0:1,ylim=0:1)

#A
legend(x=0 ,y=0.55,legend=c("No values","Non hotspots values","2.5% of FSp values"),
       pch=19,bty="n",col=c("#B6E0EE","#6BAACC","#B70F0D"),
       cex=0.72,horiz=T)

#B
legend(x=0.54,y=0.55,legend=c( text=c("No values","Non-hotspot values","Values in top 2.5% for FRic",
                                      "Values in top 2.5% for PD","Congruence zone")),
       pch=19,bty="n",col=c("#B6E0EE","#6BAACC","#B70F0D","#EDDF19","#fdae61"),
       cex=0.72,ncol=3)

#C
legend(x=0,y=0.01,legend=c("No values","Non hotspots values","2.5% of the four indices", "5% of the four indices", "10% of the four indices"),
       pch=19,bty="n",col=c("#B6E0EE","#6BAACC","red2", "#fdae61","purple"),
       cex=0.72,horiz=F, ncol=3)

#D
legend(x=0.54,y=0.01,legend=c("No values","Non hotspots values","2.5% of the three indices", "5% of the three indices", "10% of the three indices"),
       pch=19,bty="n",col=c("#B6E0EE","#6BAACC","red2", "#fdae61","purple"),
       cex=0.72,horiz=F, ncol=3)
dev.off()

#####################################################
#  3-9 Figure S13 in Supp Mat Draftman plot       #
#####################################################

shark_data<-data.frame(SR=Data_end_shark$SR, PD=Data_end_shark$PD, FRic=Data_end_shark$FRic,
                       FUn=Data_end_shark$FUn, PUn=Data_end_shark$PUn,
                       FSp=Data_end_shark$FSp, SR_treat=Data_end_shark$SR_treat, SR_Top25_EDGE=Data_end_shark$SR_Top_EDGE2,
                       SR_Top25_FUSE=Data_end_shark$SR_Top_FUSE100, SR_Top25_ED2=Data_end_shark$ED2,
                       SR_Top25_FUn=Data_end_shark$SR_Top_FUN25, SR_Top25_FSp=Data_end_shark$SR_Top_FSp25)


pdf("outputs/Figure_S13_june2023.pdf",width=7,height=7)
corrplot(cor(shark_data, use="na.or.complete"), method = 'ellipse', order = 'AOE', type = 'upper')
dev.off()


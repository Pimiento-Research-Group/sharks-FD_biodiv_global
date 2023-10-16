################################################################################
#======================== Functions Script Carto shark PD =====================#
#==================== 09 avril 2020 albouycamille@gmail.com ===================#
################################################################################

Map_the_world_3 <- function(Data_PA=grid$P_A,coord_X=grid$X_ENTIER,coord_Y=grid$Y_ENTIER,
                            col=c("blue", "red3"),breaks= seq(0,40,2),include.lowest=T,
                            xlim=c(-180,180),ylim=c(-90,90),coast=coast,names_fig="a)",
                            cex_names_fig=1,x_names=0.5,y_names=99,cex.axis.leg=0.5,cex_point=0.6,
                            bg_col="gray80",bathy=F,legend=T,digit=2,title_leg="lala",
                            adShape=FALSE,Poly_sup,ncol=2,col_coast="gray20",colorlegend=F,
                            leg_incr=2,inset_leg=c(-0.21,0.2),axisy=T,axisx=T,...){
  
  require(shape)
  plot(rnorm(1000),col="white",type="n",axes=F,xlab="",ylab="",xlim=xlim,ylim=ylim, bty="L")
  if(nchar(bg_col)!=0){rect(xleft=-181.5,ybottom=-90,xright=181.5,ytop=90,col=bg_col)}
  
  if(bathy==T){
    colour_rast <- colorRampPalette(c("gray47","gray80"))(50)
    plot(r_bathy,add=T,col=colour_rast,legend=F,xlim=c(-190,190),ylim=c(-90,90))
  }
  
  colour <- colorRampPalette(col)(length(breaks))
  richesse <- cut(Data_PA,breaks=breaks,include.lowest=include.lowest,dig.lab=digit,right=F)
  vect_col <- colour[richesse]
  vect_col[which(is.na(vect_col))] <- "gray80"
  
  points(coord_X,coord_Y,col=vect_col,pch=15,cex=cex_point) 
  
  if(is.null(coast)==F){plot(coast,col=col_coast,add=T,border="gray25",lwd=0.1)} # end of if
  if(adShape==T){plot(Poly_sup,add=T,border="#D70230",density=25,lwd=0.6,col="grey30")}
  if(legend==T){
    legend("topright",legend=levels(richesse),col=colour,ncol=ncol,pch=16,pt.cex=1.2,cex=0.5,lwd=0.2,lty=0,
           title=title_leg,text.col="black",bg= "white",inset=inset_leg,xpd=T,bty="n")
  }
  abline(h=0,lty="dotted",col="grey30")
  if(axisx==T){
    axis(side=1,line=0,cex.axis=0.5,lwd=0.35,tcl=-0.25,bg="white",labels=c("-180°","-135°","-90°","-45°","0°","45°","90°","135°","180°"),
         at=c(-179,-135,-90,-45,0,45,90,135,179),mgp=c(3,0.05,0))
  } else {
    axis(side=1,line=0,cex.axis=0.5,lwd=0.35,tcl=-0.25,bg="white",labels=c("","","","","","","","",""),
         at=c(-179,-135,-90,-45,0,45,90,135,179),mgp=c(3,0.05,0))
  }
  
  if(axisy==T){
    axis(side=2,line=0,las=2,cex.axis=0.5,lwd=0.35,tcl=-0.25,bg="white",labels=c("-75°","-50°","-25°","0°","25°","50°","75°"),
         at=c(-75,-50,-25,0,25,50,75),mgp=c(3,0.25,0))
  } else{
    axis(side=2,line=0,las=2,cex.axis=0.5,lwd=0.35,tcl=-0.25,bg="white",labels=c("","","","","","",""),
         at=c(-75,-50,-25,0,25,50,75),mgp=c(3,0.25,0))
  }
  
  mtext(names_fig,side=2,line=x_names,at=y_names,cex=cex_names_fig,bg="white",las=2)    
  
  if (colorlegend==T) {
    colorlegend(zlim=round(c(min(breaks),max(breaks)),2),zval=breaks,col=colour,main=title_leg,xpd=FALSE,main.cex=0.6,cex=0.7,posx=c(0.87, 0.9),...)
  }
  box(lwd=0.6)
} # Map_the_world_3

################################################################################
get_leg <- function (data=my_coast$FRic,breaks=c(seq(0,0.91, 0.05), 0.9189),
                     Colour=Colour,title="Functional Richness",...) {
  emptyplot()
  legend <- levels(cut(data,breaks=breaks,include.lowest=TRUE,dig.lab=4,right=F))
  jet.color <- colorRampPalette(Colour)
  colour <- jet.color(length(breaks))
  legend("center",pch=15,legend=legend,col=colour,title=title,text.col="black",...)
} # end of get_leg

##############################################################################
get_lat_grad <- function(Data,ylim,xlim,Label,Labelpos=0.3,col.poly="#42BAB1",
                         xlab=c("0","0.2","0.4","0.6"),xat=c(0.005,0.2,0.4,0.6),cex.axis=0.3){

  Lat_FR <- do.call(rbind,lapply(split(Data,f=Data[,2]),function(x){mean(x[,3],na.rm=T)}))
  Lat_FR <- data.frame(Latitude=as.numeric(rownames(Lat_FR)),as.numeric(Lat_FR),stringsAsFactors=F)
  Lat_FR <- na.omit(Lat_FR)
  
  plot(y=Lat_FR[,1],x=Lat_FR[,2],xlab="",ylab="",axes=F,type="n",col="grey10",
         cex.lab=1,mgp=c(2.1,0.8,0),xlim=xlim,ylim=ylim,pch=19,cex=0.4,xaxs="i")
  mtext(Label,side=1,line=2,at=Labelpos,cex=0.45)
    
  axis(side=2,line=0,las=2,cex.axis=cex.axis,lwd=0.35,tcl=-0.25,bg="white",mgp=c(3,0.3,0),
         labels=c("","","","","","",""),at=c(-75,-50,-25,0,25,50,75))
  axis(side=1,line=0,lwd=0.35,bg="white",labels=xlab,at=xat,mgp=c(3,0.05,0),
         cex.axis=cex.axis,tcl=-0.25)
    
  lw1 <- predict(loess(Lat_FR[,2]~as.numeric(Lat_FR[,1]),span=0.1,degree=2),se=T)
  polygon(y=c(Lat_FR[,1],rev(as.numeric(Lat_FR[,1]))),x=c(lw1$se.fit,rev(lw1$fit)),col=col.poly,border=NA)
    
  lines(lw1$fit,as.numeric(Lat_FR[,1]),col="black",lwd=0.3)
  dt <- data.frame(as.numeric(Lat_FR[,1]),lw1$fit)
  points(x=dt[which.max(dt[,2]),2],y=dt[which.max(dt[,2]),1],pch=21,bg="orange",col="black",cex=0.75,lwd=0.3)
  abline(h=0,lty=2,cex=0.2,lwd=0.3); box(lwd=0.3) 
  
} # get get_lat_grad 

###############################################################################


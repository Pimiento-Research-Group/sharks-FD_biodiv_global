 ######################################
 #
 # World mapping fonction
 #
#
######################################

Map_the_world <- function (Data_PA=grid$P_A, coord_X=grid$X_ENTIER,coord_Y=grid$Y_ENTIER,
                           col=c("White","yellow","red3"),breaks=seq(0,40,2),include.lowest=T,
                           xlim=c(-166,166),ylim=c(-83.5,83.5),zlevels=6,main_legend_cex=1.2,
                           coastline=T,main_legend="Species Richness",output=F,legend_increment=100,
                           names_fig="a)",cex_names_fig=1,...){
    require(shape)
     
    if (output==TRUE) {
      
         tiff(filename="world_map.tiff",width=6.5,height=4.5,units="in",res=300,compression="lzw",type="cairo")
         plot(rnorm(1000),col="white",type="n",axes=F,xlab="",ylab="",xlim=xlim,ylim=ylim)
         if (coastline==TRUE) {plot(coast,col="black",add=TRUE)} # end of if
        
         jet.color <-  colorRampPalette(col)
         colour <-  jet.color(length(breaks)-1)
         richesse <-  cut(Data_PA , breaks=breaks,include.lowest=include.lowest)
         
         points(coord_X,coord_Y,col=colour[richesse],pch=15,cex=0.3) 
         axis(side=1,line=0,cex.axis=0.7,lwd=0.35,tcl=-0.25,bg="white",labels=c("-180","-135","-90","-45","0","45","90","135","180"),at=c(-179,-135,-90,-45,0,45,90,135,179),mgp=c(3,0.5,0))
         axis(side=2,line=0,las=2,cex.axis=0.7,lwd=0.35,tcl=-0.25,bg="white",labels=c("-90","-45","0","45","90"),at=c(-90,-45,0,45,90),mgp=c(3,0.5,0))
         mtext(names_fig,side=2,line=0.5,at= 103,cex=cex_names_fig,bg="white",las=2)  	
         box()
         
         colorlegend(zlim=round(c(min(Data_PA),max(Data_PA))), dz=legend_increment,col=colour, zlevels=zlevels,main=main_legend,xpd=FALSE,main.cex=0.6,cex=0.7,posx=c(0.87, 0.9))
         dev.off()
         
    } else {
         plot(rnorm(1000),col="white",type="n",axes=F,xlab="",ylab="",xlim=xlim,ylim=ylim)
         jet.color <- colorRampPalette(col)
         colour <- jet.color(length(breaks)-1)
         
         richesse <- cut(Data_PA,breaks=breaks,include.lowest=include.lowest)
         points(coord_X,coord_Y,col=colour[richesse],pch=15,cex=0.3) 
         
       if (coastline==TRUE) {plot(coast,col="black",add=TRUE,lwd=0.02)} # end of if
          
         axis(side=1,line=0,cex.axis=0.8,lwd=0.35,tcl=-0.25,bg="white",labels=c("-180","-135","-90","-45","0","45","90","135","180"),at=c(-179,-135,-90,-45,0,45,90,135,179),mgp=c(3,0.5,0))
         axis(side=2,line=0,las=2,cex.axis=0.8,lwd=0.35,tcl=-0.25,bg="white",labels=c("-90","-45","0","45","90"),at=c(-90,-45,0,45,90),mgp=c(3,0.5,0))
         mtext(names_fig,side=2,line=0.3,at=87,cex=cex_names_fig,bg="white",las=2)  		
         box()
         colorlegend(zlim=round(c(min(Data_PA),max(Data_PA))),dz=legend_increment,col=colour,zlevels=zlevels,main=main_legend,xpd=FALSE,main.cex=main_legend_cex,cex=0.7,posx=c(0.92, 0.94))
    } # end of ifelse
  } # end of function Mpa_the_world

#################################################
# A modifier pour rendre plus fléxible
#  Map the hotspot and congruence between hotspot
##################################################
col = rainbow(5)

carto_congruence <- function(data="Hotspot",data2="Hotspot5",data3="Hotspot10",x="RS_tot",index2="FRiC_r",
                             index3="PD",names_fig="a",xlim=c(-166,166),ylim=c(-83.5,83.5),cex_point=0.3,
                             cex_text_fig=1,cex_legend=0.5,text=c(" 10% of the three variables","20% of the three variables",
                             "30% of the three variables"),cex_names_fig=0.1,
                             main_legend="Congruence zone between PD,FD, and SR with values in the upper"){
    
    col <- c("red","#CC00FFFF","lightskyblue3")
    w <- eval(parse(text=paste("cbind(",data,"[which(",data,"$",x,"==1 & ",data,"$",index2,"==1 & ",data,"$",index3,"==1 ),],color=col[1])",sep="")))# hot spot
    w5 <- eval(parse(text=paste("cbind(",data2,"[which(",data2,"$",x,"==1 & ",data2,"$",index2,"==1 & ",data2,"$",index3,"==1 ),],color=col[2])",sep="")))# hot spot

    plot(rnorm(1000),col="white",bg="grey",axes=F,xlab="",ylab="",xlim=xlim,ylim=ylim)
    rect(xleft=-180,ybottom =-90,xright=180,ytop=90,density=NULL,angle=45,col="#EEDDB6",border=NULL)
    
    points(x=w5[,"X"],y=w5[,"Y"],col=as.character(w5$color),pch=15,cex=0.3)
    points(x=w[,"X"],y=w[,"Y"],col=as.character(w$color),pch=15,cex=0.3)
    plot(coast,col="black",add=TRUE)
    
    axis(side=1,line=0,cex.axis=0.7,lwd=0.35,tcl=-0.25,bg="white",labels=c("-180","-135","-90","-45","0","45","90","135","180"),at=c(-179,-135,-90,-45,0,45,90,135,179),mgp=c(3, 0.5, 0))
    axis(side=2,line=0,las=2,cex.axis=0.7,lwd=0.35,tcl=-0.25,bg="white",labels=c("-90","-45","0","45","90"),at=c(-90,-45,0,45,90),mgp=c(3, 0.5, 0))
    mtext(names_fig,side=2,line=0.5,at= 103,cex=cex_names_fig,bg="white",las=2)  	
    
    box()
    legend("topright",leg=text,pch=15,col=col,bg="white",cex=cex_legend,lwd=0.2,box.lwd=0.6,lty=0,title =main_legend,ncol=2)

}# end of function carto_congruence

###################################################################################
carto_congruence_2 <- function(data="Hotspot",data2="Hotspot5",x="RS_tot",index2="FRiC_r",
                               index3="PD",names_fig="a",xlim=c(-166,166),ylim=c(-83.5,83.5),
                               cex_point=0.3,cex_text_fig=1,cex_legend=0.5,text=c(" 10% of the three variables","20% of the three variables"),
                               cex_names_fig=0.1,main_legend="Congruence zones between PD,FD, and SR with values in the upper",
                               ncol=2,bathy=F){
    
    col <- c("red","#993399")
    w <- eval(parse(text=paste("cbind(",data,"[which(",data,"$",x,"==1 & ",data,"$",index2,"==1 & ",data,"$",index3,"==1 ),],color=col[1])",sep="")))# hot spot
    w5 <- eval(parse(text=paste("cbind(",data2,"[which(",data2,"$",x,"==1 & ",data2,"$",index2,"==1 & ",data2,"$",index3,"==1 ),],color=col[2])",sep="")))# hot spot
   
    plot(rnorm(1000),col="white",bg="grey",axes=F,xlab="",ylab="",xlim=xlim,ylim=ylim)
    
    if(bathy==TRUE){plot(r_bathy,col=colour_rast,add=TRUE,legend=FALSE)}
    else{rect(xleft =-180,ybottom=-90,xright=180,ytop=90,density=NULL,angle=45,col="#EEDDB6",border=NULL)}
    
    points(x=w5[,"X"],y=w5[,"Y"],col=as.character(w5$color),pch=15,cex=0.3)
    points(x=w[,"X"],y=w[,"Y"],col=as.character(w$color),pch=15,cex=0.3)
    abline(h=0,lty="dotted",col="grey30")
    
    plot(coast,col="black",add=TRUE,lwd=0.02)
    axis(side=1,line=0,cex.axis=0.7,lwd=0.35,tcl=-0.25,bg="white",labels=c("-180","-135","-90","-45","0","45","90","135","180"),at=c(-179,-135,-90,-45,0,45,90,135,179),mgp=c(3, 0.15,0))
    axis(side=2,line=0,las=2,cex.axis=0.7,lwd=0.35,tcl=-0.25,bg="white",labels=c("-90","-45","0","45","90"),at=c(-90,-45,0,45,90),mgp=c(3,0.5,0))
    mtext(names_fig,side=2,line=0.5,at= 103,cex=cex_names_fig,bg="white",las=2)  	
    
    rect(50,-80,180,-55,col="#FFFFFF99")
    legend(x=55,y=-55,leg=text,pch=15,col=col,bg="white",pt.cex=0.9,cex=cex_legend,lwd=0.2,box.lwd=0.6,lty=0,title=main_legend,ncol=1,bty="n",x.intersp=0.5)
    
    box()
}# end of function carto_congruence

###################################################################################
get_mini_plot <- function(xlim=c(129,158),ylim=c(30,48),pos=c(0.71,0.98,0.57,0.93),col=c("red","#993399"),bathy=r_bathy){
  
  par(fig=pos,new=TRUE)
  plot(rnorm(1000),col="white",type="n",axes=F,xlab="",ylab="",xlim=xlim,ylim=ylim)
  plot(bathy,col=colour_rast,legend=FALSE,add=TRUE)
  
  par(fig=pos,new=TRUE)
  
  data="Hotspot";data2="Hotspot5";x="RS_tot";index2="FRiC_r";index3 = "PD_rel"
  w <- eval(parse(text=paste("cbind(",data,"[which(",data,"$",x,"==1 & ",data,"$",index2,"==1 & ",data,"$",index3,"==1 ),],color=col[1])",sep="")))# hot spot
  w5 <- eval(parse(text=paste("cbind(",data2,"[which(",data2,"$",x,"==1 & ",data2,"$",index2,"==1 & ",data2,"$",index3,"==1 ),],color=col[2])",sep="")))# hot spot
  
  plot(rnorm(1000),col="white",type="n",axes=F,xlab="",ylab="",xlim=xlim,ylim=ylim)
  points(x=w5[,"X"],y=w5[,"Y"],col=as.character(w5$color),pch=15,cex=0.6)
  points(x=w[,"X"],y=w[,"Y"],col=as.character(w$color),pch=15,cex=0.6)
  
  plot(coast,col="black",add=TRUE,lwd=0.3)
  box(lwd=0.8,col="white")
  
} # get_mini_plot

####################################################################################
carto_congruence_single <- function(nbvar=3,data="Hotspot",x="RS_tot",index2="FRiC_r",
                      names_fig="a",xlim=c(-166,166),ylim=c(-50,55),col_pt=c("#fdae61","#EDDF19","#B70F0D"),cex_point=0.3,cex_text_fig=1,
                      cex_legend=0.5,text=c("Non hotspots values"," 2.5% of the three variables"),
                      main_legend="",ncol=1,legend=T,
                      pos_leg="bottomleft",pt.cexleg=1,cex.axis=0.4,pos_legX=88,pos_legY=0.5,axisx=T,axisy=T,...){

    col <- col_pt
  
    if(nbvar==1){
      w <- eval(parse(text=paste("cbind(",data,"[which(",data,"$",x,"==1 ),],color=col[1])",sep="")))# hot spot
      Nothing <- eval(parse(text=paste("cbind(",data,"[which(",data,"$",x,"==0),],color='#6BAACC')",sep="")))# hot spot
    }

    plot(rnorm(1000),col="white",bg="grey",axes=F,xlab="",ylab="",xlim=xlim,ylim=ylim)
    rect(xleft =-180,ybottom=-90,xright=180,ytop=90,density=NULL,angle=45,col="#B6E0EE",border=NULL)
    points(x= Nothing[,"X"],y=Nothing[,"Y"],col=as.character(Nothing$color),pch=15,cex=cex_point)
    points(x=w[,"X"],y=w[,"Y"],col=as.character(w$color),pch=15,cex=cex_point)
    abline(h=0,lty="dotted",col="grey30")

    plot(coast,col="black",add=TRUE,lwd=0.02)
    if(axisx==T){
    axis(side=1,line=0,cex.axis=0.5,lwd=0.35,tcl=-0.25,bg="white",labels=c("-180°","-135°","-90°","-45°","0°","45°","90°","135°","180°"),
         at=c(-179,-135,-90,-45,0,45,90,135,179),mgp=c(3,0.05,0))
    }
    
    if(axisy==T){
      axis(side=2,line=0,las=2,cex.axis=0.5,lwd=0.35,tcl=-0.25,bg="white",labels=c("-75°","-50°","-25°","0°","25°","50°","75°"),
           at=c(-75,-50,-25,0,25,50,75),mgp=c(3,0.25,0))
      
    }

    mtext(names_fig,side=2,line=pos_legX,at=pos_legY,cex=cex_text_fig,bg="white",las=2)
        colour=c("#B6E0EE","#6BAACC",col_pt)
  if(legend==T){
    legend(pos_leg,leg=c("No values",text),col="black",pt.bg=colour,pch=22,cex=cex_legend,lty=NA,
           lwd=0.2,box.lwd=0.6,pt.cex=pt.cexleg,ncol=ncol,bg="#FFFFFF95",x.intersp=0.1,...)
  }      

    box(lwd=0.6)

}# end of function carto_congruence


###################################################################################
carto_congruence_t <- function(nbvar=3,data="Hotspot",data2="Hotspot5",data3="Hotspots10",x="RS_tot",index2="FRiC_r",index3="PD",index4="Threats",
                      names_fig="a",xlim=c(-166,166),ylim=c(-50,55),col_pt=c("#fdae61","#EDDF19","purple"),cex_point=0.3,cex_text_fig=1,
                      cex_legend=0.5,text=c("Non hotspots values"," 2.5% of the three variables","5% of the three variables"),
                      main_legend="Congruence zone between PD,FD, and SR with values in the upper",ncol=1,
                      pos_leg="bottomleft",pt.cexleg=1,cex.axis=0.4,pos_legX=88,pos_legY=0.5,axisx=T,axisy=T,legend=TRUE,...){

	col <- col_pt 
	if(nbvar==4){
	  	w <- eval(parse(text=paste("cbind(",data,"[which(",data,"$",x,"==1 & ",data,"$",index2,"==1 & ",data,"$",index3,"==1 &",
	  	          data,"$",index4,"==1 ),],color=col[1])",sep="")))# hot spot
	    w5 <- eval(parse(text=paste("cbind(",data2,"[which(",data2,"$",x,"==1 & ",data2,"$",index2,"==1 & ",data2,"$",index3,"==1 & ",
	               data2,"$",index4,"==1 ),],color=col[2])",sep="")))# hot spot
	    w10 <- eval(parse(text=paste("cbind(",data3,"[which(",data3,"$",x,"==1 & ",data3,"$",index2,"==1 & ",data3,"$",index3,"==1 & ",
	               data3,"$",index4,"==1 ),],color=col[3])",sep="")))# hot spot
	    Nothing <- eval(parse(text=paste("cbind(",data3,"[which(",data3,"$",x,"==0 | ",data3,"$",index2,"==0 | ",data3,"$",index3,"==0| ",data3,"$",index4,"==0),],color='#6BAACC')",sep="")))
	}
	
	if(nbvar==3){
	  w <- eval(parse(text=paste("cbind(",data,"[which(",data,"$",x,"==1 & ",data,"$",index2,"==1 & ",data,"$",index3,"==1 ),],color=col[1])",sep="")))# hot spot
	  w5 <- eval(parse(text=paste("cbind(",data2,"[which(",data2,"$",x,"==1 & ",data2,"$",index2,"==1 & ",data2,"$",index3,"==1),],color=col[2])",sep="")))# hot spot
	  w10 <- eval(parse(text=paste("cbind(",data3,"[which(",data3,"$",x,"==1 & ",data3,"$",index2,"==1 & ",data3,"$",index3,"==1),],color=col[3])",sep="")))# hot spot
	  Nothing <- eval(parse(text=paste("cbind(",data3,"[which(",data3,"$",x,"==0 | ",data3,"$",index2,"==0 | ",data3,"$",index3,"==0),],color='#6BAACC')",sep="")))# hot spot
	}
	
	if(nbvar==2){
	  w <- eval(parse(text=paste("cbind(",data,"[which(",data,"$",x,"==1 & ",data,"$",index2,"==1 ),],color=col[1])",sep="")))# hot spot
	  w5 <- eval(parse(text=paste("cbind(",data2,"[which(",data2,"$",x,"==1 & ",data2,"$",index2,"==1),],color=col[2])",sep="")))# hot spot
	  w10 <- eval(parse(text=paste("cbind(",data3,"[which(",data3,"$",x,"==1 & ",data3,"$",index2,"==1),],color=col[3])",sep="")))# hot spot
	  Nothing <- eval(parse(text=paste("cbind(",data3,"[which(",data3,"$",x,"==0 | ",data3,"$",index2,"==0),],color='#6BAACC')",sep="")))# hot spot
	}
	
	
	if(nbvar==1){
	  w <- eval(parse(text=paste("cbind(",data,"[which(",data,"$",x,"==1 ),],color=col[1])",sep="")))# hot spot
	  w5 <- eval(parse(text=paste("cbind(",data2,"[which(",data2,"$",x,"==1),],color=col[2])",sep="")))# hot spot
	  w10 <- eval(parse(text=paste("cbind(",data3,"[which(",data3,"$",x,"==1) ,],color=col[3])",sep="")))# hot spot
	  Nothing <- eval(parse(text=paste("cbind(",data3,"[which(",data3,"$",x,"==0),],color='#6BAACC')",sep="")))# hot spot
	}

	plot(rnorm(1000),col="white",bg="grey",axes=F,xlab="",ylab="",xlim=xlim,ylim=ylim)
	rect(xleft =-180,ybottom=-90,xright=180,ytop=90,density=NULL,angle=45,col="#B6E0EE",border=NULL)
	points(x= Nothing[,"X"],y=Nothing[,"Y"],col=as.character(Nothing$color),pch=15,cex=cex_point)
	points(x=w10[,"X"],y=w10[,"Y"],col=as.character(w10$color),pch=15,cex=cex_point)
	points(x=w5[,"X"],y=w5[,"Y"],col=as.character(w5$color),pch=15,cex=cex_point)
	points(x=w[,"X"],y=w[,"Y"],col=as.character(w$color),pch=15,cex=cex_point)
	abline(h=0,lty="dotted",col="grey30")

	plot(coast,col="black",add=TRUE,lwd=0.02)    
	if(axisx==T){
	axis(side=1,line=0,cex.axis=0.5,lwd=0.35,tcl=-0.25,bg="white",labels=c("-180°","-135°","-90°","-45°","0°","45°","90°","135°","180°"),
	     at=c(-179,-135,-90,-45,0,45,90,135,179),mgp=c(3,0.05,0))
	}
	
	if(axisy==T){
	  axis(side=2,line=0,las=2,cex.axis=0.5,lwd=0.35,tcl=-0.25,bg="white",labels=c("-75°","-50°","-25°","0°","25°","50°","75°"),
	       at=c(-75,-50,-25,0,25,50,75),mgp=c(3,0.25,0))
	  
	}

	mtext(names_fig,side=2,line=pos_legX,at=pos_legY,cex=cex_text_fig,bg="white",las=2)
		colour=c("#B6E0EE","#6BAACC",col_pt)
	if(legend==TRUE){
  legend(pos_leg,leg=c("No values",text),col="black",pt.bg=colour,pch=22,cex=cex_legend,lty=NA,
        lwd=0.2,box.lwd=0.6,pt.cex=pt.cexleg,ncol=ncol,bg="#FFFFFF95",x.intersp=0.1,...)
	}
	box(lwd=0.6)

}# end of function carto_congruence

###################################################################################
carto_congruence_t2<- function(nbvar=3,data="Hotspot",data2="Hotspot5",data3="Hotspots10",x="RS_tot",index2="FRiC_r",index3="PD",index4="Threats",
                               names_fig="a",xlim=c(-166,166),ylim=c(-50,55),col_pt=c("#fdae61","#EDDF19","purple"),cex_point=0.3,cex_text_fig=1,
                               cex_legend=0.5,text=c("Non hotspots values"," 2.5% of the three variables","5% of the three variables"),
                               main_legend="Congruence zone between PD,FD, and SR with values in the upper",ncol=1,
                               pos_leg="bottomleft",pt.cexleg=1,cex.axis=0.4,pos_legX=88,pos_legY=0.5,axisx=T,axisy=T,legend=TRUE,...){
  
  col <- col_pt 
  if(nbvar==4){
    w <- eval(parse(text=paste("cbind(",data,"[which(",data,"$",x,"==1 & ",data,"$",index2,"==1 & ",data,"$",index3,"==1 &",
                               data,"$",index4,"==1 ),],color=col[1])",sep="")))# hot spot
    w5 <- eval(parse(text=paste("cbind(",data2,"[which(",data2,"$",x,"==1 & ",data2,"$",index2,"==1 & ",data2,"$",index3,"==1 & ",
                                data2,"$",index4,"==1 ),],color=col[2])",sep="")))# hot spot
    w10 <- eval(parse(text=paste("cbind(",data3,"[which(",data3,"$",x,"==1 & ",data3,"$",index2,"==1 & ",data3,"$",index3,"==1 & ",
                                 data3,"$",index4,"==1 ),],color=col[3])",sep="")))# hot spot
    Nothing <- eval(parse(text=paste("cbind(",data3,"[which(",data3,"$",x,"==0 | ",data3,"$",index2,"==0 | ",data3,"$",index3,"==0| ",data3,"$",index4,"==0),],color='#6BAACC')",sep="")))
  }
  
  if(nbvar==3){
    w <- eval(parse(text=paste("cbind(",data,"[which(",data,"$",x,"==1 & ",data,"$",index2,"==1 & ",data,"$",index3,"==1 ),],color=col[1])",sep="")))# hot spot
    w5 <- eval(parse(text=paste("cbind(",data2,"[which(",data2,"$",x,"==1 & ",data2,"$",index2,"==1 & ",data2,"$",index3,"==1),],color=col[2])",sep="")))# hot spot
    w10 <- eval(parse(text=paste("cbind(",data3,"[which(",data3,"$",x,"==1 & ",data3,"$",index2,"==1 & ",data3,"$",index3,"==1),],color=col[3])",sep="")))# hot spot
    Nothing <- eval(parse(text=paste("cbind(",data3,"[which(",data3,"$",x,"==0 | ",data3,"$",index2,"==0 | ",data3,"$",index3,"==0),],color='#6BAACC')",sep="")))# hot spot
  }
  
  if(nbvar==2){
    w <- eval(parse(text=paste("cbind(",data,"[which(",data,"$",x,"==1 & ",data,"$",index2,"==1 ),],color=col[1])",sep="")))# hot spot
    w5 <- eval(parse(text=paste("cbind(",data2,"[which(",data2,"$",x,"==1 & ",data2,"$",index2,"==1),],color=col[2])",sep="")))# hot spot
    w10 <- eval(parse(text=paste("cbind(",data3,"[which(",data3,"$",x,"==1 & ",data3,"$",index2,"==1),],color=col[3])",sep="")))# hot spot
    Nothing <- eval(parse(text=paste("cbind(",data3,"[which(",data3,"$",x,"==0 | ",data3,"$",index2,"==0),],color='#6BAACC')",sep="")))# hot spot
  }
  
  
  if(nbvar==1){
    w <- eval(parse(text=paste("cbind(",data,"[which(",data,"$",x,"==1 ),],color=col[1])",sep="")))# hot spot
    w5 <- eval(parse(text=paste("cbind(",data2,"[which(",data2,"$",x,"==1),],color=col[2])",sep="")))# hot spot
    w10 <- eval(parse(text=paste("cbind(",data3,"[which(",data3,"$",x,"==1) ,],color=col[3])",sep="")))# hot spot
    Nothing <- eval(parse(text=paste("cbind(",data3,"[which(",data3,"$",x,"==0),],color='#6BAACC')",sep="")))# hot spot
  }
  
  plot(rnorm(1000),col="white",bg="grey",axes=F,xlab="",ylab="",xlim=xlim,ylim=ylim)
  rect(xleft =-180,ybottom=-90,xright=180,ytop=90,density=NULL,angle=45,col="#B6E0EE",border=NULL)
  points(x= Nothing[,"X"],y=Nothing[,"Y"],col=as.character(Nothing$color),pch=15,cex=cex_point)
  points(x=w10[,"X"],y=w10[,"Y"],col=as.character(w10$color),pch=15,cex=cex_point)
  points(x=w5[,"X"],y=w5[,"Y"],col=as.character(w5$color),pch=15,cex=cex_point)
  points(x=w[,"X"],y=w[,"Y"],col=as.character(w$color),pch=15,cex=cex_point)
  abline(h=0,lty="dotted",col="grey30")
  
  plot(coast,col="black",add=TRUE,lwd=0.02)    
  if(axisx==T){
    axis(side=1,line=0,cex.axis=0.5,lwd=0.35,tcl=-0.25,bg="white",labels=c("-180°","-135°","-90°","-45°","0°","45°","90°","135°","180°"),
         at=c(-179,-135,-90,-45,0,45,90,135,179),mgp=c(3,0.05,0))
  }
  
  if(axisy==T){
    axis(side=2,line=0,las=2,cex.axis=0.5,lwd=0.35,tcl=-0.25,bg="white",labels=c("-75°","-50°","-25°","0°","25°","50°","75°"),
         at=c(-75,-50,-25,0,25,50,75),mgp=c(3,0.25,0))
    
  }
  
  mtext(names_fig,side=2,line=pos_legX,at=pos_legY,cex=cex_text_fig,bg="white",las=2)
  colour=c("#B6E0EE","#6BAACC",col_pt)
  if(legend==TRUE){
    legend(pos_leg,leg=c("No values",text),col="black",pt.bg=colour,pch=22,cex=cex_legend,lty=NA,
           lwd=0.2,box.lwd=0.6,pt.cex=pt.cexleg,ncol=ncol,bg="#FFFFFF95",x.intersp=0.1,...)
  }
  box(lwd=0.6)
  
}# end of function carto_congruence


###################################################################################
#function plot
carto_function <- function(y=Hotspot,x="RS_tot",index2="FRiC_new",names_fig="a",pos_legY=87,pos_legX=0.5,
                           cex_point=0.3,cex_text_fig=1,cex_legend=0.5,cex.axis=0.4,pt.cexleg=1,
                           text=c("Non hotspots values for both variables","Values in top 2.5% for FRIC",
                                  "Values in top 2.5% for species richness","Congruence zone"),
                           xlim=c(-166,166),pos_leg="topright",ylim=c(-83.5,83.5),ncol=2,legend=TRUE,...){
    
  w <- eval(parse(text=paste("cbind(y[which(y$",x,"==1 & y$",index2,"==1),],color='#fdae61')",sep="")))# congruence fdae61 orange 
  yy <-eval(parse(text=paste("cbind(y[which(y$",x,"==0 & y$",index2,"==1),],color='#B70F0D')",sep="")))# hotspot FRic #B70F0D rouge
  z <- eval(parse(text=paste("cbind(y[which(y$",x,"==1 & y$",index2,"==0),],color='#EDDF19')",sep=""))) # hot spot richness # jaune
  x <- eval(parse(text=paste("cbind(y[which(y$",x,"==0 & y$",index2,"==0),],color='#6BAACC')",sep=""))) # rien

	plot(rnorm(1000),col="white",bg="grey",axes=F,xlab="",ylab="",xlim=xlim,ylim=ylim)
	rect(xleft =-180,ybottom=-90,xright=180,ytop=90,density=NULL,angle=45,col="#B6E0EE",border=NULL)

  points(x$X,x$Y,col=as.character(x$color),pch=15,cex=cex_point)
 	points(yy$X,yy$Y,col=as.character(yy$color),pch=15,cex=cex_point)
	points(w$X,w$Y,col=as.character(w$color),pch=15,cex=cex_point)
	points(z$X,z$Y,col=as.character(z$color),pch=15,cex=cex_point)
 
	plot(coast,col="black",add=T,lwd=0.02)   
	abline(h=0,lty="dotted",col="grey30")
	axis(side=1,line=0,cex.axis=0.5,lwd=0.35,tcl=-0.25,bg="white",labels=c("-180°","-135°","-90°","-45°","0°","45°","90°","135°","180°"),
	     at=c(-179,-135,-90,-45,0,45,90,135,179),mgp=c(3,0.05,0))
	axis(side=2,line=0,las=2,cex.axis=0.5,lwd=0.35,tcl=-0.25,bg="white",labels=c("-75°","-50°","-25°","0°","25°","50°","75°"),
	     at=c(-75,-50,-25,0,25,50,75),mgp=c(3,0.25,0))
  mtext(names_fig,side=2,line=pos_legX,at=pos_legY,cex=cex_text_fig,bg="white",las=2)
  box(lwd=0.6)
 
  colour=c("#B6E0EE","#6BAACC","#B70F0D","#EDDF19","#fdae61")
  if(legend==TRUE){
  legend(pos_leg,leg=c("No values",text),col="black",pt.bg=colour,pch=22,cex=cex_legend,lty=NA,
         lwd=0.2,box.lwd=0.6,pt.cex=pt.cexleg,ncol=ncol,bg="#FFFFFF95",x.intersp=0.1,...)
         }
    
}# end of function carto_function()

###################################################################################
Map_the_world_small <- function(Data_PA=grid$P_A, coord_X = grid$X_ENTIER,coord_Y=grid$Y_ENTIER,
                                col=c("White","yellow","red3"),breaks=seq(0,40,2),include.lowest=T,
                                xlim=c(-166,166),ylim=c(-83.5,83.5),zlevels=6,coastline=T,
                                main_legend="Species Richness",output=F,legend_increment=100,names_fig="a)",
                                cex_names_fig=1,xpt=-116.5,ypt=30.5,cex.pt=0.3,legend=T,main="USA",...){
    require(shape)
     
    if (output==TRUE) {
         	tiff(filename = "world_map.tiff",width = 6.5, height = 4.5, units = "in", res=300,compression = "lzw")
         
        #par(mar=c(1, 2, 2,5))
         
         	plot(rnorm(1000),col="white",type="n",axes=F,xlab="",ylab="",xlim=xlim,ylim=ylim,main=main)
         
         	jet.color <-  colorRampPalette(col)
         	colour <-  jet.color(length(breaks)-1)
         
         	richesse <-  cut(Data_PA , breaks=breaks,include.lowest=include.lowest)
         
         	points(coord_X,coord_Y,col = colour[richesse], pch = 15,cex = cex.pt) 
         
         		if (coastline == TRUE) {
             		plot(coast,col="black",add=TRUE)
             		} # end of if
         		
         	axis(side=1,line=0,cex.axis=0.7,lwd=0.35,tcl=-0.25,bg="white",labels="",at=seq(xlim[1],xlim[2],20),mgp=c(3, 0.5, 0))
         	#axis(side=2,line=0,las=2,cex.axis=0.7,lwd=0.35,tcl=-0.25,bg="white",labels="",at=c(-90,-45,0,45,90),mgp=c(3, 0.5, 0))
         
         #mtext(names_fig,side=2,line=0.5,at= 100,cex=cex_names_fig,bg="white",las=2)  	
        
        box()
         
         	#colorlegend(zlim=round(c(0,max(Data_PA))), dz=legend_increment,col=colour, zlevels=zlevels,main=main_legend,xpd=FALSE,main.cex=0.6,cex=0.7,posx=c(0.87, 0.9))
         dev.off()
         
         } else {
             	
             #par(mar=c(1, 2, 2,5))
             
             	plot(rnorm(1000),col="white",type="n",axes=F,xlab="",ylab="",xlim=xlim,ylim=ylim,main=main)
             	
             	jet.color <-  colorRampPalette(col)
             	colour <-  jet.color(length(breaks)-1)
             
             	richesse <-  cut(Data_PA , breaks=breaks,include.lowest=include.lowest)
             
             	points(coord_X,coord_Y,col = colour[richesse], pch = 15,cex =  cex.pt) 
             
             		if (coastline == TRUE) {
                 		plot(coast,col="black",add=TRUE)
                 		} # end of if

                    points(x=xpt,y=ypt,cex=3,pch=1,col="red",lwd=1.5)
                    points(x=xpt,y=ypt,cex=3,pch=1,col="white",lwd=2.5)
             
             	axis(side=1,line=0,cex.axis=0.8,lwd=0.35,tcl=-0.25,bg="white",labels=seq(xlim[1],xlim[2],20),at=seq(xlim[1],xlim[2],20),mgp=c(3, 0.5, 0))
                axis(side=2,line=0,las=2,cex.axis=0.8,lwd=0.35,tcl=-0.25,bg="white",labels=seq(ylim[1],ylim[2],20),at=seq(ylim[1],ylim[2],20),mgp=c(3, 0.5, 0))
             
             #mtext(names_fig,side=2,line=0.5,at= 100,cex=cex_names_fig,bg="white",las=2)  		
             box()
             
             if(legend==TRUE){
               colorlegend(zlim=round(c(min(Data_PA),max(Data_PA))), dz=legend_increment,col=colour,zlevels=zlevels, main=main_legend,xpd=FALSE,main.cex=0.9,cex=0.65,posx=c(0.85, 0.89))
             }

             
             } # end of ifelse
     } # end of function Mpa_the_world
 
###################################################################################
 get_pcoa_graph <- function(grid=grille_MAR_MAMM,data=pcoa_coord,data_func_glob=FE_vertices_global$d2D,
                          data_func_loc= vertices_Axis23_MM ,id=43755,ax=c(1,2),xlab="PC 1",ylab="PC 2",
                          colin="#F9E2A0"){
   
   vert <- data[data_func_glob,]
   vert_1 <- rbind(vert,vert[1,])
   
   vert_loc <- data[data_func_loc[[id]],]
   vert_loc_1 <- rbind(vert_loc,vert_loc[1,])
   
   plot(rnorm(100),col="white",ylim=c(-0.4,0.3),xlim=c(-0.4,0.35),type="n",axes=F,xlab=xlab,ylab=ylab,mgp=c(1.25,0.25,0.5),font.lab=2)
   polygon(data.frame(vert_1[,ax[1]],vert_1[,ax[2]]),col="grey88")
   
   points(vert_loc[,ax[1]],vert_loc[,ax[2]],col=colin,pch=19)
   points(vert_loc[,ax[1]],vert_loc[,ax[2]],col="black",pch=1)
   
   polygon(data.frame(vert_loc_1[,ax[1]],vert_loc_1[,ax[2]]),col=colin)
   
   points(x=data[data_func_glob,ax[1]],y=data[data_func_glob,ax[2]],xlab="",pch=16,col="black")
   
   data_loc <- pcoa_coord[names(grid[as.character(id),which(grid[as.character(id),]==TRUE)]),]
   
   points(data_loc[,ax[1]],data_loc[,ax[2]],pch=21,col=colin,lty=2)
   points(data_loc[,ax[1]],data_loc[,ax[2]],pch=1,col="black",lty=2)
   
   axis(side=1,cex.axis=0.8,at = seq(-0.5,0.5,0.2),lwd=0.35,tcl=-0.25,bg="white",labels=as.character(seq(-0.5,0.5,0.2)),mgp=c(0,0.25,0))
   axis(side=2,cex.axis=0.8,at = seq(-0.5,0.5,0.2),lwd=0.35,tcl=-0.25,bg="white",labels=as.character(seq(-0.5,0.5,0.2)),mgp=c(0,0.25,0))
   
   box(lwd=0.55)
 } # end of function

###################################################################################

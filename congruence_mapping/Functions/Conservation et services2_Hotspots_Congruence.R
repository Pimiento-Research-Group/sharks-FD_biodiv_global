################FONCTION Cons_Serv_Permut##################################################################
# Principe: Pour chaque variable, tri des 5% valeurs les plus fortes qui sont remplacées par 1 (hotspots),#
# les autres valeurs sont remplacées par 0                                                                #
#  Calcul du nombre de hotspots communs (observés) entre variables prises 2 à 2                           #
#  Variables prise 2 à 2, permutations d'une des 2 variables et calcul du nombre de hotspots communs      #                                                                                               #
#  Calcul de la p value                                                                                   #
#                                                                                                         #
# Input: matrice objets (e.g. sites, cellules) x variables (e.g. réserves, endémiques, pressions)         #
#        nombre de permutations désirées (nbPermut)                                                       #
#                                                                                                         #
# Output:                                                                                                 #
#    - automatique: nombre de lignes, nombre de colonnes, nombre hotspots observés par variable, tableau  #
# résultats (nombre double hotspots observés, nombre double hotspots attendu, nb theorique double hotspots#
# >= observé, p value                                                                                     #
#    - à la demande: $Mat.Hotspot donne la matrice binaire des hotspots observés                          #
###########################################################################################################
         
Cons_Serv_Permut=function(x,nbPermut,seuil=10/100){

	library(gtools)

	#####################################################################################################
	##Fonction donnant la plus faible valeur des z plus fortes valeurs (hotspots) de chaque variable###  
	fonc5=function(y,z){
		oo=order(y,decreasing=TRUE)
		x2=y[oo]
		val50=x2[z]
	}# end fonction fonc5
	#####################################################################################################

	nblgn=dim(x)[1]#nombre lignes du fichier d'origine
	nbcol=dim(x)[2]#nombre colonnes du fichier d'origine

	#Nombre lignes= seuil nombre lignes du fichier d'origine
	Nb_ligne=ceiling(nblgn*seuil)

	x01=x
	#Remplacement, pour chaque variable, des valeurs des hotspots par 1 et par 0 pour les non hotspots####
	for(i in 1:nbcol){

		v=fonc5(x[,i],Nb_ligne)

		x01[x01[,i]<v,i]=0
		x01[x01[,i]!=0,i]=1

	}# end for i
	#######################################################################################################

	#Combinaisons des n° de colonnes prises 2 à 2
	comb=combinations(nbcol,2,c(1:nbcol))
	#Nombre de combinaisons
	comblgn=dim(comb)[1]
	ResHot=matrix(0,comblgn,7,dimnames=list(NULL,c("Var1","Var2","Nb.hotspots.communs.Obs","Espere=NiNj/Nt","Theo>=Obs","p_value","Sign")))
	ResHot=as.data.frame(ResHot)

	for (j in 1:comblgn){
		#Extraction des colonnes 2 à 2
		xV1V2=x01[,c(comb[j,1],comb[j,2])]
		#Calcul valeur esperee
		ninj=apply(xV1V2,2,sum)
		Esp=(ninj[1]*ninj[2])/nblgn
		#Nombre de hotspots communs
		concord.hot=sum(xV1V2[,1]*xV1V2[,2])
		#Affectations valeurs dans tableau resultats
		ResHot[j,1:2]=colnames(xV1V2)
		ResHot[j,3]=concord.hot
		ResHot[j,4]=Esp
		Att=(concord.hot-Esp)/Esp

		vHot=vector("numeric")
			##TEST######################################################################################################
			for(k in 1:nbPermut){
				permut=sample(xV1V2[,1])
				hot.com=sum(permut*xV1V2[,2])#Nb hotspots communs entre 2 variables
				vHot[k]=hot.com

			}# end for k
		TsupO=length(vHot[vHot>=concord.hot])#Nb permutations où hotspots communs >= nb hotspots communs observe
		pValue=(TsupO+1)/(nbPermut+1)
		############################################################################################################
		ResHot[j,5:6]=c(TsupO,pValue)

	}# end for j

	ResHot[,7]="NS"
	ResHot[ResHot[,6]<0.05|ResHot[,6]>0.95,7]="*"
	ResHot[ResHot[,6]<0.01|ResHot[,6]>0.99,7]="**"
	ResHot[ResHot[,6]<0.001|ResHot[,6]>0.999,7]="***"

	nbhot=apply(x01,2,sum)#nombre hotspots par variable
	res=list(Nb.lignes=nblgn,Nb.colonnes=nbcol,Nb.Hotspots=nbhot,Resultats=ResHot,Mat.Hotspot=x01)
	print(res[1:4])
	return(res)

}# end fonction Cons_Serv_Permut

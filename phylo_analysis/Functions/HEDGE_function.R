### HED2: function to calculate HEDGE, LEDGE and HED scores 
### which also returns the pruned tree used for calculating scores and the scores by node

HED2 <- function(phyl, proba, subtree=FALSE, tol=1e-8){
  
  # phyl is an object of class hclust, phylo or phylo4
  # proba is a vector with names corresponding to the tips of tree
  # the function itself put the species in the same order in the phylogenetic tree and in proba
  # If proba does not contain all tips of phyl, then phyl is pruned and a new root is defined at the most recent common ancestor of the retained species as the root of the pruned tree. If subtree is FALSE, the branch between the new root and the root of phyl is kept; if subtree is TRUE, it is removed. 
  if(is.null(names(proba)) | any(is.na(names(proba)))) stop("Vector proba must have names corresponding to the tips of the phylogenetic tree")
  nsp <- length(proba)
  if(nsp==1) stop("At least two tips should be selected in proba")
  arg.phyl <- .checkphyloarg(phyl)
  tre <- arg.phyl$phyl.phylo
  if(!is.rooted(tre)) stop("phyl must be a rooted phylogenetic tree")
  if(any(proba < -tol) | any(is.na(proba))) stop("Values in proba must be known and positive")
  if(!all(names(proba)%in%tre$tip.label)) stop("missing tips in the phylogenetic tree in comparison with names of vector proba")
  if(subtree){
    tre <- drop.tip(tre, tip=tre$tip.label[!tre$tip.label%in%names(proba)])
  } else
  {
    tre <- drop.tip(tre, tip=tre$tip.label[!tre$tip.label%in%names(proba)], root.edge = tre$Nnode)
  } 
  tre4 <- as(tre, "phylo4")
  if(!hasNodeLabels(tre4)){
    nodeLabels(tre4) <- names(nodeLabels(tre4))
  } else{
    e <- nodeLabels(tre4)
    e[is.na(e)] <- names(e[is.na(e)])
    nodeLabels(tre4) <- e
  }
  proba <- proba[tipLabels(tre4)]
  if(nNodes(tre4)==1){
    des <- descendants(tre4, nodeLabels(tre4), type="tips")
    prob <- prod(proba[names(des)])
  } else{
    des <- descendants(tre4, nodeLabels(tre4), type="tips")
    prob <- unlist(lapply(des, function(x) prod(proba[x])))
  }
  valnodes <- prob*edgeLength(tre4, nodeLabels(tre4))
  valtips <- proba*edgeLength(tre4, names(proba))
  names(valtips) <- names(proba)
  names(valnodes) <- nodeLabels(tre4)
  vals <- c(valnodes, valtips)
  # score for the tips (= the species)
  hedge <- unlist(lapply(ancestors(tre4, tipLabels(tre4), "ALL"), function(x) sum(vals[names(x)], na.rm=TRUE)))
  ledge <- hedge/proba[names(hedge)]*(1-proba[names(hedge)])
  hed <- hedge + ledge
  res <- cbind.data.frame(HEDGE=hedge, LEDGE=ledge, HED=hed)
  rownames(res) <- names(hedge)
  # scores for all the nodes (including tips)
  Nhedge <- unlist(lapply(ancestors(tre4, getNode(tre4), "ALL"), function(x) sum(vals[names(x)], na.rm=TRUE)))
  Nproba <- c(proba, prob)
  Nledge <- Nhedge/Nproba*(1-Nproba)
  Nhed <- Nhedge + Nledge
  Nres <- cbind.data.frame(HEDGE=Nhedge, LEDGE=Nledge, HED=Nhed)
  rownames(Nres) <- names(Nhedge)
  # list of outputs
  reslist <- list(tre4, res, Nres)
  names(reslist) <- c("pruned.tree", "scores", "node.scores")
  return(reslist)
}
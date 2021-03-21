# lengthRatio
length_ratio <- function (tree.phylo,d) {
 
  leaves <- 1:nrow(d)
  Nleaves <- length(leaves)
  nodes  <- (Nleaves+1):(Nleaves+tree.phylo$Nnode)
  #internal edges 
  Iedges <- which(tree.phylo$edge[,1]>Nleaves & tree.phylo$edge[,2]>Nleaves)
  Tedges<-matrix(0,length(Iedges),3)
  #compute for each Iedge X1, X2
  for (i in 1:length(Iedges)) {
    node.1 <- tree.phylo$edge[Iedges[i],2]
    node.2 <- tree.phylo$edge[Iedges[i],1]
    # rk for each internal edge, node 1 (cd) < node 2 (ab)
    # ie. node.ab is the node with the smallest size, i.e.the closest to a leaf
    X1 <- getDescendants(tree.phylo, node=node.1, curr=NULL)
    X1 <- sort(X1[which(X1<(Nleaves+1))],decreasing=FALSE)
    X2 <- setdiff(1:Nleaves,X1)
    X2 <- sort(X2,decreasing=FALSE)
    
    num=mean(as.matrix(d[X1,X2]))
    inX1=d[X1,X1]
    inX1=inX1[upper.tri(inX1)]
    inX2=d[X2,X2]
    inX2=inX2[upper.tri(inX2)]
    denom=(sum(inX1)+sum(inX2))/(length(inX1)+length(inX2))
    Tedges[i,1] <- num/denom
    Tedges[i,2:3]<-c(length(X2),length(X1))
  }
  Idedges=tree.phylo$edge[Iedges,]
  if (length(Iedges)==1)  Idedges=t(Idedges)
  lengthratio<-cbind(Idedges,Tedges)
  return(lengthratio)
}
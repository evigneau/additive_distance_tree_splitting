# proportion of well-designed triples of objects/items
prop_wd_trip <- function (tree.phylo,d) {
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
    Tedges[i,1] <- trip_elemen_2(d,X1,X2)
    Tedges[i,2:3]<-c(length(X2),length(X1))
  }
  propwdtrip<-cbind(tree.phylo$edge[Iedges,],Tedges)
  return(propwdtrip)
}

expand.grid.df <- function(...) Reduce(function(...) merge(..., by=NULL), list(...))

trip_elemen_2<-function(d,A,B){
  
  cA=t(combn(A,2,simplify=TRUE))
  cB=t(combn(B,2,simplify=TRUE))
  combinaisonsAB=expand.grid.df(A,cB)    
  combinaisonsBA=expand.grid.df(B,cA)   
  
  mat_areteAB<-matrix(ncol=4,nrow=nrow(combinaisonsAB))
  for (i in 1:nrow(mat_areteAB)){
    a<-combinaisonsAB[i,1]
    b1<-combinaisonsAB[i,2]
    b2<-combinaisonsAB[i,3]
    C1<-0
    C1<-(d[b1,b2]<min(d[a,b1],d[a,b2]) )
    mat_areteAB[i,]<-c(a,b1,b2,C1)
  }
  mat_areteBA<-matrix(ncol=4,nrow=nrow(combinaisonsBA))
  for (i in 1:nrow(mat_areteBA)){
    b<-combinaisonsBA[i,1]
    a1<-combinaisonsBA[i,2]
    a2<-combinaisonsBA[i,3]
    C1<-0
    C1<-(d[a1,a2]<min(d[a1,b],d[a2,b]))
    mat_areteBA[i,]<-c(a,b1,b2,C1)
  }
  
  taux_triplets<-sum(c(mat_areteAB[,4],mat_areteBA[,4])==1)/length(c(mat_areteAB[,4],mat_areteBA[,4]))*100 
  return(taux_triplets)
}
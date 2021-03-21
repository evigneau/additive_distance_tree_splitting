# search for subtrees using binary splitting

dichotomize <- function(desc,d,crit,crit.ratiodiam,crit.siz) {

   dbis=d[desc,desc]

   ############ remark :if there are 3 or less objets => no internal edge => length_ratio can not be evaluated
  if (nrow(dbis) >3) {
    tree.phylo.bis<-ape::nj(as.dist(dbis))
    #plot(tree.phylo.bis)
    #edgelabels(tree.phylo.bis$edge[,1],frame="none",cex=0.8)
    if (crit=="lengthratio")  matcrit<-length_ratio(tree.phylo.bis,dbis)
    if (crit=="proptriplets") matcrit<-prop_wd_trip(tree.phylo.bis,dbis)
    # if (crit=="lengthedge")   matcrit<-cbind(tree.phylo.bis$edge,tree.phylo.bis$edge.length)[Iedges,]
    criterion=matcrit
    criterion<-criterion[order(criterion[,3],decreasing=TRUE),] # ordered criterion
    if(is.vector(criterion)) criterion=t(criterion)
    step=1
    if (criterion[step,1]<criterion[step,2]) {
      n.g <- criterion[step,2]
      n.d <- criterion[step,1]
    }    else {
      n.g <- criterion[step,1]
      n.d <- criterion[step,2]
    }
    desc.g <- phytools::getDescendants(tree.phylo.bis,n.g)
    desc.g=sort(desc.g[which(desc.g<(length(tree.phylo.bis$tip.label)+1))],decreasing=FALSE)
    num.g=desc[desc.g]
        # print("=====================================")
        # print(tree.phylo.bis$tip.label[desc.g])
        # print("-------------------------------------")
    desc.d=setdiff(1:length(tree.phylo.bis$tip.label),desc.g)
    num.d=desc[desc.d]
        # print(tree.phylo.bis$tip.label[desc.d])
        # print("=====================================")

    splitting.status=c(TRUE,TRUE)
    dg=d[num.g,num.g]
    dg=dg[upper.tri(dg,diag=FALSE)]
    if((max(dg)/max(d[upper.tri(d,diag=FALSE)])<crit.ratiodiam)|(length(num.g)<=crit.siz)) splitting.status[1]=FALSE
    dd=d[num.d,num.d]
    dd=dd[upper.tri(dd,diag=FALSE)]
    if((max(dd)/max(d[upper.tri(d,diag=FALSE)])<crit.ratiodiam)|(length(num.d)<=crit.siz)) splitting.status[2]=FALSE
    rdiam=c(max(dg)/max(d[upper.tri(d,diag=FALSE)]),max(dd)/max(d[upper.tri(d,diag=FALSE)]))
    siz=c(length(num.g),length(num.d))
  }  # end if nrow(dbis)>3

  if (nrow(dbis)==3) {
    pos=which(dbis[lower.tri(dbis)]==min(dbis[lower.tri(dbis)]))
    num.g=which((lower.tri(dbis)==TRUE) ,arr.ind=TRUE)[pos,]
    num.d=setdiff(1:3,num.g)
    dg=dbis[num.g[1],num.g[2]]
    dd=0
    num.g=desc[num.g]
    num.d=desc[num.d]
    splitting.status=c(TRUE,FALSE)
    if(dg/max(d[upper.tri(d,diag=FALSE)])<crit.ratiodiam) splitting.status[1]=FALSE
    rdiam=c(dg/max(d[upper.tri(d,diag=FALSE)]),0)
    siz=c(length(num.g),length(num.d))
  }
  if (nrow(dbis)==2) {
    num.g=desc[1]
    num.d=desc[2]
    splitting.status=c(FALSE,FALSE)
    rdiam=c(0,0)
    siz=c(1,1)
  }

  return(list(num.g=num.g,num.d=num.d,splitting.status=splitting.status,rdiam=rdiam,siz=siz))
}

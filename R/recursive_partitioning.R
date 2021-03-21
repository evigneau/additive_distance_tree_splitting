#' recursive partitioning
#'
#' @param d : distance matrix for the whole set of objects/items.
#' @param nameobj : the label for the objects/items. If NULL (by default) dimnames of d are used.
#' @param crit : criterion for splitting, either "lengthratio" (by default) either "proptriplets")
#' @param crit.ratiodiam : above the provided threshold the subtree is candidate for splitting. \cr
#'                         The criterion considered is the ratio of the diameter of the subset in a subtree \cr
#'                        to diameter of the whole sets of objects (relative diameter). \cr
#'                        The diameter is defined as the maximum pairwise distance within the set of objects; \cr
#' @param crit.siz : size of a subtree (in absolute value), \cr
#'                   greater (>) to this threshold (3, by default), the subtree is candidate for splitting.
#'
#' @return \item{partitions}{ Results of the partitions at each step}
#'
#' @references Koenig, L., Cariou, V., Symoneaux, R., Coulon-Leroy, C., Vigneau, E. (2021). Food Quality and Preference, 89, 104137(2021).
#'
#'@import ape
#'@import phytools
#'@importFrom stats as.dist
#'@importFrom utils combn
#'
#' @examples data(DistAnimal)
#' tree=ape::nj(DistAnimal)
#' ape::plot.phylo(tree,type="unrooted")
#' recursive_partitioning(d=DistAnimal)
#'
#' @export


recursive_partitioning<-function(d,nameobj=NULL, crit="lengthratio",crit.ratiodiam=0.6,crit.siz=3)
{

  n=nrow(as.matrix(d))
  partitions=list()
  resav= list()
  if (is.null(nameobj)) {
    desc=dimnames(as.matrix(d))[[1]]
  } else {
    desc=nameobj
  }

  ngp=0
  # First step
  s=1
  partitions[[s+1]]=list()
  cas=0
  res=dichotomize(desc,d=as.matrix(d),crit=crit,crit.ratiodiam=crit.ratiodiam,crit.siz=crit.siz)
  partitions[[s]][[ngp+1]]=res$num.g
  partitions[[s]][[ngp+2]]=res$num.d
  if (res$splitting.status[1]==FALSE) {
    cas=cas+1
    partitions[[s+1]][[ngp+cas]]=res$num.g
  }
  if (res$splitting.status[2]==FALSE) {
    cas=cas+1
    partitions[[s+1]][[ngp+cas]]=res$num.d
  }
  ngp=ngp+cas
  resav[[s]]=res
  allsplittingstatus=resav[[s]]$splitting.status
  allrdiam=resav[[s]]$rdiam;       #print(allrdiam)
  allsiz=resav[[s]]$siz;           #print(allsiz)

  # next steps, until the splitting status is FALSE for each subtree
  while (sum(allsplittingstatus==FALSE)<length(allsplittingstatus)) {
    s=s+1;
    partitions[[s+1]]=list()
    inc=0
    L=length(resav);
    for (i in 1:L) {
      if (resav[[i]]$splitting.status[1]) {
        desc=resav[[i]]$num.g
        res=dichotomize(desc,d=as.matrix(d),crit=crit,crit.ratiodiam=crit.ratiodiam,crit.siz=crit.siz)
        partitions[[s]][[inc+1]]=res$num.g
        partitions[[s]][[inc+2]]=res$num.d
        inc=inc+2
        resav[[length(resav)+1]]=res
      } else {
        partitions[[s]][[inc+1]]=resav[[i]][[1]]
        inc=inc+1
        if (length(resav[[i]])>(2+2)) {
          tempo=resav[[i]][-2]
          tempo$splitting.status=tempo$splitting.status[-2]
          tempo$rdiam=tempo$rdiam[-2]
          tempo$siz=tempo$siz[-2]
          resav[[length(resav)+1]]=tempo
        } else {
          resav[[length(resav)+1]]=resav[[i]]
        }
      }

      if (length(resav[[i]])>(2+2)) {
        if (resav[[i]]$splitting.status[2]) {
          desc=resav[[i]]$num.d
          res=dichotomize(desc,d=as.matrix(d),crit=crit,crit.ratiodiam=crit.ratiodiam,crit.siz=crit.siz)
          partitions[[s]][[inc+1]]=res$num.g
          partitions[[s]][[inc+2]]=res$num.d
          inc=inc+2
          resav[[length(resav)+1]]=res
        } else {
          partitions[[s]][[inc+1]]=resav[[i]]$num.d
          inc=inc+1
          if (length(resav[[i]])>2) {
            tempo=resav[[i]][-1]
            tempo$splitting.status=tempo$splitting.status[-1]
            tempo$rdiam=tempo$rdiam[-1]
            tempo$siz=tempo$siz[-1]
            resav[[length(resav)+1]]=tempo
          } else {
            resav[[length(resav)+1]]=resav[[i]]
          }
        }}
    }  # fin boucle i
    resav=resav[-(1:L)]
    allsplittingstatus=NULL
    allrdiam=NULL
    allsiz=NULL
    for (i in 1:length(resav)) {
      allsplittingstatus=c(allsplittingstatus,resav[[i]]$splitting.status)
      allrdiam=c(allrdiam,resav[[i]]$rdiam)
      allsiz=c(allsiz,resav[[i]]$siz)
    }
     #print(allsplittingstatus);  print(allrdiam); print(allsiz)
  } # end while

  partitions=partitions[-(s+1)]
  return(partitions)
}




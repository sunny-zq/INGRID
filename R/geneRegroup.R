#' @title geneRegroup function
#'
#' @description  redefining the gene set membership of common genes. First, the gene remains to be as a member if it is a core member of the pathway. Second, if the gene is not a core member of any among the pathways (i.e., it maps to more than one pathways), then this gene is re-assigned to one of the gene sets based on the k-medoids algorithm minimizing the binary distance between genes within cluster distance.
#'
#' @export
#' @import foreach cluster methods
#' @param plist user provided gene lists organized in pathways
#'
#'
#' @examples
#' data(TCGA)
#' geneRegroup.result=geneRegroup(plist=TCGA$pathList )
#' prefilter.results=prefilter( data=TCGA$geneexpr, time=TCGA$t, status=TCGA$d,
#'                              plist=geneRegroup.result@gset )



geneRegroup<- function(plist=TCGA$pathList){

  # names is the unique genes in the plist pathways
  names<-unique(as.vector(unlist(plist)))

  pmat<-foreach(i=1:length(names(plist)), .combine='cbind')%do%{
    ifelse( names%in%plist[[i]], 1, 0 )
  }

  rownames(pmat)=names
  smat<-pmat[rowSums(pmat)>1,]
  umat<-pmat[rowSums(pmat)==1,]

  colnames(pmat)=colnames(smat)=colnames(umat)=names(plist)

  paths<-apply(umat,1,function(x){which(x==1)})


  ##### Added          ############################################################
  ##### the following step solved the issue of there is no unique gene to a pathway

  temp1<- as.data.frame(table(paths))
  temp2<- foreach(i=1:length(plist),.combine = "c")%do%{
    ifelse( (c(1:length(plist))%in%temp1$paths)[i], temp1$Freq[sum((c(1:length(plist))%in%temp1$paths)[1:i])], NA )
  }
  mat<-cbind(sapply(plist,length),temp2)

  ################################################################################
  colnames(mat)<-c("total genes","unique genes")

  d=as.dist(dist(smat, "binary"))
  si=rep(0,10)
  for(k in 2:10){
    pam_k<-pam(d, k,diss=TRUE)
    si[k]<-summary(silhouette(pam_k))$avg.width
  }

  pam.fit<-pam(d, 2, diss=TRUE)
  plot(pam.fit)
  gset1<-gset2<-list()
  for(i in 1:length(pam.fit$medoids)){
    idx<-as.vector(pam.fit$clustering)
    gset1[[i]]<-names(pam.fit$clustering)[idx==i]
  }
  names(gset1)<-pam.fit$medoids

  for(i in 1:length(names(plist))){
    idx<-as.vector(paths)
    gset2[[i]]<-names(paths)[idx==i]
  }
  names(gset2)<-names(plist)

  gset<-c(gset1, gset2)
  gset<-gset[lapply(gset,length)>0]

  methods::new( "RegroupGene",
                gset = gset,
                inputdata = plist )
}

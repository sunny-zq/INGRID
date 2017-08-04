#' @title Gene Regrouping summary
#'
#' @description Summary of gene Regrouping results.
#'
#' @aliases show,geneRegroup-method
#' @param object output of geneRegroup function

setMethod(
  f="show",
  signature="RegroupGene",
  definition=function( object ) {

    plist<-object@inputdata
    geneset.before <- names(plist)
    n.path.before <- length(geneset.before)
    path.index <- sprintf("pathway_%s",seq(1:n.path.before))

    path.before<-as.data.frame(names(object@inputdata))
    rownames(path.before)<-path.index
    colnames(path.before)<-"Pathway Name"

    n_newset<- sum(ifelse( names(object@gset)%in%names(plist), 0, 1 ))
    idx <- which(ifelse( names(object@gset)%in%names(plist), 0, 1 )==1)
    table<- matrix(rep(0,n_newset*n.path.before),ncol=n.path.before)
    for (j in 1:n_newset){
      names<- unlist(object@gset[idx[j]])
      pmat<-sapply(1:n.path.before,function(i){ifelse( names%in%plist[[i]], 1, 0 )})
      table[j,]<-colSums(pmat)
    }

    colnames(table)<-path.index
    rownames(table)<-names(object@gset)[idx]
    cat( "Summary: Gene-Regrouping results (class: RegroupGene)\n" )
    cat( "--------------------------------------------------\n" )
    cat( "Gene sets before gene regrouping\n" )
    str(object@inputdata)
    cat( "--------------------------------------------------\n" )
    cat( "Gene sets after gene regrouping\n" )
    str(object@gset)
    cat( "--------------------------------------------------\n" )
    cat( "Pathway names before gene regrouping\n" )
    print(path.before)
    cat( "--------------------------------------------------\n" )
    cat( "Number of genes in new gene set that are from the old pathway\n" )
    print(table)

  }
)

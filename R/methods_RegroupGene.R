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
    path.index <- sprintf("gene_set_%s",seq(1:n.path.before))

    path.before<-as.data.frame(names(object@inputdata))
    rownames(path.before)<-path.index
    colnames(path.before)<-"gene_set_name"

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
    cat( "Summary: Gene regrouping results (class: RegroupGene)\n" )
    cat( "--------------------------------------------------\n" )
    cat( "Gene sets before the gene regrouping\n" )
    str(object@inputdata)
    cat( "--------------------------------------------------\n" )
    cat( "Gene sets after the gene regrouping\n" )
    str(object@gset)
    cat( "--------------------------------------------------\n" )
    cat( "Gene set names before the gene regrouping\n" )
    print(path.before)
    cat( "--------------------------------------------------\n" )
    cat( "Number of genes in each of the new gene set vs. each of the original gene sets\n" )
    print(table)

  }
)

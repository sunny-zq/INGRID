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
    path.index <- sprintf("s%s",seq(1:n.path.before))
    
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
    for(i in 1:length(idx)){
      rownames(table)[i]<-paste0("gene_set_",nrow(path.before)+i)
    }
    
    cat( "Summary: Gene regrouping results (class: RegroupGene)\n" )
    cat( "--------------------------------------------------\n" )
    cat( "Gene sets before the gene regrouping\n" )
    str(object@inputdata)
    cat( "--------------------------------------------------\n" )
    cat( "Gene sets after the gene regrouping\n" )
    str(object@gset)
    cat( "--------------------------------------------------\n" )
    cat( "Comparison of the gene set before and after gene regrouping\n" )
    print(table)
    cat( "\n" )
    cat( "where\n" )
    print(path.before)
    
    
    
  }
)
#' @title prefilter function
#'
#' @description a supervised prefiltering by fitting a cox regression model on the expression measures of each mRNA.  Only genes with p values smaller than a pre-selected cut-off point were included in further analysis. Here we chose p=0.5 as default cut-off point.
#'
#' @export
#' @import foreach survival methods
#' @param data A N by P dataframe, where N is number of observations, P is number of genes
#' @param time time passed to selectGene's time argument.
#' @param status status selectGene's status argument.
#' @param p.cut p.cut is the pvalue cutoff point for prefiltering
#' @param plist output gene list from gene regrouping step
#'
#'
#' @examples
#' data(TCGA)
#' data(TCGA)
#' geneRegroup.result=geneRegroup( plist=TCGA$pathList )
#' prefilter.results=prefilter( data=TCGA$geneexpr, time=TCGA$t, status=TCGA$d,
#'                              plist=geneRegroup.result@gset )

prefilter <- function( data, time, status, p.cut=0.5, plist ){
  if (length(which(is.na(time)==TRUE))>0  ){
    stop( "There is missing value in the time variable. Cannot proceed!\n" )
  } else if  ( length(which(is.na(status)==TRUE))>0){
    stop( "There is missing value in the status variable. Cannot proceed!\n" )
  } else if(sum(apply(data,  2,  function(x){length(which(is.na(x)==TRUE)) })) >0 ){
    stop( "There is missing value in the expression data. Cannot proceed!\n" )
  } else {
    pvals <- foreach(i=1:ncol(data), .combine='c') %do% {
      if ( i%%100 == 0 ) message( paste("Preprocessing genes ",i,"/",ncol(data),"...",sep="") )

      summary(coxph( Surv(time, status)~data[,i] ))$coef[5]
    }
    message("Done!")

  # newdat=data[,pvals<=p.cut]

  xlist<- lapply( plist,function(x){
    idx<- which( (colnames(data)%in%x)==T & (pvals<=p.cut)==T )
    data[,idx,drop=FALSE]
  } )

  nc<-lapply(xlist,ncol)
  plist=plist[which(nc!=0)]
  xlist=xlist[which(nc!=0)]
  methods::new( "Prefiltered",
                xlist = xlist,
                inputdata = list( data = data, time = time, status = status, pathway = names(xlist) )
  )
  }
}

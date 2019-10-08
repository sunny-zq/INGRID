#' @title checkData summary
#'
#' @description Summary of checkData results.
#'
#' @aliases show,CheckedData-method
#'
#' @param object output of prefilter function

setMethod(
  f="show",
  signature="CheckedData",
  definition=function( object ) {
    data.before <- unlist(object@inputdata$data)
    data.after <-  unlist(object@newdata$data)
    # t<-object@inputdata$time
    # t.after<-object@newdata$t
    cat( "Summary: CheckData results (class: CheckedData)\n" )
    cat( "--------------------------------------------------\n" )
    cat("Summary of data before CheckData:\n")
    cat("          min:", round(min(data.before),2), "\n")
    cat(" 1st quantile:", round(quantile(data.before,0.01),2), "\n")
    cat("25th quantile:", round(quantile(data.before,0.25),2), "\n")
    cat("       median:", round(median(data.before),2), "\n")
    cat("75th quantile:", round(quantile(data.before,0.75),2), "\n")
    cat("99th quantile:", round(quantile(data.before,0.99),2), "\n")
    cat("          max:", round(max(data.before),2), "\n")
    cat("         mean:", round(mean(data.before),2), "\n")
    # cat("  sample size:", length(t), "\n")
    cat( "--------------------------------------------------\n" )
    cat("Summary of data after CheckData:\n")
    cat("          min:", round(min(data.after),2), "\n")
    cat(" 1st quantile:", round(quantile(data.after,0.01),2), "\n")
    cat("25th quantile:", round(quantile(data.after,0.25),2), "\n")
    cat("       median:", round(median(data.after),2), "\n")
    cat("75th quantile:", round(quantile(data.after,0.75),2), "\n")
    cat("99th quantile:", round(quantile(data.after,0.99),2), "\n")
    cat("          max:", round(max(data.after),2), "\n")
    cat("         mean:", round(mean(data.after),2), "\n")
    # cat("  sample size:", length(t.after), "\n")
    cat( "--------------------------------------------------\n" )

  }
)

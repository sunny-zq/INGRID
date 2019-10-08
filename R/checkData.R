#' @title checkData function
#'
#' @description provide EDA for the data
#'
#' @export
#' @import methods ggplot2
#' @param data A N by P dataframe, where N is number of observations, P is number of genes
#'
#'
#' @examples
#' data(TCGA)
#' checkData.results=checkData( data=TCGA$geneexpr )

checkData <- function(data) {
  # if (length(which(is.na(time)==TRUE))>0 |  length(which(is.na(status)==TRUE))>0 | sum(apply(data,  2,  function(x){length(which(is.na(x)==TRUE)) })) >0 ){
  if ( sum(apply(data,  2,  function(x){length(which(is.na(x)==TRUE)) })) >0 ){
    stop( "There is missing value in the data. Cannot proceed!\n" )
  } else {
    # temp2 = length(time)
    # temp3 = length(status)
    data <- as.data.frame(data)
    temp <- unlist(data)
    temp2 <- as.data.frame(temp)

    temp_outlier<-temp
    temp_outlier[which(temp>quantile(temp,0.99))]<-quantile(temp,0.99)
    temp_outlier[which(temp<quantile(temp,0.01))]<-quantile(temp,0.01)
    # temp_outlier2<-as.data.frame(temp_outlier)

    # devtools::install_github("thomasp85/patchwork")
    # library(patchwork)
    # p1<-ggplot(temp2,aes(y=temp))+geom_boxplot(fill="white",color="black",width=2)+
    #   labs(title="Original Data")+theme_classic()+theme(text=element_text(size=14,  family="serif"))+   theme(plot.title = element_text(size=18))+theme(plot.title = element_text(hjust = 0.5))
    # p2<-ggplot(temp2,aes(x=temp))+geom_histogram(fill="white",color="black",bins=30)+geom_vline(aes(xintercept=mean(temp)),color="red",linetype="dashed")+
    #   labs(title="Original Data",x="",y="Count")+theme_classic() +
    #   theme(text=element_text(size=14,  family="serif"))+   theme(plot.title = element_text(size=18))+theme(plot.title = element_text(hjust = 0.5))
    # p3<-ggplot(temp_outlier2,aes(y=temp_outlier))+geom_boxplot(fill="white",color="black",width=2)+
    #   labs(title="Without outliers")+theme_classic()+theme(text=element_text(size=14,  family="serif"))+   theme(plot.title = element_text(size=18))+theme(plot.title = element_text(hjust = 0.5))
    # p4<-ggplot(temp_outlier2,aes(x=temp_outlier))+geom_histogram(fill="white",color="black",bins=30)+geom_vline(aes(xintercept=mean(temp_outlier)),color="red",linetype="dashed")+
    #   labs(title="Without outliers",x="",y="Count")+theme_classic() +
    #   theme(text=element_text(size=14,  family="serif"))+   theme(plot.title = element_text(size=18))+theme(plot.title = element_text(hjust = 0.5))
    # p5<-ggplot(temp2,aes(y=log(temp)))+geom_boxplot(fill="white",color="black",width=2)+
    #   labs(title="Original Data")+theme_classic()+theme(text=element_text(size=14,  family="serif"))+   theme(plot.title = element_text(size=18))+theme(plot.title = element_text(hjust = 0.5))
    # p6<-ggplot(temp2,aes(x=log(temp)))+geom_histogram(fill="white",color="black",bins=30)+geom_vline(aes(xintercept=mean(log(temp))),color="red",linetype="dashed")+
    #   labs(title="Original Data",x="",y="Count")+theme_classic() +
    #   theme(text=element_text(size=14,  family="serif"))+   theme(plot.title = element_text(size=18))+theme(plot.title = element_text(hjust = 0.5))
    #
    # pA <- p1 + p2 + p3 + p4 + p5 + p6 + plot_layout(ncol = 2)
    # pA
    par(mfrow=c(3, 2), oma=c(0,0,3,0))
    boxplot(temp, main="Original Data")
    hist(temp, main="Original Data", xlab=NA, ylab=NA)

    boxplot(temp_outlier, main="Without outliers")
    hist(temp_outlier, main="Without outliers", xlab=NA, ylab=NA)
    temp_log<-temp
    temp_log[which(temp_log<0)]<-0

    boxplot(log(temp_log[which(temp_log>0)]), main="log-transformation")
    hist(log(temp_log), main="log-transformation", xlab=NA, ylab=NA)
    title("Data check", outer=TRUE)
    par(mfrow=c(1, 1))

    response <- readline(prompt="Do you want to replace outliers with 99th quantile for upper bound(1st quantile for lower bound)? [yes/no]: ")
    if(response == "y" | response == "yes"){
      cat("Outliers successfully removed", "\n")
      temp1<-apply(data,2,function(x){x[which(x>quantile(x,0.99))]<-quantile(x,0.99)
                                    return(x)})
      response <- readline(prompt="Do you want to do log-transformation? [yes/no]: ")
      if(response == "y" | response == "yes"){
        cat("Data is log-transformed", "\n")
        methods::new( "CheckedData",
                      newdata = list( data= as.data.frame(apply(temp1,2,log))),
                      inputdata = list( data = data )
        )
      } else{
        cat("Nothing changed", "\n")
        methods::new( "CheckedData",
                      newdata = list(data = as.data.frame(temp1)),
                      inputdata = list( data = data)
        )
      }
    } else{
      cat("Nothing changed", "\n")
      response <- readline(prompt="Do you want to do log-transformation? [yes/no]: ")
      if(response == "y" | response == "yes"){
        cat("Data is log-transformed", "\n")
        methods::new( "CheckedData",
                      newdata = list(data= as.data.frame(apply(data,2,log))),
                      inputdata = list( data = data )
        )
      } else{
        cat("Nothing changed", "\n")
        methods::new( "CheckedData",
                      newdata = list(data= as.data.frame(data)),
                      inputdata = list( data = data )
        )
      }
    }
  }
}



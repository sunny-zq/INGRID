% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/checkData.R
\name{checkData}
\alias{checkData}
\title{checkData function}
\usage{
checkData(data)
}
\arguments{
\item{data}{A N by P dataframe, where N is number of observations, P is number of genes}
}
\description{
provide EDA for the data
}
\examples{
data(TCGA)
checkData.results=checkData( data=TCGA$geneexpr )
}

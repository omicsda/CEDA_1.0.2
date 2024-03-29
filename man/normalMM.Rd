% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CEDA.R
\name{normalMM}
\alias{normalMM}
\title{Performing empirical Bayes modeling on limma results}
\usage{
normalMM(data, theta0)
}
\arguments{
\item{data}{A numeric matrix containing limma results and log2 gene 
expression levels that has a column nameed 'lfc' and a column 
named 'exp.level.log2'}

\item{theta0}{Standard deviation of log2 fold changes under permutations}
}
\value{
A numeric matrix containing limma results, RNA expression levels,
  posterior log2 fold ratio, log p-values, and estimates of mixture model
}
\description{
This function perform an empirical Bayes modeling on log fold ratios
and return the posterior log fold ratios.
}
\examples{
\donttest{
nmm.fit <- normalMM(data,theta0)
}

}

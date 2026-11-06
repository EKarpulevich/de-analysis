\name{Hobotnica}
\alias{Hobotnica}
\title{Hobotnica algorithm}
\usage{
Hobotnica(distMatrix, annotation)
}
\description{
Calculates Hobotnica score on dataset given grouping
}
\examples{
data(iris)
distMatrix = dist(iris[, -5])
annotation = iris[, 5]
Hobotnica(distMatrix, annotation)
}

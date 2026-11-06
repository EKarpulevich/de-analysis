\name{GenerateDistributionByPermutations}
\alias{GenerateDistributionByPermutations}
\title{Parallel generation of Hobotnica scores by permuting initial annotation of dataset.}
\usage{
GenerateDistributionByPermutations(distMatrix, annotation, nPermutations, nCores)
}
\description{
Generates Hobotnica scores by randomly permuting the annotation given matrix of distances, annotation and number of permutations to generate.
Note that by default program is executed in a non-parallel manner; set nCores = 0 to use maximal number of cores - 1 or nCores > 1 to choose number of cores on your own.
nCores = 1 corresponds to non-parallel execution.
}
\examples{
data(iris)
distMatrix = dist(iris[, -5])
annotation = iris[, 5]
distMatrix = as.matrix(distMatrix)
GenerateDistributionByPermutations(distMatrix, annotation, 300)
}

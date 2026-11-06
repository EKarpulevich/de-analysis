\name{GenerateRandomSignatures}
\alias{GenerateRandomSignatures}
\title{Random signatures generation}
\usage{
GenerateRandomSignatures(dataset, signatureLength, nSignatures, nCores)
}
\description{
Generates random signatures given dataset, length of each signature (one number, all generates signatures are of equal length) and number of signatures.
}
\examples{
data(iris)
dataset = iris[, -5]

GenerateRandomSignatures(dataset, 100, 200)


}




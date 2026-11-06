

GenerateRandomSignatures <- function(dataset, signatureLength=100, nSignatures=0, nCores=1) {
    if (length(dim(dataset)) != 2) {
        stop("The dataset dim length should be equal 2, stopping.")
    }
    if (signatureLength <= 0) {
        stop("Signature length should be greater than zero, stopping.")
    }
    if (signatureLength > dim(dataset)[1]) {
        stop ("signatureLength shoudld be less or equal to size of dataset, stopping.")
    }
    if (nSignatures == 0) {
        nSignatures = dim(dataset)[1]*10
    }


    if (nCores == 0){
        nCores = detectCores() - 1
    }

    cl <- parallel::makeForkCluster(nCores)
    doParallel::registerDoParallel(cl)

    randomGeneSignatures <- foreach (i=1:nSignatures) %dopar% {
        sample(rownames(dataset), signatureLength, replace=FALSE)
    }
    parallel::stopCluster(cl)
    return(randomGeneSignatures)



}



GenerateDistributionByPermutations <- function(distMatrix, annotation, nPermutations=0, nCores=1) {
    # 
    if (length(dim(distMatrix)) != 2) {
        stop("The distMatrix dim length should be equal 2, stopping.")
    }
    if (dim(distMatrix)[1] != dim(distMatrix)[2] ) {
        stop("distMatrix should be a square matrix, stopping")
    }
    if (length(signature) <= 0) {
        stop("Signature length should be greater than zero, stopping.")
    }

    if (length(annotation) <= 0) {
        stop("Annotation length should be greater than zero, stopping")
    }

    if (length(annotation)  != dim(distMatrix)[2]) {
        stop("Annotation length should be equal number of columns in dataset, stopping")
    }

    if (nPermutations == 0) {
        nPermutations = dim(distMatrix)[1]*5
    }
    if (nCores == 1) {
        scores <- list()
        for (i in 1:nPermutations) {
            permutedAnnotation <- sample(annotation, length(annotation), replace=FALSE)
            scores[[i]] <- Hobotnica(distMatrix, permutedAnnotation)
        }

    } else {
    if (nCores == 0) {
        nCores = detectCores() - 1
    }
    cl <- parallel::makeForkCluster(nCores)
    doParallel::registerDoParallel(cl)

    scores <- foreach (i = 1:nPermutations) %dopar% {

        permutedAnnotation <- sample(annotation, length(annotation), replace=FALSE)
        Hobotnica(distMatrix, permutedAnnotation)

    
    }
    

    }


    return (scores)


}




LengthPlotter <- function(dataset, annotation, rangedGenes,  distFunction=dist, minLength=10, maxLength=200, name=NULL,  nCores=1) {

    if ((length(rangedGenes) != dim(dataset)[1] && length(rangedGenes) < maxLength) || length(rangedGenes) <= minLength) {
        stop("lenght of rangedGenes should be equal to number of genes in dataset or equal or greater than maxLength and greater than minLength, stopping.") 
    }
    if (minLength >= maxLength) {
        stop("maxLength should be greater than minLength, stopping.")
    }
    if (nCores <= 0) {
        stop("nCores should be greater or equal to zero, stopping.")
    }
    if (name == NULL) {
        name = paste("Hobotnica",  name, paste(minLength, maxLength, sep=":"), sep=" ")
    }


    if (nCores == 1) {
    scores <- list()
    for (len in minLength:maxLength) {
        datasetCut <- dataset[1:len, ]
        distMatrix <- distFunction(datasetCut)
        scores[[len]] <- Hobotnica(dataset, annotation)

    }
    plot <- qplot(minLength:maxLength, unlist(scores), main=name) + labs(x="Signature length", y="Score")
    return (plot)

    } else {
        if (nCores == 0) {
            nCores = detectCores() -1
        } 
        scores <- foreach(i = 1:(maxLength-minLength)) %dopar% {
                datasetCut <- dataset[1:len, ]
                distMatrix <- distFunction(datasetCut)
                Hobotnica(dataset, annotation)

            }

       names(scores) <-  minLength:maxLength
       plot <- qplot(minLength:maxLength, unlist(scores), main=name) + labs(x="Signature length", y="Score")
       return (plot)
        
    }
}


MDSPlotter <- function(distMatrix, annotation, name = NULL) {
    if (dim(distMatrix)[1] != dim(distMatrix)[2]) {
        stop("distMatrix must be a square matrix, stopping.")
    }
    if (is.null(name)) {
        name <- "MDS Plot"

    }
    fit <- cmdscale(distMatrix,k=2, list.=TRUE)
    x <- fit$points[,1]
    y <- fit$points[,2]
    plot <- qplot(x, y, xlab="Coordinate 1", ylab="Coordinate 2", 
                main=name, colour=annotation)+ theme(legend.position = "none", plot.title = element_text(hjust = 0.5, size=16, face="bold"))
    
    return(plot)
}

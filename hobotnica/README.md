# Hobotnica
## Description
Hobotnica is a package that contains implementation of Hobotnica algorithm accompanied with  several  
helper functions. The algorithm measures how well the dataset was divided onto subgroups (given the  
matrix of distances between samples,or distance matrix, and a vector of labels characterizing the  
belonging of each sample in dataseet to one single group, or annotation) by a single real-valued coefficient  
which is close to null when the annotation is random and close to one when annotation is correct.

## Installation
The current recommended way is cloning this repo and executing:  
`R CMD INSTALL --no-multiarch --with-keep.source .`  
in it.

## Data 
All data and scripts related to application note in F1000Research are located in `data` branch

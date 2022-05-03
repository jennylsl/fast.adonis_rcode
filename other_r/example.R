library(vegan)
library(plyr)
library(parallel)
set.seed(10010)
source("fast.adonis.R")

# no weights
# create a distance matrix
dim<- 4
D <- dist(matrix(rnorm(10^dim*2), ncol=10))
# transform to A
A <- -0.5*as.matrix(D)^2
# dim(A)
# D<- NULL
# create a matrix for independent variables
X<- data.frame(matrix(rnorm(10^dim*2), ncol=10))
# fast adonis
fit1 <- fast.adonis(A ~ X1, data=X,permutations = 50,boot.times = 10)

# with weights
weights <- abs(rnorm(dim(A)[1]))
fit2 <- fast.adonis(A~X1, data=X,permutations = 50, boot.times = 10, weights = weights)



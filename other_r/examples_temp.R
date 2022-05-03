rm(list = ls())
library(vegan)
data(dune)
data(dune.env)
source("fa4.R")

A <- -0.5*as.matrix(vegdist(dune))^2
## default test by terms

# fit1<-adonis(dune ~ Management, data = dune.env)
# fit1
# fat.fit <- fast.adonis(A~Management, data = dune.env,
#                        boot.times = 10)
# fat.fit
# fat.fit <- fast.adonis(A~Management, data = dune.env,
#                        boot.times =0,permutations = 0)
# fat.fit
# fat.fit <- fast.adonis(A~Management, data = dune.env,
#                        boot.times =0,permutations = NULL)
# fat.fit
# fat.fit <- fast.adonis(A~Management, data = dune.env,
#                        boot.times =NULL,permutations = NULL)
# fat.fit

# adonis2(dune ~ Management+A1+Moisture, data = dune.env,by="margin")
# fat.fit <- fast.adonis(A~Management+A1+Moisture, data = dune.env,
#                        boot.times = 10, by="margin")



fat.fit <- fast.adonis(A~Management+A1+Moisture, data = dune.env,boot.times = 1,permutations = 1)
fat.fit
fat.fit2 <- fast.adonis(A~Management+A1+Moisture, data = dune.env,boot.times = 1,permutations = 2)
fat.fit2
adonis(dune ~ Management+A1+Moisture, data = dune.env)

order_list<-list(c(3,2,1),c(2,3,1))

fat.fit <- fast.adonis(A ~ Management+A1+Moisture, data = dune.env,
                       boot.times = 10, order_list = order_list,permutations = 1)
fat.fit
weights <- rep(c(1:4),5)
fat.fit <- fast.adonis(A ~ Management+A1+Moisture, data = dune.env,
                       boot.times = 5, order_list = order_list,permutations = 1,
                       weights = weights)
fat.fit
fat.fit$tab.add
adonis(dune ~ Moisture+A1+Management, data = dune.env)
adonis(dune ~ A1+Moisture+Management, data = dune.env)
dim(A)
order_list<- list(c(1,2,3),c(3,2,1))


formula <- A ~ Management+A1
permutations=999
boot.times= 10
boot.se = "WCB"
weights=NULL
boot.sample.size<- NULL
num_orders=NULL
order_list=NULL
by="terms"
data <- dune.env
parallel = getOption("mc.cores")
fff<-fast.adonis(A ~ Management+A1, data = dune.env,boot.times = 10)
fff$aov.tab


library(vegan)
library(plyr)
library(parallel)
set.seed(10010)

# no weights
# create a distance matrix
dim<- 4
D <- dist(matrix(rnorm(10^dim*2), ncol=10))

# transform to A
A <- -0.5*as.matrix(D)^2
dim(A)
# D<- NULL
# create a matrix for independent variables
X<- data.frame(matrix(rnorm(10^dim*2), ncol=10))

fit <- fast.adonis(A ~ X1, data=X,permutations = 10,boot.times = 10)
fit_ado<- adonis(D~X1, data=X,permutations = 10)
fit_ado2<- adonis2(D~X1, data=X,permutations = 10)
fit$aov.tab
weights <- abs(rnorm(dim(A)[1]))
time1<- sys.time()
fit <- fast.adonis(A~X1, data=X,permutations = 10,boot.times = 10)
fit_ado<- adonis(D~X1, data=X,permutations = 10)
fit_ado2<- adonis2(D~X1, data=X,permutations = 10)
time2 <- sys.time


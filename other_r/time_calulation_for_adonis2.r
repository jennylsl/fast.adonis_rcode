rm(list=ls())

library(vegan)
library(ape)
library(rsq)
data(dune)
data(dune.env)
################################################################################
##Generate a large dataset n = 9000 for a sampled weights w
################################################################################
N_set <- c(1000,5000,10000)
p_set <- c(1,20,50,100)
param_Set <- expand.grid(N_set,p_set)
time_length <-numeric( dim(param_Set)[1])
# param_Set_new <- cbind(param_Set,time_length)
# write.csv(param_Set_new,"mycodes.adonis2_pemut50.csv")
for(ind in 1: dim(param_Set)[1]){
  # ind<-14
  set.seed(10010+ind)
  p<- param_Set[ind,2]
  n <- param_Set[ind,1]
  
  n_set = n
  w.matrix=matrix(NA,2,20)
  w.matrix=t(apply(w.matrix, 1, function(x)x=c(sample(1:300,20))))
  w.matrix=floor(w.matrix*n_set/rowSums(w.matrix))
  X <- mvtnorm::rmvnorm(n_set,mean = rep(0,100),sigma = diag(rep(1,100)))
  X<- as.data.frame(X)
  ##build a large matrix
  w=w.matrix[1,]
  
  dune.env.large=dune.env[rep(row.names(dune.env),w),]
  dune.large=dune[rep(row.names(dune),w),]
  dune.env.large$Moisture <- as.numeric(dune.env.large$Moisture )
  d=vegdist(dune.large)
  X <- as.data.frame(X)[1:dim(dune.large)[1],]
  formu_par1 <- paste( "V",1:p," + ",collapse = "",sep="")
  formula_input = as.formula(paste0("d ~ ",substr(formu_par1,1,nchar(formu_par1)-2)))
  tim1=Sys.time()
  (adonis2(formula_input , data = X,permutations = 50, by = "terms"))
  
  tim2=Sys.time()
  time_length[ind]<- difftime(tim2, tim1, units='mins')
  print(c(ind,param_Set[ind,1],param_Set[ind,2],sum(w),time_length[ind]))
}

# use 1000 samples, v1 0.250
# use 31 s for d=10000
# use 27 min for d=30000
# use 2.219 min for d=10000
# 2.39
# d=vegdist(dune.large)
# tim1=Sys.time()
# (adonis(d ~ V1+V2+V3+V4+V5+V6+V7+V8+V9+V10 , data = X,permutations = 2))
# 
# tim2=Sys.time()
# tim2-tim1
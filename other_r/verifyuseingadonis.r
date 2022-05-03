library(vegan)
X<- pheno_inuse[,c(175,93,174)]

colnames(X)
ind_p<- 1
X_new <- X[!is.na(X[,ind_p]),]
D_new<- D[!is.na(X[,ind_p]),!is.na(X[,ind_p])]
t1<- Sys.time()
adonis(D_new ~ well_id,data=X_new, permutations = 50)
t22<- Sys.time()
colnames(pheno_inuse)
t22-t1

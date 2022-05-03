rm(list=ls())
require(data.table)
library(readxl)
library(tidyverse)
pacman::p_load(lubridate)



Weighted <-as.data.frame(fread("weightednormalized.txt"))
# colnames(UnWeighted)[1:10]
# UnWeighted_D <- UnWeighted[,-1]
Weighted_D <- Weighted[,-1]
# D<- Weighted_D
# clean the distance matrix with NA from X
load("pheno_inuse.Rdata")

result_weighted <- read.csv(file = "result_AGP_Weighted_shi.csv",header = T)
result_order <- result_weighted[,1]

rm_1elemnt.fun <- function(data){
  # data<- pheno_agp_cleaned_all
  # x<- data[,1]
  
  ind_data <- apply(data, 2, function(x){length(unique(x))})
  data_return <- data[,ind_data>1]
  return(data_return)
}
pheno_agp_nv <- pheno_inuse[,result_order]
colnames(pheno_agp_nv )
# check the correlations and find there are many highly correlated value
pick_order <- c(1:29,31:48,50:55,57:59,61,65:69,71:81,83:87,
                89,91,94,96,98:110,112:114,116,117,119:129)
# pheno_agp_nv_clean <-  rm_1elemnt.fun(pheno_agp_nv)
pheno_agp_nv_be <- pheno_agp_nv[,pick_order]
pheno_agp_nv_clean <- pheno_agp_nv_be[complete.cases(pheno_agp_nv_be),]
D_clean <- Weighted_D[complete.cases(pheno_agp_nv_be),complete.cases(pheno_agp_nv_be)]

# load("AGP.X.Rdata")
# load("agp.X.Rdata")
D<- D_clean
A <- as.matrix(-0.5*as.matrix(D)^2)
dim(A)
X<- pheno_agp_nv_clean
levels(X[,50])
rm(Weighted)
rm(Weighted_D)
rm(pheno_agp_nv)
# rm(pheno_agp_nv_be)
rm(pheno_agp_nv_clean)
rm(pheno_inuse)


tim_simu1=Sys.time()
# result_collect <- matrix(NA,dim(pheno_inuse)[2],2)
# for(ind_posi_Var in 1:dim(X)[2]){

boot.times<- 500
boot.sample.size <- dim(D)[1]
tim_simu1=Sys.time()
num_orders<- 0

if(dim(as.matrix(X))[2]==1){
  CM=0
}
if(CM == 1|CM == "all"){
  ind.covariate <- match(attr(terms(formula),"term.labels"),colnames(X))
}
set.seed(10010)
if(CM==0){
  num_orders<- 0
}
if((CM=="all")&(num_orders>0)){
  order_list_ori <- list(sample(1:p,size=p))
  for(ind_order in 1: (num_orders-1)){
    order_temp <- sample(1:p,size=p)
    order_list_ori <- c(order_list_ori,list(order_temp))
  }
}


# remove NA in data
ind.NA <-!is.na(diag(as.matrix(D))) &!(apply(as.matrix(X), 1, function(x)any(is.na(x)))) & !is.na(w)

D.new <- as.matrix(D)[ind.NA,ind.NA]
X.new <- as.data.frame(X)[ind.NA,]
w.new <- w[ind.NA]
# dim(D.new)
# dim(X.new)
# number of total variables in design matrix
num.variable <- dim(as.matrix(X.new))[2]
# first variable assigmnet for each columne
ind.col <- 1:num.variable
xnam <- colnames(as.matrix(X))
## if weight is missing or weight are equal to 1 then,
if(any(lapply((X.new),class) %in% c("charactor"))) {
  break;print("not all imputing variables are factors or numbers")}

# create dummy variables for factor variables
if( any(lapply((X.new),class) %in% "factor") & any(class(X.new) %in%c("factor"))) 
{ # update the variable assigments
  ind.col <-attr( model.matrix(~X.new),"assign")[-1]
  X.new <- model.matrix(~X.new)[,-1]}
if(any(sapply((X.new),class) %in% "factor") & any(!class(X.new) %in%c("factor"))) {
  xnam <- names(X.new)
  (fmla <- as.formula(paste("~ ", paste(xnam, collapse= "+"))))
  # update the variable assigments
  ind.col <-attr( model.matrix(fmla,data=X.new),"assign")[-1]
  X.new <- model.matrix(fmla,data=X.new)[,-1]
}
dim(X.new)
set.seed(10010)
IND.matrix <- matrix(NA, boot.times, boot.sample.size)
IND.matrix <- t(apply(IND.matrix, 1, function(x) {
  x <- sample(1:length(D.new[1,]), 
              size=boot.sample.size, 
              prob = w.new/sum(w.new), 
              replace = T)
}))
# t_001<- Sys.time()
# fre<-data.frame(table(IND.matrix[1,]))

# return a set of R^2 ,each row is one time boot, each colume referred to different design matrix

boot.R2.set <-matrix(NA,boot.times,case_when((CM == 0) ~ 1,
                                             (CM == 1) ~ 2,
                                             (CM == "all") ~ num.variable*2))
if( CM =="all"){
  R2_condi_boot_set <- matrix(NA,boot.times, num.variable-1)
}
# t_002<- Sys.time()

if(num_orders>0){
  order_list <-order_list_ori
  
  R2_add_boot_set <- matrix(NA,num_orders*boot.times, num.variable-1)
  # dim(R2_add_boot_set)
}
# 2257 diml:to original# 2.6s
for( ind.boot in 1:boot.times)
{
  # ind.boot=3
  t2<- Sys.time()
  ind.temp.ori <- IND.matrix[ind.boot,]
  # ind.temp.ori <- 1:dim(D.new)[1]
  
  ind.freq.table <- data.frame(table(ind.temp.ori))
  ind.temp <- as.numeric(as.character(ind.freq.table[,1]))
  w.new.boot <- ind.freq.table[,2]
  # t_003<- Sys.time()
  D.boot <- D.new[ind.temp,ind.temp]
  
  # check the design matrix to confirm no all the same columes
  X.boot.noncheck <- as.matrix(X.new)[ind.temp,]
  X.check.func <- function(X){
    X.boot.check <- apply(X, 2, function(vec){
      length(unique(vec))>1
    })
    return(X.boot.check)
  }
  ind_remove<-0
  if(is.null(dim(X.boot.noncheck))){
    X.boot<- X.boot.noncheck
    
  }
  if(!is.null(dim(X.boot.noncheck))){
    X.boot.temp<- X.boot.noncheck[,X.check.func(X.boot.noncheck)]
    # dim(X.boot)
    ## check the correlation to see wether it could be solve
    X.temp <- as.matrix(cbind(1,X.boot.temp))
    ind_remove<-tryCatch(solve(t(X.temp*(w.new.boot))%*%X.temp),
                         error = function(e) {
                           cor.mat.temp <- cor(X.temp[,-1])
                           diag(cor.mat.temp)<-0
                           cor.mat.temp[upper.tri(cor.mat.temp)]<-0
                           DD<- which(cor.mat.temp>.3,arr.ind = 1)
                           DD<- as.data.frame(cbind(DD,cor.mat.temp[DD]))
                           DD <- DD[order(DD$V3,decreasing = T),]
                           # ind <- DD[!duplicated(DD[,2]),2]
                           ind <- DD[!duplicated(DD[,1]),1]
                           length_ind <- length(ind)
                           order_ind <- c(1,5*(1:max(1,floor(length_ind/5))))
                           length_order_ind <- length(order_ind)
                           ind_delete<-0
                           for(ind_order in 1:length_order_ind){
                             # ind_order<-6
                             ind_delete <- ind[c(1:(order_ind[ind_order]))]
                             X.temp_try <- X.temp[,-(ind_delete+1)]
                             sucess_fun <- function(X.temp_try){
                               tt<-solve(t(X.temp_try*(w.new.boot))%*%X.temp_try)
                               if(exists("tt")){ skip_to_next <- TRUE }
                               return(skip_to_next)
                             }
                             skip_to_next<- FALSE
                             tryCatch({
                               # rm(skip_to_next)
                               skip_to_next<- sucess_fun(X.temp_try);
                               return(ind_delete)
                               # rm(ind_remove)
                               if(skip_to_next){
                                 ind_remove <- ind_delete
                                 if(exists("ind_remove")){break}
                               }
                             },
                             error = function(e){
                               
                             })
                             
                             
                           }
                         })
    # ind_delete_new <- ind_delete
    if(length(dim(ind_remove))==0){
      X.boot <- as.matrix(cbind(X.boot.temp[,-ind_remove]))
      
    }
    if(length(dim(ind_remove))>0){
      X.boot <- X.boot.temp
      
    }
   
  }
  # dim(X.boot)
  # sum(X.check.func(X.boot.noncheck))
  # which(X.check.func(X.boot.noncheck)!=1)
  # t_004<- Sys.time()
  # CM 0
  # n: number of samples
  n <-  sum(w.new.boot)
  # n.ori <- dim(D.boot.ori)[1]
  
  # compute A=-0.5*D^2
  A <- -0.5*D.boot^2
  # t_005<- Sys.time()
  # prepare a new design matrix by adding a colume 1s implying the intercept
  # dim(X.new)
  X.n <- as.matrix(cbind(1,X.boot))
  
  # X.n <- as.matrix(cbind(1,X.boot))
  # t_006<- Sys.time()
  XX <- t(X.n*(w.new.boot))%*%X.n
  # dim(XX)
  H.p1 <- tcrossprod(X.n,solve(XX))
  # dim(H.p1.new)
  # t_007<- Sys.time()
  # H.P2A <-X^tA
  H.p2A<- tcrossprod(t(X.n*w.new.boot),A)
  # t_008<- Sys.time()
  # HA= X(X^tX)^{-1}X^tA
  t.HA <- sum(H.p1*w.new.boot*t(H.p2A))
  # t_009<- Sys.time()
  # HAK =X(X^tX)^{-1}X^tAK
  # tr(AK)=tr(HAK)=tr(HKA)=1/n*tr(HKAK)
  t.AK <-sum(colSums(A*w.new.boot)*w.new.boot)
  # t_010<- Sys.time()
  # SSTO: trace(G) = -1/n*tr(AK) = -1/n*sum(A)
  # R^2 = SSR/SSTO
  R2 <- (t.HA-1/n*t.AK)/(-1/n*t.AK)
  boot.R2.set[ind.boot,1] <- R2
  # t_011<- Sys.time()
  # identify the variables or dummy variables whose have only 1 value.
  if(length(dim(ind_remove))==0){
    ind.col.new <- ind.col[-(which(X.check.func(X.boot.noncheck)!=1))][-ind_remove]
    
  }
  if(length(dim(ind_remove))>0){
    ind.col.new <- ind.col[-(which(X.check.func(X.boot.noncheck)!=1))]
    
  }
  # dim(X.boot)
  if( CM == 1 | CM == "all"){
    # aa<- as.data.frame(table(ind.col.new))
    R2.select.fun <- function(ind.col.new,
                              ind.covariate,
                              X.n,
                              H.p2A,
                              t.AK,XX,
                              w.new.boot){
      
      
      # if after booting, the target variable
      # ind.covariate<- c(105,101)
      if(any(!ind.covariate %in% ind.col.new)){
        R2.model <-NA
      } else {

        H.p2A.model <- H.p2A[c(TRUE, (ind.col.new %in% ind.covariate)),]
        
       
        model.X.n <- X.n[,c(TRUE, (ind.col.new %in% ind.covariate))]
       
        # 
        # H=X(X^tX)^{-1}X^t
        # H.p=X(X^tX)^{-1}
        # H <- X.n%*%solve(t(X.n)%*%X.n)%*%t(X.n)
        XX.model <- XX[c(TRUE, (ind.col.new %in% ind.covariate)),c(TRUE, (ind.col.new %in% ind.covariate))]
        H.p1.model<- tcrossprod(model.X.n,solve(XX.model))
        
        # HA= X(X^tX)^{-1}X^tA
        t.HA.model <- sum(H.p1.model*w.new.boot*t(H.p2A.model))
        
        R2.model <- (t.HA.model-1/n*t.AK)/(-1/n*t.AK)
      }
      
      
      return(R2.model)
    }
    # ind.covariate<- 105
    R2.model <- R2.select.fun(ind.col.new,ind.covariate,
                              X.n, H.p2A, t.AK,XX,
                              w.new.boot)
    boot.R2.set[ind.boot,2] <- R2.model
    
  }
  
  # CM: all, return single(p)and sequential regression(p-1), total 2p-1 regressions.
  if( CM =="all"){
    # all.X.boot <- X.boot[,ind.col %in% ind.covariate]
    
    # create combination of variables put in the design matrix
    num.list1 <-lapply(1:  num.variable, function(i) i)
    num.list2 <-lapply(1:  num.variable, function(i) 1:i)[-c(1,num.variable)]
    
    num.list <- c(num.list1,num.list2)
    if(num_orders>0){
      for(ind_order in 1:num_orders){
        # ind_order<-2
        num.list.temp <- lapply(1: num.variable, function(i)order_list[[ind_order]][1:i])[-c(1,num.variable)]
        
        assign(paste0("num.list",(ind_order+2)),num.list.temp)
        num.list <- c(num.list,eval(parse(text=paste0("num.list",ind_order+2))))
      }
    }
    # number of groups of variables in one time of boot (p single + p-2 sequential)
    var.total <- length(num.list)
    
    # R2 set for one single boot
    R2.single.boot <- lapply(num.list,function(ind.list) {
      R2.model <- R2.select.fun(ind.col.new, ind.list,
                                X.n, H.p2A, t.AK,XX,
                                w.new.boot)
    } )
    boot.R2.set[ind.boot,-c(1:2)] <- unlist( R2.single.boot)[1:(num.variable*2-2)]
    # then conditional R2
    R2.condi_fake <- boot.R2.set[ind.boot,c(3,c(2+num.variable+1:(num.variable-2)),1)]
    R2_condi_boot_set[ind.boot,] <- R2.condi_fake[-1]-R2.condi_fake[-length(R2.condi_fake)]
    if(num_orders>0){
      for(ind_order in 1:num_orders){
        # ind_order<-1
        R2.condi_fake_add<- numeric(num.variable)
        R2.condi_fake_add[1] <- unlist( R2.single.boot)[order_list[[ind_order]][1]]
        R2.condi_fake_add[-c(1,(num.variable))]<-  unlist( R2.single.boot)[(2*num.variable-2+1+(ind_order-1)*(num.variable-2)):(2*num.variable-2+(num.variable-2)+(ind_order-1)*(num.variable-2))]
        R2.condi_fake_add[(num.variable)]<- boot.R2.set[ind.boot,1]
        R2_add_boot_set[(ind_order-1)*boot.times+ind.boot,]<- R2.condi_fake_add[-1]-R2.condi_fake_add[-length(R2.condi_fake_add)]
      }
      
    }
  }
  t3<-Sys.time()
  
  print(c(ind.boot,t3-t2))
  # print(c(ind.boot,t_011-t_010,t_010-t_009,t_009-t_008,
  #         t_008-t_007,t_007-t_006, t_006-t_005,t_005-t_004,t_004-t_003,
  #         t_003-t2,t2-t_002,t_002-t_001,t_001-t_000,t_000-t_0000))
}
tim_simu2=Sys.time()
time_length <- difftime(tim_simu2, tim_simu1, units='mins')

ind_na_boot <- boot.R2.set[complete.cases(boot.R2.set),]
dim(ind_na_boot)
ind_na_boot_condi <- R2_condi_boot_set[complete.cases(R2_condi_boot_set),]
tim_end <- Sys.time()

## if weight is missing or weight are equal to 1 then,
# marginal
# load("AGP.X.Rdata")
result_collect_al_boot <- as.data.frame(cbind(ind_na_boot[,c(3:(2+dim(X)[2]))]))
se_result <- apply(result_collect_al_boot, 2, function(x){sd(x)})
result_collect_al_boot_se <- as.data.frame(cbind(colnames(X),se_result))
colnames(result_collect_al_boot_se)<-c("Variable","SD")
# write.csv(result_collect_al,file="result_AGP.csv")

write.csv(result_collect_al_boot_se,file="AGP/result_AGP_Weighted_boot.csv")
dim(ind_na_boot_condi)
result_collect_condi_al_boot <- as.data.frame(cbind(
                                               cbind(ind_na_boot[,3],ind_na_boot_condi)))
se_result_condi <- apply(result_collect_condi_al_boot , 2, function(x){sd(x)})
result_collect_al_boot_condi_se <- as.data.frame(cbind(colnames(X),se_result_condi))
colnames(result_collect_al_boot_condi_se)<-c("Variable","SD")
write.csv(result_collect_al_boot_condi_se,file="AGP/result_AGP_Weighted_boot_condi.csv")


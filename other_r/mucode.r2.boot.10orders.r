rm(list=ls())
library(vegan)
library(ape)
library(rsq)
library(data.table)
library(readxl)
library(plyr)
require(dplyr)
data(dune)
data(dune.env)
################################################################################
##Generate a large dataset n = 9000 for a sampled weights w
################################################################################
N_set <- c(20000)
# p_set <- c(10)
p_set <- c(1,20,50,100)
# 20000 10 : 19.171
param_Set <- expand.grid(N_set,p_set)
time_length <-numeric( dim(param_Set)[1])
# param_Set_new <- cbind(param_Set,time_length)
# write.csv(param_Set_new,"mycodes.boot10orders.csv")
# write.csv(param_Set_new,"mycodes.boot10.csv")
for(ind in 1: dim(param_Set)[1]){
  # for(ind in  20){
  # ind<-2
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
  D<-as.matrix(d)
  formula= ~ V1
  X<- X[,1:p]
  
  CM<-"all"
  w<-rep(1,dim(D)[1])
  permutations <-50
  boot.times<- 10
  boot.sample.size <- dim(as.matrix(X))[1]
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
  if(CM=="all" &(num_orders!=0)){
    order_list_ori <- list(sample(1:p,size=p))
    for(ind_order in 1: (num_orders-1)){
      order_temp <- sample(1:p,size=p)
      order_list_ori <- c(order_list_ori,list(order_temp))
    }
  }
  
  
  # remove NA in data
  ind.NA <-!is.na(diag(as.matrix(D))) &!(apply(as.matrix(X), 1, function(x)any(is.na(x)))) & !is.na(w)
  
  D.new <- as.matrix(D)[ind.NA,ind.NA]
  rm(D)
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
  # library(plyr)
  
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
    # ind.boot=1
    # t2<- Sys.time()
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
    
    if(is.null(dim(X.boot.noncheck))){
      X.boot<- X.boot.noncheck
      
    }
    if(!is.null(dim(X.boot.noncheck))){
      X.boot<- X.boot.noncheck[,X.check.func(X.boot.noncheck)]
      
    }
    # t_004<- Sys.time()
    # CM 0
    # n: number of samples
    n <-  sum(w.new.boot)
    # n.ori <- dim(D.boot.ori)[1]
    
    # compute A=-0.5*D^2
    A <- -0.5*D.boot^2
    # t_005<- Sys.time()
    # prepare a new design matrix by adding a colume 1s implying the intercept
    X.n <- as.matrix(cbind(1,X.boot))
    # X.n <- as.matrix(cbind(1,X.boot))
    # t_006<- Sys.time()
    H.p1 <- tcrossprod(X.n,solve(t(X.n*(w.new.boot))%*%X.n))
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

    if( CM == 1 | CM == "all"){
      R2.select.fun <- function(ind.col,
                                ind.covariate,
                                X.boot.noncheck,
                                H.p2A,
                                t.AK,
                                w.new.boot){
        
        model.X.boot.check <- X.check.func(X.boot.noncheck)
        
        if(all(model.X.boot.check&(ind.col %in% ind.covariate)==1)){
          R2.model <-NA
        } else {
          if(sum(model.X.boot.check-1)<0){
            ind.delet<- which(model.X.boot.check==0)
            H.p2A.model <- H.p2A[c(TRUE, (model.X.boot.check&(ind.col %in% ind.covariate))[-ind.delet]),]
          }
          
          if(sum(model.X.boot.check-1)==0){
            H.p2A.model <- H.p2A[c(TRUE, (model.X.boot.check&(ind.col %in% ind.covariate))),]
          }
          
          model.X.boot<- X.boot.noncheck[,(model.X.boot.check&(ind.col %in% ind.covariate))]
          
          # a new design matrix
          model.X.n <- as.matrix(cbind(1,model.X.boot))
          
          # H=X(X^tX)^{-1}X^t
          # H.p=X(X^tX)^{-1}
          # H <- X.n%*%solve(t(X.n)%*%X.n)%*%t(X.n)
          H.p1.model<- tcrossprod(model.X.n,solve(t(model.X.n*(w.new.boot))%*%model.X.n))
          
          # HA= X(X^tX)^{-1}X^tA
          t.HA.model <- sum(H.p1.model*w.new.boot*t(H.p2A.model))
          
          R2.model <- (t.HA.model-1/n*t.AK)/(-1/n*t.AK)
        }
        
        return(R2.model)
      }
      R2.model <- R2.select.fun(ind.col,ind.covariate, X.boot.noncheck, H.p2A, t.AK, w.new.boot)
      boot.R2.set[ind.boot,2] <- R2.model
      
    }
    
    # CM: all, return single(p)and sequential regression(p-1), total 2p-1 regressions.
    if( CM =="all"){
      # all.X.boot <- X.boot[,ind.col %in% ind.covariate]
      
      # create combination of variables put in the design matrix
      num.list1 <-lapply(1:  num.variable, function(i) i)
      num.list2 <-lapply(1:  num.variable, function(i) 1:i)[-c(1,num.variable)]
      
      num.list<- c(num.list1,num.list2)
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
        R2.model <- R2.select.fun(ind.col, ind.list, X.boot.noncheck, H.p2A, t.AK, w.new.boot)
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
    # t3<-Sys.time()
    
    # print(c(ind.boot,t2-t0,t3-t2))
    # print(c(ind.boot,t_011-t_010,t_010-t_009,t_009-t_008,
    #         t_008-t_007,t_007-t_006, t_006-t_005,t_005-t_004,t_004-t_003,
    #         t_003-t2,t2-t_002,t_002-t_001,t_001-t_000,t_000-t_0000))
  }
  tim_simu2=Sys.time()
  time_length[ind]<- difftime(tim_simu2, tim_simu1, units='mins')
  print(c(ind,param_Set[ind,1],param_Set[ind,2],sum(w),time_length[ind]))
}

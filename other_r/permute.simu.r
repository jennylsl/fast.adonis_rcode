rm(list=ls())
library(vegan)
library(ape)
library(rsq)
library(data.table)
library(readxl)
require(dplyr)
data(dune)
data(dune.env)
################################################################################
##Generate a large dataset n = 9000 for a sampled weights w
################################################################################
N_set <- c(10000)
# N_set <- c(20000)
# p_set <- c(1,20,50,100)
p_set <- c(100)
# p_set <- c(1,10,20)
# p_set <- c(50,100)
# 20000 10 : 19.171
param_Set <- expand.grid(N_set,p_set)
time_length <-numeric( dim(param_Set)[1])
# param_Set_new <- cbind(param_Set,time_length)
# write.csv(param_Set_new,"mycodes.permute10orders.csv")
# write.csv(param_Set_new,"mycodes.permute50orders.csv")
# write.csv(param_Set_new,"mycodes.permute2050orders.csv")
# write.csv(param_Set_new,"mycodes.permute50new.csv")
for(ind in 1: dim(param_Set)[1]){
  # for(ind in  20){
    # ind<-3
  set.seed(10010+ind)
  p<- param_Set[ind,2]
  # p<- 1
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
  
  tim1=Sys.time()
  
  if(dim(as.matrix(X))[2]==1){
    CM=0
  }
  if(CM == 1|CM == "all"){
    ind.covariate <- match(attr(terms(formula),"term.labels"),colnames(X))
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
  set.seed(10010)
  num_orders<- 10
  # p<-5
  if(CM==0){
    num_orders<- 0
  }
  if(CM=="all" &(num_orders>0)){
    order_list_ori <- list(sample(1:p,size=p))
    for(ind_order in 1: (num_orders-1)){
      order_temp <- sample(1:p,size=p)
      order_list_ori <- c(order_list_ori,list(order_temp))
    }
  }
  
  if(all(w==1)|missing(w)){
    R2_set <- matrix(NA,1,case_when((CM == 0) ~ 1,
                                    (CM == 1) ~ 2,
                                    (CM == "all") ~ num.variable*2))
    
    if(num_orders>0){
      order_list <-order_list_ori
     
      R2_add_set <- matrix(NA, num_orders, num.variable-1)
    }
    
    ## permutation
    D.boot <- D.new
    # dim(D.boot)
    
    # check the design matrix to confirm no all the same columes
    X.boot <- as.matrix(X.new)
    # CM 0
    # n: number of samples
    n <-  dim(D.boot)[1]
    
    # compute A=-0.5*D^2
    A <- -0.5*D.boot^2
    
    # prepare a new design matrix by adding a colume 1s implying the intercept
    X.n <- as.matrix(cbind(1,X.boot))
    # X.n <- as.matrix(cbind(1,X.boot))
    H.p1 <- tcrossprod(X.n,solve(t(X.n)%*%X.n))
    # dim(H.p1.new)
    # H.P2A <-X^tA
    H.p2A<- tcrossprod(t(X.n),A)
    
    # HA= X(X^tX)^{-1}X^tA
    t.HA <- sum(H.p1*t(H.p2A))
    # HAK =X(X^tX)^{-1}X^tAK
    # tr(AK)=tr(HAK)=tr(HKA)=1/n*tr(HKAK)
    t.AK <- sum(A)
    # SSTO: trace(G) = -1/n*tr(AK) = -1/n*sum(A)
    # R^2 = SSR/SSTO
    R2 <- (t.HA-1/n*t.AK)/(-1/n*t.AK)
    R2_set[1] <- R2
    R2.condi_true<- numeric(num.variable-1)
    if( CM == 1 | CM == "all"){
      R2.select.fun <- function(ind.col,
                                ind.covariate,
                                X.n,
                                H.p2A,
                                t.AK){
        # ind.covariate=1
        H.p2A.model <- H.p2A[c(TRUE, (ind.col %in% ind.covariate)),]
        model.X.n<- X.n[,c(TRUE,(ind.col %in% ind.covariate))]
        H.p1.model<- tcrossprod(model.X.n,solve(t(model.X.n)%*%model.X.n))
        
        # HA= X(X^tX)^{-1}X^tA
        t.HA.model <- sum(H.p1.model*t(H.p2A.model))
        
        R2.model <- (t.HA.model-1/n*t.AK)/(-1/n*t.AK)
        
        return(R2.model)
      }
      R2.model <- R2.select.fun(ind.col,ind.covariate,
                                X.n, H.p2A, t.AK)
      R2_set[2] <- R2.model
      
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
        R2.model <- R2.select.fun(ind.col, ind.list, X.n,
                                  H.p2A, t.AK)
      } )
      R2_set[-c(1:2)] <- unlist( R2.single.boot)[1:(num.variable*2-2)]
      # then conditional R2
      R2.condi_fake <- R2_set[c(3,c(2+num.variable+1:(num.variable-2)),1)]
      R2.condi_true <- R2.condi_fake[-1]-R2.condi_fake[-length(R2.condi_fake)]
      if(num_orders>0){
        for(ind_order in 1:num_orders){
          # ind_order<-1
          R2.condi_fake_add<- numeric(num.variable)
          R2.condi_fake_add[1] <- unlist( R2.single.boot)[order_list[[ind_order]][1]]
          R2.condi_fake_add[-c(1,(num.variable))]<-  unlist( R2.single.boot)[(2*num.variable-2+1+(ind_order-1)*(num.variable-2)):(2*num.variable-2+(num.variable-2)+(ind_order-1)*(num.variable-2))]
          R2.condi_fake_add[(num.variable)]<- R2_set[1]
          R2_add_set[ind_order,]<- R2.condi_fake_add[-1]-R2.condi_fake_add[-length(R2.condi_fake_add)]
        }
        
      }
    }
    
    
    if(permutations>0){
      pern<- function(permutations,n){
        p <- vegan:::getPermuteMatrix(permutations, n, strata =NULL)
        return(p)
      }
      perm_mat<-pern(permutations = permutations,1:n)
      
      permute.r2.fun <-function(permutations,perm_mat, 
                                X.n, A, t.AK, H.p1, n,
                                CM,order_list,num.list){
        # CM="all"
        R2_set_perm <- matrix(NA,permutations,case_when((CM == 0) ~ 1,
                                                        (CM == 1) ~ 2,
                                                        (CM == "all") ~ num.variable*2))
        if(!missing(order_list)|!missing(num.list)|CM=="all"|num_orders!=0){
          
          R2_add_set_perm <- matrix(NA,num_orders*permutations, num.variable-1)
        }
        ## permutation
        if( CM =="all"){
          R2_condi_set_perm <- matrix(NA,permutations,num.variable-1)
        }
        # t1<- Sys.time()
        for(ind_perm in 1:permutations){
          
          # ind_perm=1
          permut.ind <- perm_mat[ind_perm,]
          # dim(X)
          X.n.perm<- X.n[permut.ind,]
          H.p1.perm <- tcrossprod(X.n.perm,solve(t(X.n.perm)%*%X.n.perm))
          # A_perm<-A[permut.ind,permut.ind]
          
          
          # H.p2A.perm<- (t(X.n.perm)%*%A)
          
          H.p2A.perm<- tcrossprod(t(X.n.perm),A)
          
          # HA= X(X^tX)^{-1}X^tA
          t.HA.perm <- sum(H.p1.perm*t(H.p2A.perm))
          # SSTO: trace(G) = -1/n*tr(AK) = -1/n*sum(A)
          # R^2 = SSR/SSTO
          R2_perm <- (t.HA.perm-1/n*t.AK)/(-1/n*t.AK)
          R2_set_perm[ind_perm,1] <- R2_perm
          
          if( CM == 1 | CM == "all"){
            R2.select.fun <- function(ind.col,
                                      ind.covariate,
                                      X.n,
                                      H.p2A,
                                      t.AK){
              # ind.covariate=3
              H.p2A.model <- H.p2A[c(TRUE, (ind.col %in% ind.covariate)),]
              
              model.X.n<- X.n[,c(TRUE,(ind.col %in% ind.covariate))]
              
              H.p1.model<- (model.X.n%*%solve(t(model.X.n)%*%model.X.n))
              
              t.HA.model <- sum(H.p1.model*t(H.p2A.model))
              
              R2.model <- (t.HA.model-1/n*t.AK)/(-1/n*t.AK)
              return(R2.model)
            }
            R2.model_perm <- R2.select.fun(ind.col,ind.covariate,
                                           X.n.perm, H.p2A.perm, t.AK)
            R2_set_perm[ind_perm,2] <- R2.model_perm
            
          }
          
          # CM: all, return single(p)and sequential regression(p-1), total 2p-1 regressions.
          if( CM =="all"){
            
            # all.X.boot <- X.boot[,ind.col %in% ind.covariate]
            
            # create combination of variables put in the design matrix
            # num.list1 <-lapply(1:  num.variable, function(i) i)
            # num.list2 <-lapply(1:  num.variable, function(i) 1:i)[-c(1,num.variable)]
            # num.list<- c(num.list1,num.list2)
            # 
            
            
            # number of groups of variables in one time of boot (p single + p-2 sequential)
            var.total <- length(num.list)
            
            # R2 set for one single boot
            R2.single.boot_perm <- lapply(num.list,function(ind.list) {
              R2.model_perm <- R2.select.fun(ind.col, ind.list, 
                                             X.n.perm, H.p2A.perm, t.AK)
            } )
            R2_set_perm[ind_perm,-c(1:2)] <- unlist( R2.single.boot_perm)[1:(num.variable*2-2)]
            R2.condi_fake <- R2_set_perm[ind_perm,c(3,c(2+num.variable+1:(num.variable-2)),1)]
            R2_condi_set_perm[ind_perm,] <- R2.condi_fake[-1]-R2.condi_fake[-length(R2.condi_fake)]
            
            if(num_orders>0){
              for(ind_order in 1:num_orders){
                # ind_order<-2
                R2.condi_fake_add_perm<- numeric(num.variable)
                R2.condi_fake_add_perm[1] <- unlist( R2.single.boot)[order_list[[ind_order]][1]]
                R2.condi_fake_add_perm[-c(1,(num.variable))]<-  unlist( R2.single.boot)[(2*num.variable-2+1+(ind_order-1)*(num.variable-2)):(2*num.variable-2+(num.variable-2)+(ind_order-1)*(num.variable-2))]
                R2.condi_fake_add_perm[(num.variable)]<- R2_set[1]
                R2_add_set_perm[(ind_order-1)*permutations+ind_perm,]<- R2.condi_fake_add[-1]-R2.condi_fake_add[-length(R2.condi_fake_add)]
              }
              
            }
          }
        } 
        if(CM==0|CM==1){
          perm_return <- R2_set_perm
        }
        if(CM=="all"){
          perm_return <- list(R2_set_perm,R2_condi_set_perm)
          temp <- rowSums(t(R2_set_perm)-as.vector(R2_set)>0)/permutations
          temp_condi <- rowSums(t(R2_condi_set_perm)-as.vector(R2.condi_true)>0)/permutations
          
        }
        if(CM=="all" & !(missing(order_list))){
          perm_return <- list(R2_set_perm,R2_condi_set_perm,R2_add_set_perm)
        }
        return( perm_return )
        
      }
      
      perm_result <- permute.r2.fun(permutations,perm_mat, 
                                    X.n, A, t.AK, H.p1, n,
                                    CM=CM,order_list,num.list)
      
      
    }
  }
  tim2=Sys.time()
  time_length[ind]<- difftime(tim2, tim1, units='mins')
  print(c(ind,param_Set[ind,1],param_Set[ind,2],sum(w),time_length[ind]))
}

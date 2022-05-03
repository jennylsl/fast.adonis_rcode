library(data.table)
library(vegan)
library(readxl)
require(dplyr)
# install.packages("readxl")
# rm(list=ls())
############################################################################################
# D.AARP <- fread("aarp/AARPdistMatBrayCurtis.txt")[,-1]
# D.AHS <- fread("ahs/AHSdistMatBrayCurtis.txt")[,-1]
D.PLCO <- fread("plco/PLCOdistMatBrayCurtis.txt")[,-1]


# load phenotype for AARP, AHS, PLCO
# phenoAARP <- read.csv("aarp/pheno_alphaAARP.csv",header = 1)
# phenoAHS <- read.csv("ahs/pheno_alphaAHS.csv",header = 1)
phenoPLCO <- read.csv("plco/pheno_alphaPLCO.csv",header = 1)
pcaPLCO <- read_excel("plco/PLCO.alpha.diversity.xlsx",na="--")
colnames(pcaPLCO)<-c("id","osa","pd","shannon")
# clean data with order sample ID,

# phenoAHS <- phenoAHS[match(colnames(D.AHS),phenoAHS$Sample.ID),]
# phenoAARP <- phenoAARP[match(colnames(D.AARP),phenoAARP$Sample.ID),]
phenoPLCO <- phenoPLCO[match(colnames(D.PLCO),phenoPLCO$Sample.ID),]

## cluster built from weights
#

Compute.var.R2.Distance.Matrix.Boot.new <- function(D=D.PLCO, 
                                                    X=X, 
                                                    w=phenoPLCO$weight, 
                                                    formula= ~ race7_new , 
                                                    boot.times=10,
                                                    boot.sample.size=10000,
                                                    CM=0){
  D<-as.matrix(d)
  boot.times=10
  formula= ~ V1
  X<- X[,1:10]

  CM<-"all"
  w<-rep(1,dim(D)[1])
  boot.sample.size=dim(D)[1]
  permutations <-50

  ## 
  # t0<-Sys.time()
  # X <- phenoPLCO[,c("bmi_curr_new","race7_new","smk_new","educat_new",
  #                   "age_at_collection","drinker_status_dhq_new",
  #                   "sex","cig_years","diabetes_f","y0_anycancer")]
  # # X <-pcaPLCO
  # X$smk_new <- factor(X$smk_new)
  # 
  # w<-phenoPLCO$weight
  # CM="all"
  # boot.sample.size = 10000
  # boot.times = 10
  # CM="all"
  # when there is one variable in design matrix X, automaticly set CM=0
  # if(dim(X)[2] == 1){
  #   CM=0
  # }

  if(CM == 1|CM == "all"){
    ind.covariate <- match(attr(terms(formula),"term.labels"),colnames(X))
  }
  if(dim(as.matrix(X))[2]==1){
    CM=0
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
  
  if(all(w==1)|missing(w)){
    R2_set <- matrix(NA,1,case_when((CM == 0) ~ 1,
                                    (CM == 1) ~ 2,
                                    (CM == "all") ~ num.variable*2))
    num_orders<- 10
    
    if(num_orders>0){
      order_list <- list(sample(1:10,size=10),
                         sample(1:10,size=10),
                         sample(1:10,size=10),
                         sample(1:10,size=10),
                         sample(1:10,size=10),
                         sample(1:10,size=10),
                         sample(1:10,size=10),
                         sample(1:10,size=10),
                         sample(1:10,size=10),
                         sample(1:10,size=10))
      R2_add_set <- matrix(NA,num_orders, num.variable-1)
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
          R2.condi_fake_add[-c(1,(num.variable))]<-  unlist( R2.single.boot)[(2+num.variable+1+(ind_order-1)*(num.variable-2)):(2+num.variable+8+(ind_order-1)*(num.variable-2))]
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
                                X.n, A, t.AK, H.p1, n,CM,order_list){
        # CM="all"
        R2_set_perm <- matrix(NA,permutations,case_when((CM == 0) ~ 1,
                                                        (CM == 1) ~ 2,
                                                        (CM == "all") ~ num.variable*2))
         if(!missing(order_list)){
          
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
                R2.condi_fake_add_perm[-c(1,(num.variable))]<-  unlist( R2.single.boot)[(2+num.variable+1+(ind_order-1)*(num.variable-2)):(2+num.variable+8+(ind_order-1)*(num.variable-2))]
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
        }
        if(CM=="all" & !(missing(order_list))){
          perm_return <- list(R2_set_perm,R2_condi_set_perm,R2_add_set_perm)
        }
        return( perm_return )

      }
      
      perm_result <- permute.r2.fun(permutations,perm_mat, 
                                    X.n, A, t.AK, H.p1, n,
                                    CM="all",order_list)

      
    }
  }

  # difftime(t1,t0,units = "mins")
  # prepare X design matrix with categorical variables
  # screen variables that are not classified as "factor" and "numeric"
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
  ########## prepare boot matrix #############
  # prepare a matrix, boot.times * boot.samples: each row is a set of indicator of IDs
  # boot.times=2
  # boot.sample.size=2600
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
  library(plyr)
  
  # return a set of R^2 ,each row is one time boot, each colume referred to different design matrix

  boot.R2.set <-matrix(NA,boot.times,case_when((CM == 0) ~ 1,
                                               (CM == 1) ~ 2,
                                               (CM == "all") ~ num.variable*2))
  
  # t_002<- Sys.time()
  ############# Next to calculate R2 for different modes within one boot############
  
  # 2257 diml:to original# 2.6s
  for( ind.boot in 1:boot.times)
  {
    # ind.boot=1
    t2<- Sys.time()
    # ind.temp.ori <- IND.matrix[ind.boot,]
    ind.temp.ori <- 1:dim(D.new)[1]
    
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
      
      # number of groups of variables in one time of boot (p single + p-2 sequential)
      var.total <- length(num.list)
      
      # R2 set for one single boot
      R2.single.boot <- lapply(num.list,function(ind.list) {
        R2.model <- R2.select.fun(ind.col, ind.list, X.boot.noncheck, H.p2A, t.AK, w.new.boot)
      } )
      boot.R2.set[ind.boot,-c(1:2)] <- unlist( R2.single.boot)
    }
    # t3<-Sys.time()
    
    # print(c(ind.boot,t2-t0,t3-t2))
    # print(c(ind.boot,t_011-t_010,t_010-t_009,t_009-t_008,
    #         t_008-t_007,t_007-t_006, t_006-t_005,t_005-t_004,t_004-t_003,
    #         t_003-t2,t2-t_002,t_002-t_001,t_001-t_000,t_000-t_0000))
  }
  
  # t2.b2<-Sys.time()
  # 
  # t3<-Sys.time()
  # t3-t2.b2
  
  ############## summary the result ################
  result.table1 <- data.frame("WS Mean R2"=colMeans(boot.R2.set),
                              "WS SE"=apply(boot.R2.set, 2, sd),
                              "WS CI"=apply(boot.R2.set, 2,function(colume)paste(round(quantile(colume,probs=c(0.05,0.95)),5),collapse = "~")))
  
  ########### output#############
  
  # CM=0 build default model formula D~X
  model.default <- as.formula(
    paste(deparse(substitute(D)), formula[[2]], sep=" ~ ")
  )
  
  # output matrix
  output.list <- list(c(model.default,result.table1[1,]))
  
  # CM=1
  if(CM==1|CM=="all"){
    
    model.specified <- formula
    list.formula <- c(model.default,formula)
    
    output.list <- lapply((1:length(list.formula)),function(l){
      c(list.formula[[l]],result.table1[l,])
    })
  }
  
  # CM=all
  if(CM=="all"){
    
    model.sequence <- (lapply(num.list,function(x) {
      as.formula(paste(" ",paste(xnam[x], collapse = " + "), sep = " ~ "))
    })  )
    list.formula.all <-c(list.formula, model.sequence)
    
    output.list <- lapply((1:length(list.formula.all)),function(l){
      c(list.formula.all[[l]],result.table1[l,])
    })
  }
  return(output.list)
}
D<-as.matrix(d)
p<- dim(D)[1]
boot.times=10
formula= ~ A1
X<- X
CM<-"all"
w<-rep(1,dim(D)[1])
boot.sample.size=dim(D)[1]
t0<-Sys.time()
t4<- Sys.time()
Compute.var.R2.Distance.Matrix.Boot.new(D =as.matrix(d),
                                        X= X,
                                        w=rep(1,p),
                                        formula= ~ V1,
                                        boot.times = 3,
                                        boot.sample.size = p,
                                        CM="all")
t5<- Sys.time()
t5-t4
# use 6 s for d=10000
# use 2.4 min for d=30000
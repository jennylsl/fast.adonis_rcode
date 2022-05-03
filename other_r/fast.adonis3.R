
#### 
## `fast.adonis` is a R function that is used for analysis of variance using distance matrices. It can be used
## for both weighted and unweighted data. For unweighted data, it performs similarly with R functions 'adonis' 
## and 'adonis2' from package(vegan). Besides a permutation test, `fast.adonis` provides a standard error for R2.
## The standard error can by computed by within cluster bootstraps or bootstrapping samples by weights. For unweighted data,
## `fast.adonis` does not provide permutation tests.
####
# formula : Model Formula. The LHS must be a matrix that is computed from a dissimilarity matrix.
#           For example, given a Bray Curtis dissimilarity matrix D, LHS= -0.5*D^2. The RHS defines the
#           independent variables. These can be continuous variables or factors.
# data    : the data frame for the independent variables.
# permutations : number of permutations
# boot.times : times for bootstrapping
# boot.se : the name of any method for bootstraps. Bootstraps are used to calculate standard errors for R2.
#           'boot.se = WCB' is within cluster bootstrap. Clusters are first created by weights.'boot.se = AS' is 
#           booting samples directly by weights. By default, 'boot.se = WCB'.
# boot.sample.size :  sample size for bootstrap
# weights : If weights is NULL or a vector of 1s, `fast.adonis` conducts both permutations and bootstraps.
#           Otherwise, fast.adonis` does not conduct permutations
# order_list: a list of different orders of input variables. For example, an input formula is RHS ~ A+B+C. 
#             The order_list can be list(c(3,2,1),c(2,3,1)). The output will be three tabs. 
#             The first is the analysis table for formula, RHS ~ A+B+C. The second is the table for formula, RHS ~ C+B+A.
#             The last is the table for formula, RHS ~ B+C+A.
# by      : 'by=terms' will assess significance for each term (sequentially from first to last). 'by=margin' will assess
#           the marginal effects of the terms (each marginal term analysis in a model without all other variables). 
#           If 'by=margin' is set, "order_list" will be set to NULL.
# parallel : number of parallel processes. With parallel = 1 uses ordinary, non-parallel processing. The parallel processing
#           is done with the parallel package.
####
#  Reference: Fitting multivariate modesl to community data: A comment on distance-based redundancy analysis.Ecology 2001
#             SSR:trace(H hat(Y)hat(Y') H)=trace(HGH)=trace(HA)-1/n*trace(HAK)-1/n*trace(HKA)+1/n^2*trace(HKAK)
#             G: Gower's centered matrix by centering the elements of A; A: A=-0.5*D^2; K: all elements are 1, n.samples * n.samples;
#             H: idempotent "hat" matrix H=X(X'X)^(-1)X'.
# F test:     old version: (sum( G * t(H) )/df.Exp)/(sum( G * t(I-H.snterm) )/df.Res). df.Exp and df.Res are degrees of freedom.
#             The new F test: ((trace(HA)-1/n*trace(KA))/df.Exp)/((-trace(KA))/df.Res)
`fast.adonis`<-
  function(formula, data=NULL, permutations=999, boot.times= 100, boot.se = "WCB",
           boot.sample.size= NULL, weights= NULL, order_list=NULL,
           by="terms", parallel = getOption("mc.cores")
           )
  {
    # clean parameters
    # ## we accept only by = "terms" or "margin" 
    A<- -0.5*as.matrix(vegdist(dune))^2
    formula <- A~Management+A1+Moisture

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
    # order_list<-list(c(3,2,1),c(2,3,1))
    # 
    # 
    set.seed(10010)
    source("suppli.code.R")
    
    if ( !is.null(by) ){
      by <- match.arg(by, c("terms", "margin"))}
    if ( !is.null( boot.se) ){
      boot.se <- match.arg( boot.se, c("WCB", "AS"))}
    
    ## Set first parallel processing for all terms
    if ( is.null(parallel) ){
      parallel <- 1}
    isParal <-  parallel > 1
    isMulticore <- .Platform$OS.type == "unix"
    
    ## evaluate data
    Terms <- terms(formula, data = data)
    lhs <- formula[[2]]
    lhs <- eval(lhs, data, parent.frame()) 
    formula[[2]] <- NULL                # to force evaluation
    rhs.frame <- model.frame(formula, data, drop.unused.levels = TRUE) # to get the data frame of rhs
    rhs <- model.matrix(formula, rhs.frame) # and finally the model.matrix
    grps <- attr(rhs,"assign")
    qrhs <- qr(rhs)
    options(contrasts=options()$contrasts)
    
    ## Take care of aliased variables and pivoting in rhs
    rhs <- rhs[, qrhs$pivot, drop=FALSE]
    rhs <- rhs[, 1:qrhs$rank, drop=FALSE]
    grps <- grps[qrhs$pivot][1:qrhs$rank]
    ind.col <- grps[-1]
    u.grps <- unique(grps)
    nterms <- length(u.grps) - 1
    
    # clean order_list and num_orders
    if(!is.null(order_list) & by%in% "terms" & nterms > 1){
      num_orders <- length(order_list)}else{
        num_orders <- NULL
        order_list <- NULL
      }
    if ( nterms < 1 ){
      stop("right-hand-side of formula has no usable terms")}
    
    n <- nrow(lhs) # sample size
    if ( is.null(boot.sample.size) ){
      boot.sample.size <- n}
    ## generate boot samples if weights is not null
    if ( is.null(weights) ){
      weights <- rep(1,n)}
    # degree of freedom
    df.Exp <- sapply(u.grps[-1], function(i) sum(grps==i) )
    df.Res <- n - qrhs$rank
    
    # compute the unbiased R2 for whole population by weights
    ##  Sequential computing: create combination of variables put in the design matrix
    # add a function num.list_generate 
    num.list <- num.list_gener(nterms, num_orders, order_list)

    ## t.AK : trace of matrix calculated by A%*%K. A is the lhs of the formula.
    ##      K is a n x n matrix that elements are 1s. The output of t.AK is unchanged for permutation.
    t.AK <- sum(colSums(lhs*weights)*weights)
    
    ## compute R2
    R2.original <- R2.calc(lhs, rhs, weights, nterms, ind.col,
                           t.AK, num.list, num_orders, SS=TRUE)
    temp.dim <- ifelse(is.null(num_orders)|is.null(order_list),nterms*2,c(nterms*2+num_orders*nterms))
    
    ## bootstrap to see its standard derivation
    if ( !is.null( boot.times) ){
      Ind.matrix <- matrix(NA, boot.times, n)
      
      if ( boot.se %in% c("WCB") ){
        Ind.matrix <- WCB.sample.fun(Ind.matrix, weights, n)}
      if ( boot.se %in% c("AS") ){
        Ind.matrix <- AS.sample.fun(Ind.matrix, weights, n)}
      boot.R2.set <- matrix(NA, boot.times, temp.dim )
      
      # parallel computing
      if ( isParal && isMulticore ){
        boot.R2.set <- t(matrix(unlist(mclapply(1:boot.times,function(ind_boot){
          boot.fun(ind_boot, weights, rhs, lhs, ind.col, Ind.matrix,
                   nterms, num.list, num_orders)},mc.cores = parallel)), nrow=temp.dim, ncol=boot.times))
      }else {
        for( ind.boot in 1:boot.times ){
          # ind.boot<-1
          boot.R2.set[ind.boot,] <- boot.fun(ind.boot, weights, rhs, lhs,ind.col, Ind.matrix, nterms,
                                             num.list, num_orders)}
      }
      # SE
      SE.Mat <- apply(boot.R2.set, 2, sd)
    }else { # no permutations
      SE.Mat <- rep(NA, temp.dim)
    }
    # Permutations
    ## require R package Vegan
    require(vegan)
    ## we do permutations only when samples do not require weights
    
    if ( permutations>0 & all(weights==1) ){
      # generate permutation indicators
      # set.seed(10010)
      perm_mat <- pern(permutations = permutations,1:n)
      R2_set_perm <- matrix(NA,permutations,temp.dim+2)        
      weights.perm <- rep(1,n)
      
      ## permutation
      if ( isParal && isMulticore ){
        R2_set_perm <- t(matrix(unlist(mclapply(1:permutations,function(ind_perm){
          permut.ind <- perm_mat[ind_perm,]
          lhs.perm<- lhs[permut.ind,permut.ind]

          R2.calc(lhs.perm, rhs, weights.perm, nterms, ind.col, t.AK,
                  num.list, num_orders, SS=TRUE)[[2]]},
          mc.cores = parallel)),nrow=temp.dim+2,ncol=permutations))
      }else{ # create R2 perm table
        for( ind_perm in 1:permutations ){
          # ind_perm=2
          permut.ind <- perm_mat[ind_perm,]
          lhs.perm<- lhs[permut.ind,permut.ind]
          R2_set_perm[ind_perm,] <- R2.calc(lhs.perm, rhs, weights.perm, nterms,
                                            ind.col, t.AK, num.list, num_orders, SS=TRUE)[[2]]
        }

        df.Exp.l<- c(sum(df.Exp), df.Exp[1], df.Exp[-nterms], df.Exp[-1])
        if(!is.null(num_orders)){
          
          for(ind_orders in 1:num_orders){
            # ind_orders<-1
            orders <- order_list[[ind_orders]]
            df.Exp.l<- c(df.Exp.l,df.Exp[orders])
          }
        }
        # F test 
        F.Mod.perm <- t(t(R2_set_perm[,1:temp.dim])/(df.Exp.l))/(R2_set_perm[,c(temp.dim+1)]/df.Res)
        
      }
    } else { # no permutations
      F.Mod.perm <- rep(NA, temp.dim)
    }
    
    if ( nterms==1 ){
      # indicator for positions for selected variables
      ind_posi <- c(1,(temp.dim+1),(temp.dim+2))
    }
    if ( by=="terms" & nterms>=2 ){
      ind_posi <- c(2,(1+nterms+1:(nterms-1)),(temp.dim+1),(temp.dim+2))
    }
    if ( by=="margin" & nterms>=2 ){
      ind_posi <- c(3:(nterms+1),2*nterms,(temp.dim+1),(temp.dim+2))
    }

    SumsOfSqs <- R2.original[[2]][ind_posi]
    F.Mod <- SumsOfSqs[-c(length(SumsOfSqs)-1,length(SumsOfSqs))]/(df.Exp)/(SumsOfSqs[c(length(SumsOfSqs)-1)]/df.Res)
    P <- (rowSums(t(F.Mod.perm[,ind_posi[-c(length(ind_posi)-1,length(ind_posi))]])>=as.vector(F.Mod))+1)/(permutations+1)
    SE <- SE.Mat[ind_posi[-c(length(ind_posi)-1,length(ind_posi))]]
    
    # output
    tab <- data.frame(Df = c(df.Exp, df.Res, n-1),
                      SumsOfSqs = SumsOfSqs,
                      R2 = SumsOfSqs/SumsOfSqs[length(SumsOfSqs)],
                      se.R2 = c(SE, NA, NA),
                      F.Model = c(F.Mod, NA,NA),
                      P = c(P, NA, NA))
    
    rownames(tab) <- c(attr(attr(rhs.frame, "terms"), "term.labels")[u.grps],
                       "Residuals", "Total")
    colnames(tab)[ncol(tab)] <- "Pr(>F)"
    if ( by=="terms" & all(weights==1) ){
      attr(tab, "heading") <- c(vegan:::howHead(attr(perm_mat, "control")),
                                "Terms added sequentially (first to last)\n")
    }
    if ( by=="margin" & all(weights==1) ){
      attr(tab, "heading") <- c(vegan:::howHead(attr(perm_mat, "control")),
                                "Marginal effects of terms\n")
    }
    if ( any(weights!=1)){
      attr(tab, "heading") <- c("Weights included (no permutation)")
    }
    
    class(tab) <- c("anova", class(tab))
    out <- list(aov.tab = tab,
                call = match.call(),
                model.matrix = rhs,
                terms = Terms)
    
    out$call$formula <- Terms
    
    # if additional orders are input
    if((!is.null(num_orders)) & (!is.null(order_list))){
      tab1 <- list(tab)
      for(ind_orders_tab in 1:num_orders){
        # ind_orders_tab<-1
        orders <- order_list[[ind_orders_tab]]
        ind_posi_ord <- c((2*nterms+1:(nterms)+(ind_orders_tab-1)*(nterms)),temp.dim+1,temp.dim+2)
        
        SumsOfSqs.ad <- R2.original[[2]][ind_posi_ord]
        df.Exp.ad <- sapply(u.grps[-1], function(i) sum(grps==i) )[orders]
        F.Mod.ad <- SumsOfSqs.ad[-c(length(SumsOfSqs.ad)-1,length(SumsOfSqs.ad))]/(df.Exp.ad)/(SumsOfSqs.ad[c(length(SumsOfSqs.ad)-1)]/df.Res)
        P.ad <- (rowSums(t(F.Mod.perm[,ind_posi_ord[-c((length(ind_posi_ord)-1),length(ind_posi_ord))]])>=as.vector(F.Mod.ad))+1)/(permutations+1)
        SE.ad <- SE.Mat[ind_posi_ord[-c((length(ind_posi_ord)-1),length(ind_posi_ord))]]
        # 
        tab.ad <- data.frame(Df = c(df.Exp.ad, df.Res, n-1),
                             SumsOfSqs = SumsOfSqs.ad,
                             R2 = SumsOfSqs.ad/SumsOfSqs.ad[length(SumsOfSqs.ad)],
                             se.R2 = c(SE.ad, NA, NA),
                             F.Model = c(F.Mod.ad, NA,NA),
                             P = c(P.ad, NA, NA))
        rownames(tab.ad) <- c(attr(attr(rhs.frame, "terms"), "term.labels")[u.grps][orders],
                              "Residuals", "Total")
        colnames(tab.ad)[ncol(tab.ad)] <- "Pr(>F)"
        attr(tab.ad, "heading") <- c(vegan:::howHead(attr(perm_mat, "control")),
                                  "Terms added sequentially (first to last)\n")
        class(tab.ad) <- c("anova", class(tab.ad))

        tab1 <- c(tab1,list(tab.ad))
      }
      tab1<- tab1[-1]
      
      out <- list(aov.tab = tab, call = match.call(), model.matrix = rhs,
                  terms = Terms, tab.add=tab1)
    }
    
    class(out) <- "adonis"
    out
  }


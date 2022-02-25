<<<<<<< HEAD
## fast.adonis is a R function that is used for analysis of variance using distance matrices
# 
`fast.adonis`<-
  function(formula, data=NULL, permutations=999, boot.times= 100, boot.se = "WCB",
           boot.sample.size= 2000, weights= NULL, num_orders=NULL, order_list=NULL,
           by="terms",...
=======
`fast.adonis` <-
  function(formula, data=NULL, permutations=999, boot.times= 100, boot.se = "WCB",
           boot.sample.size= 2000, weights= NULL, num_orders=NULL, order_list=NULL,
           by="terms"
>>>>>>> ef20608ebb44871f24734d4576a97eb6b300bf6d
           )
  {
    # clean parameters
    ## we accept only by = "terms", "margin" or NULL
    set.seed(10010)
    
    if (!is.null(by)){
      by <- match.arg(by, c("terms", "margin"))}
    if(is.null(boot.sample.size)){
      boot.sample.size <- n}
    if(!is.null( boot.se)){
      boot.se <- match.arg( boot.se, c("WCB", "AS"))}
    
    ## evaluate data
    Terms <- terms(formula, data = data)
    lhs <- formula[[2]]
    lhs <- eval(lhs, data, parent.frame()) 
    formula[[2]] <- NULL                # to force evaluation
    rhs.frame <- model.frame(formula, data, drop.unused.levels = TRUE) # to get the data frame of rhs
    rhs <- model.matrix(formula, rhs.frame) # and finally the model.matrix
    grps <- attr(rhs,"assign")
    qrhs <- qr(rhs)
    
    ## Take care of aliased variables and pivoting in rhs
    rhs <- rhs[, qrhs$pivot, drop=FALSE]
    rhs <- rhs[, 1:qrhs$rank, drop=FALSE]
    grps <- grps[qrhs$pivot][1:qrhs$rank]
    ind.col <- grps[-1]
    u.grps <- unique(grps)
    nterms <- length(u.grps) - 1
    if (nterms < 1){
      stop("right-hand-side of formula has no usable terms")}

    n <- nrow(lhs) # sample size
    
    ## generate boot samples if weights is not null
    if(is.null(weights)){
      weights <- rep(1,n)}

    # compute the unbiased R2 for whole population by weights
    ##  Sequential computing: create combination of variables put in the design matrix
    # add a funtion num.list_generate 
    num.list <- num.list_gener(nterms, num_orders, order_list)

    ## t.AK : trace of matrix calculated by A%*%K. A is the lhs of the formula.
    ##      K is a n x n matrix that elements are 1s. The output of t.AK is unchanged for permutation.
    t.AK <- sum(colSums(lhs*weights)*weights)
    
    ## compute R2
    R2.original <- R2.calc(lhs, rhs, weights, nterms, t.AK, num.list, num_orders, SS=TRUE)
    
    ## bootstrap to see its standard derivation
    if(!is.null( boot.times)){
      Ind.matrix <- matrix(NA, boot.times, n)
      if(boot.se %in% c("WCB")){
        Ind.matrix <- WCB.sample.fun(Ind.matrix, weights, n)}
      if(boot.se %in% c("AS")){
        Ind.matrix <- AS.sample.fun(Ind.matrix, weights, n)}
      boot.R2.set <- matrix(NA, boot.times, nterms*2)
      for( ind.boot in 1:boot.times){
        boot.R2.set[ind.boot,] <- boot.fun(ind.boot, weights, rhs, lhs)}
      SD.Mat <- apply(boot.R2.set, 2, sd)
      }
    
    # Permutations
    ## require R package Vegan
    require(vegan)
    ## we do permutations only when samples do not require weights
<<<<<<< HEAD
    
    if(permutations>0 & all(weights==1)){
      # generate permutation indicators
      perm_mat <- pern(permutations = permutations,1:n)
      
      # create R2 perm table
      R2_set_perm <- matrix(NA,permutations,nterms*2)
      
=======
    if(permutations>0 & all(weights==1)){
      # generate permutation indicators
      perm_mat <- pern(permutations = permutations,1:n)
      # create R2 perm table
      R2_set_perm <- matrix(NA,permutations,nterms*2)
>>>>>>> ef20608ebb44871f24734d4576a97eb6b300bf6d
      # order list for 
      if(!missing(order_list)&!is.null(order_list)){
        R2_add_set_perm <- matrix(NA,num_orders*permutations, nterms-1)}
      ## permutation
      for(ind_perm in 1:permutations){
        # ind_perm=1
        permut.ind <- perm_mat[ind_perm,]
        
        lhs.perm<- lhs[permut.ind,permut.ind]
        weights.perm <- rep(1,n)
        # 
        R2_set_perm[ind_perm,] <- R2.calc(lhs.perm, rhs, weights.perm, nterms, t.AK, num.list, num_orders)
        
      }
      perm_return <- list(R2_set_perm)
      # p values
      
      P <- ((rowSums(t(R2_set_perm)-as.vector(R2.original[[1]])>0)+1)/(permutations+1))
      # (rowSums(t(f.perms) >= F.Mod - EPS)+1)/(permutations+1)

<<<<<<< HEAD
    } else { # no permutations
       P <- rep(NA, nterms)
    }
    if(nterms==1 & (is.null(num_orders)) ){
=======
    }
    else { # no permutations
       P <- rep(NA, nterms)
    }
    if(nterms==1 & (is.null(num_orders)) ){
      
>>>>>>> ef20608ebb44871f24734d4576a97eb6b300bf6d
      SumsOfSqs <- R2.original[[2]]
      df.Exp <- 1
      df.Res <- n - 1
      F.Mod <- SumsOfSqs[-c(length(SumsOfSqs)-1,length(SumsOfSqs))]/(df.Exp)/(SumsOfSqs[c(length(SumsOfSqs)-1)]/df.Res)
<<<<<<< HEAD
=======
      R2.t <- R2.original[[1]]
>>>>>>> ef20608ebb44871f24734d4576a97eb6b300bf6d
      tab <- data.frame(Df = c(df.Exp, df.Res, n-1),
                        SumsOfSqs = SumsOfSqs,
                        F.Model = c(F.Mod, NA,NA),
                        R2 = SumsOfSqs/SumsOfSqs[length(SumsOfSqs)],
<<<<<<< HEAD
                        se.R2 = c(SD.Mat,NA, NA),
=======
>>>>>>> ef20608ebb44871f24734d4576a97eb6b300bf6d
                        P = c(P, NA, NA))
    }
    if(by=="terms" & nterms>=2 & (is.null(num_orders)) ){
      SumsOfSqs <- R2.original[[2]][c(2,(1+nterms+1:(nterms+1)))]
      df.Exp <- sapply(u.grps[-1], function(i) sum(grps==i) )
      df.Res <- n - qrhs$rank
      F.Mod <- SumsOfSqs[-c(length(SumsOfSqs)-1,length(SumsOfSqs))]/(df.Exp)/(SumsOfSqs[c(length(SumsOfSqs)-1)]/df.Res)
<<<<<<< HEAD
=======
      R2.t <- R2.original[[1]][c(2,(1+nterms+1:(nterms-1)))]
>>>>>>> ef20608ebb44871f24734d4576a97eb6b300bf6d
      tab <- data.frame(Df = c(df.Exp, df.Res, n-1),
                        SumsOfSqs = SumsOfSqs,
                        F.Model = c(F.Mod, NA,NA),
                        R2 = SumsOfSqs/SumsOfSqs[length(SumsOfSqs)],
<<<<<<< HEAD
                        se.R2 = c(SD.Mat[c(2,(1+nterms+1:(nterms-1)))],NA, NA),
=======
>>>>>>> ef20608ebb44871f24734d4576a97eb6b300bf6d
                        P = c(P[c(2,(1+nterms+1:(nterms-1)))], NA, NA))

    }
    if(by=="margin" & nterms>=2 & (is.null(num_orders)) ){
<<<<<<< HEAD
      SumsOfSqs <- R2.original[[2]][c(2:(1+nterms),(length(R2.original[[2]])-1),length(R2.original[[2]]))]
      df.Exp <- sapply(u.grps[-1], function(i) sum(grps==i) )
      df.Res <- n - qrhs$rank
      F.Mod <- SumsOfSqs[-c(length(SumsOfSqs)-1,length(SumsOfSqs))]/(df.Exp)/(SumsOfSqs[c(length(SumsOfSqs)-1)]/df.Res)
=======
      SumsOfSqs <- R2.original[[2]][c(2:(1+nterms))]
      df.Exp <- sapply(u.grps[-1], function(i) sum(grps==i) )
      df.Res <- n - qrhs$rank
      F.Mod <- SumsOfSqs[-c(length(SumsOfSqs)-1,length(SumsOfSqs))]/(df.Exp)/(SumsOfSqs[c(length(SumsOfSqs)-1)]/df.Res)
      R2.t <- R2.original[[1]][c(2:(1+nterms))]
>>>>>>> ef20608ebb44871f24734d4576a97eb6b300bf6d
      tab <- data.frame(Df = c(df.Exp, df.Res, n-1),
                        SumsOfSqs = SumsOfSqs,
                        F.Model = c(F.Mod, NA,NA),
                        R2 = SumsOfSqs/SumsOfSqs[length(SumsOfSqs)],
<<<<<<< HEAD
                        se.R2 = c(SD.Mat[c(2:(1+nterms))],NA,NA),
                        P = c(P[c(2:(1+nterms))], NA, NA))
      
    }
    
=======
                        P = c(P[c(2:(1+nterms))], NA, NA))
      
    }
>>>>>>> ef20608ebb44871f24734d4576a97eb6b300bf6d
    rownames(tab) <- c(attr(attr(rhs.frame, "terms"), "term.labels")[u.grps],
                       "Residuals", "Total")
    rownames(tab) <- c(attr(attr(rhs.frame, "terms"), "term.labels")[u.grps],
                       "Residuals", "Total")
    colnames(tab)[ncol(tab)] <- "Pr(>F)"
<<<<<<< HEAD
    if(by=="terms" & all(weights==1)){
      attr(tab, "heading") <- c(vegan:::howHead(attr(perm_mat, "control")),
                                "Terms added sequentially (first to last)\n")
    }
    if(by=="margin" & all(weights==1)){
      attr(tab, "heading") <- c(vegan:::howHead(attr(perm_mat, "control")),
                                "Terms added marginaly (one term for each analysis)\n")
    }
    if(any(weights!=1)){
      attr(tab, "heading") <- c("Weights included (no permutation)")
    }
    
    
=======
    attr(tab, "heading") <- c(howHead(attr(p, "control")),
                              "Terms added sequentially (first to last)\n")
>>>>>>> ef20608ebb44871f24734d4576a97eb6b300bf6d
    class(tab) <- c("anova", class(tab))
    out <- list(aov.tab = tab, call = match.call(),
                model.matrix = rhs, terms = Terms)
    class(out) <- "fast.adonis"
    out
  }


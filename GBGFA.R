############################################################################
#                                                                          #
# Implementation of Gene-wise Prior Bayesian Group Factor Analysis (GBGFA) #
# This implementation is based on the CCAGFA R package:                    #
# https://cran.r-project.org/web/packages/CCAGFA/index.html                #
#                                                                          #
############################################################################

#############################################################################
#
# A wrapper for running the GFA model Nrep times
# and choosing the final model based on the best
# lower bound. This is the recommended way of applying
# the algorithm.
#
#############################################################################
GBGFAexperiment <- function(X,K,opts,Nrep=50) {
  
  if(Nrep==1) {
    return(GBGFA(X,K,opts))
  }
  
  model <- GBGFA(X,K,opts)  
  lb <- tail(model$cost,1)
  
  for(rep in 2:Nrep){
    model_new <- GBGFA(X,K,opts)
    lb_new <- tail(model_new$cost,1)
    
    if(lb_new > lb){
      lb <- lb_new
      model <- model_new
    }
    
    if(opts$verbose==1) {
      print(paste("Run ",rep,"/",Nrep,": ",length(model$cost)," iterations with final cost ",lb[rep],sep=""))
    }
  }
  
  return(model)
}

#############################################################################
#
# The main function for GBGFA
#
# Inputs:
#   X    : List of M data matrices. X[[m]] is a matrix with
#          N rows (samples) and D_m columns (features). The
#          samples need to be co-occurring.
#          NOTE: All of these should be centered, so that the mean
#                of each feature is zero
#          NOTE: The algorithm is roughly invariant to the scale
#                of the data, but extreme values should be avoided.
#                Data with roughly unit variance or similar scale
#                is recommended.
#   K    : The number of components
#   opts : List of options (see function getDefaultOpts())
#
# Output:
# The trained model, which is a list that contains the following elements:
#   Y    : The mean of the latent variables; N times K matrix (in paper transposed)
#   covY : The covariance of the latent variables; K times K matrix
#   YY   : The second moments YY^T; K times K matrix
#
#   W    : List of the mean projections; list of length M of D_i times K matrices
#   covW : List of the covariances of the projections; list of length M, of lists of lengths D[1], D[2],...,D[M] of K times K matrices
#   WW   : List of the second moments WW^T; list of length M of K times K matrices
#
#   tau  : The mean precisions (inverse variance, so 1/tau gives the
#          variances denoted by sigma in the paper); M-element vector
#
#   gp_gamma: The mean precisions of the projection weights, the
#          variances of the ARD prior; vector of length Dsu (Dsu time 1) 
#
#   cost : Vector collecting the variational lower bounds for each
#          iteration
#
#   D    : Data dimensionalities; M-element vector
#   Dsu  : Number of unique features (genes) across all views
#
#   indecies  : Unique identifiers of features within each view; list of M elements, each a D[m] x 2 matrix
#   datavar   : The total variance in the data sets, needed for
#               GFAtrim()
#   addednoise: The level of extra noise as in opts$addednoise
#
#############################################################################
GBGFA <- function(X,K,opts) {
  # Check that data is centered
  if(!all(unlist(lapply(X,colMeans)) < 1e-7) & opts$verbose == 2){
    print("Warning: The data does not have zero mean.")
  }
  
  #
  # Store dimensionalities of data sets 
  #
  M <- length(X)                      # The number of views
  D <- sapply(X,ncol)                 # Collect the number of features in vector D
  Ds <- sum(D)                        # Total number of features
  N <- nrow(X[[1]])                   # The number of samples
  datavar <- vector()                 # The total variance of the data, needed
  for(m in 1:M) {                     #     for scaling in the initialization
    datavar[m]=sum(apply(X[[m]],2,var)) #     and for GFAtrim()
  }
  
  #
  # Get view-indecies (m and d) of each feature
  #
  indecies <- lapply(1:M, function(m){
    ids <- colnames(X[[m]])
    
    # Check if any features are missing identifiers
    if(is.null(ids) | sum(is.na(ids)) > 0){
      stop("Error: Feature identifiers ( colnames(X[[m]]) ) are required...\n")
    }
    
    # Check if any feature identifiers are duplicated
    if(length(ids[duplicated(ids)]) != 0){
      stop("Error: Feature identifiers are not unique...\n")
    }
    
    res <- cbind(m=rep(m, length(ids)), d=seq(1:length(ids)))
    rownames(res) <- ids
    return(res)             
    
  })
  
  ms <- split(do.call('rbind', indecies)[,1],        # unique feature names with their corresponding view i
              rownames(do.call('rbind', indecies)))  # ndecies (m)
  ds <- split(do.call('rbind', indecies)[,2],        # unique feature names with their corresponding feature 
              rownames(do.call('rbind', indecies)))  # indecies in each view (d)
  ufts <- unique(names(ms))                          # unique feature names
  Dsu <- length(ufts)                                # total number of unique features
  freq <- sapply(ms, length)                         # count of views for each feature (gene)
  
  
  # Some constants for speeding up the computation
  const <- - N*Ds/2*log(2*pi) # Constant factors for the lower bound
  Yconst <- sapply(X,function(x){sum(x^2)})
  id <- rep(1,K)                       # Vector of ones for fast matrix calculations
  alpha_0 <- opts$prior.alpha_0        # Easier access for hyperprior values
  beta_0 <- opts$prior.beta_0
  alpha_0t <- opts$prior.alpha_0t
  beta_0t <- opts$prior.beta_0t
  
  #
  # Initialize the model randomly; other initializations could
  # be done, but overdispersed random initialization is quite good.http://belltown.fhcrc.org/rstudio/clear.cache.gif
  # 
  
  # Latent variables
  Y <- matrix(rnorm(N*K,0,1),N,K) # The mean 
  #cat(Y, "\n")
  covY <- diag(1,K)               # The covariance
  YY <- covY + covY*N             # The second moments
  
  # ARD and noise parameters
  gp_gamma <- rep(1,Dsu)             # The mean of the ARD precisions
  loggp_gamma <- rep(1, Dsu)         # The mean of <\log gp_gamma >
  b_gamma <- rep(1, Dsu)            # The parameters of the Gamma distribution
  a_gamma <- alpha_0 + K*freq/2
  
  tau <- rep(opts$init.tau,M)     # The mean noise precisions
  a_tau <- alpha_0t + N*D/2       # The parameters of the Gamma distribution  
  b_tau <- rep(0,M)               #     for the noise precisions
  
  # Alpha needs to be initialized to match the data scale
  gp_gamma <- sapply(ufts, function(ft){
    # Sum over the view-specific components of the prior
    tmp <- sapply(ms[[ft]], function(m){
      K*D[m]/(datavar[m]-1/tau[m])
    })
    return(sum(tmp))
  })
  
  # The projections
  # No need to initialize the projections randomly, since their updating
  # step is the first one; just define the variables here
  W <- vector("list",length=M)    # The means
  covW <- vector("list",length=M) # The covariances
  WW <- vector("list",length=M)   # The second moments
  for(m in 1:M) {
    W[[m]] <- matrix(0,D[m],K)
    
    rownames(W[[m]]) <- colnames(X[[m]])
  
    # now list of Dm matrices of size KxK
    covW[[m]] <- lapply(1:D[m], function(i){
      return(diag(1,K))
    })
    
    # or at once:
    WW[[m]] <- crossprod(W[[m]]) + Reduce('+', covW[[m]])
  }
  
  # Rotation parameters
  R <- diag(K)      # The rotation matrix
  Rinv <- diag(K)   # Its inverse
  r <- as.vector(R) # Vectorized version of R
  
  # parameter list for the optimization function (see ?optim)
  par <- list(K=K,D=D,Ds=Ds,N=N,WW=WW,YY=YY,M=M)
  
  cost <- vector()  # For storing the lower bounds
  
  #
  # The main loop
  #
  
  for(iter in 1:opts$iter.max){ 
    print(paste("Iteration: ", iter, sep=""))
    #
    # Update the projections
    #
    
    for(m in 1:M){
      WW[[m]] <- diag(0,K)
      # Efficient and robust way of computing
      # solve(diag(gp_gamma) + tau * YY^T)
      for(d in 1:D[m]){
        # Get the feature name to access its prior
        i <- which(colnames(X[[m]])[d] == ufts) # index of this feature in gp_gamma (same as index in ufts)  
        
        covW[[m]][[d]] <- chol2inv(chol(YY*tau[m]+diag(gp_gamma[i],K)))
        
        W[[m]][d,] <- crossprod(X[[m]][ , d],Y)%*%covW[[m]][[d]]*tau[m]    
        
      }
      
      WW[[m]] <- crossprod(W[[m]]) + Reduce('+', covW[[m]])
      
    }
    
    # 
    # Update the latent variables
    #
    
    # Efficient and robust way of computing
    # solve(diag(1,K) + tau * WW^t)
    covY <- diag(1,K)
    for(m in 1:M) {
      covY <- covY + tau[m]*WW[[m]]
    }
    
    covY <- chol2inv(chol(covY))
    
    Y <- X[[1]]%*%W[[1]]*tau[1]
    for(m in 2:M){
      Y <- Y + X[[m]]%*%W[[m]]*tau[m]
    }
    Y <- Y%*%covY
    YY <- crossprod(Y) + N*covY
    
    # Update gp_gamma, the ARD parameters
    
    b_gamma <- sapply(1:Dsu, function(i){
      res <- beta_0
      
      for(j in 1:length(ms[[i]])){
        m <- ms[[i]][j]
        d <- ds[[i]][j]
        
        res <- res + (crossprod(W[[m]][d,]) + sum(diag(covW[[m]][[d]])))/2       
      }
      
      return(res)
    })
    
    gp_gamma <- a_gamma/b_gamma
    
    #
    # Update tau, the noise precisions
    #
    for(m in 1:M) {
      b_tau[m] <- beta_0t + ( Yconst[m] + sum(WW[[m]]*YY) - 2*sum(Y*(X[[m]]%*%W[[m]])))/2      
    }
    tau <- a_tau/b_tau
    
    #
    # Calculate the lower bound.
    # Consists of calculating the likelihood term and 
    # KL-divergences between the factorization and the priors
    #
    
    # The likelihood
    logtau <- digamma(a_tau) - log(b_tau) 
    loggp_gamma <- digamma(a_gamma) - log(b_gamma) # updated
    lb.p <- const + sum(N*D/2*logtau) - sum( (b_tau - beta_0t)*tau ) # unchanged
    lb <- lb.p
    
    # E[ ln p(Y)] - E[ ln q(Y) ] , in paper X, unchanged
    lb.px <- - sum(diag(YY))/2   
    lb.qx <- - N*sum(log(svd(covY,nu=0,nv=0)$d))/2 - N*K/2 #N*K*(1+log(2*pi))/2 # - N*K/2
    lb <- lb + lb.px - lb.qx
    
    # E[ ln p(W)] - E[ ln q(W)]
    lb.pw <- K*sum(freq*loggp_gamma)/2  
    for(m in 1:M){              # this sums the second over m the second element and multiplies it by its respective W[[m]]
      mask <- match(colnames(X[[m]]), ufts) # subset gp_gamma/loggp_gamma by the genes in this view
      
      tmp <- sapply(1:nrow(W[[m]]), function(ridx){
        sum(W[[m]][ridx,]^2 + diag(covW[[m]][[ridx]]))
      })
      lb.pw <- lb.pw - sum(tmp*gp_gamma[mask])/2  # gp_gamma[mask] is of length Dm; WW[[m]] is of dim KxK
    }
    
    #lb.qw <- -Ds*K*(1+log(2*pi))/2
    lb.qw <- 0 
    for(m in 1:M){
      for(d in 1:D[m]){
        lb.qw <- lb.qw - sum(log(svd(covW[[m]][[d]],nu=0,nv=0)$d))/2 
      }
    }
    
    lb <- lb + lb.pw - lb.qw
    
    # E[ln p(gp_gamma)] - E[ln q(gp_gamma)]
    lb.pa <- Dsu*( -lgamma(alpha_0) + alpha_0*log(beta_0) ) + (alpha_0-1)*sum(loggp_gamma) - beta_0*sum(gp_gamma)
    lb.qa <- -sum(lgamma(a_gamma)) + sum(a_gamma* log(b_gamma) ) + sum((a_gamma-1)*loggp_gamma) - sum(a_gamma) #sum(b_gamma*gp_gamma)
    lb <- lb + lb.pa - lb.qa
    
    # E[ln p(tau)] - E[ln q(tau)]
    lb.pt <- M*( -lgamma(alpha_0t) + alpha_0t*log(beta_0t) ) + (alpha_0t-1)*sum(logtau) - beta_0t*sum(tau)
    lb.qt <- -sum(lgamma(a_tau)) + sum(a_tau*log(b_tau)) + sum((a_tau-1)*logtau) - sum(a_tau)
    lb <- lb + lb.pt - lb.qt
    
    # Store the cost function 
    cost[iter] <- lb
    
    if(opts$verbose==2) {
      print(paste("Iteration:",iter,"/ cost:",lb))
    }
    
    # Convergence if the relative change in cost is small enough
    if(iter>1){
      diff <- cost[iter]-cost[iter-1]
      if(diff < 0){
        stop("ERRORRR.. lb not monotonically increasing ...\n")
      }
      if( abs(diff)/abs(cost[iter]) < opts$iter.crit | iter == opts$iter.max ){ 
        break
      }
    }
    
  } # the main loop of the algorithm ends
  
  # Add a tiny amount of noise on top of the latent variables,
  # to supress possible artificial structure in components that
  # have effectively been turned off
  Y <- Y + opts$addednoise*matrix(rnorm(N*K,0,1),N,K) %*% chol(covY)
  
  # return the output of the model as a list
  list(W=W,covW=covW,YY=YY,WW=WW,Y=Y,covY=covY,tau=tau,gp_gamma=gp_gamma,cost=cost,D=D,K=K,addednoise=opts$addednoise,datavar=datavar)
}

E <- function(r,par) {
  #
  # Evaluates the (negative) cost function value wrt the transformation
  # matrix R used in the generic optimization routine
  #
  R <- matrix(r,par$K)
  eS <- svd(R)
  
  val <- -sum(par$YY*(tcrossprod(eS$u*outer(rep(1,par$K),1/eS$d^2))),eS$u)/2
  val <- val + (par$Ds-par$N)*sum(log(eS$d))
  for(m in 1:par$M) {
    val <- val - par$D[m]*sum( log( colSums(R*(par$WW[[m]]%*%R)) ))/2
  }
  return(-val)
}

gradE <- function(r,par) {
  #
  # Evaluates the (negative) gradient of the cost function E()
  #
  R <- matrix(r,par$K)
  eS <- svd(R)
  Rinv <- tcrossprod(eS$v*outer(rep(1,par$K),1/eS$d), eS$u)
  gr <- as.vector( tcrossprod( tcrossprod( eS$u*outer(rep(1,par$K),1/eS$d^2) , eS$u )%*%par$YY + diag(par$Ds - par$N,par$K), Rinv) )
  
  tmp1 <- par$WW[[1]]%*%R
  tmp2 <- 1/colSums(R*tmp1)
  tmp1 <- par$D[1]*as.vector( tmp1*outer(rep(1,par$K),tmp2) )
  gr <- gr - tmp1
  for(m in 2:par$M){
    tmp1 <- par$WW[[m]]%*%R
    tmp2 <- 1/colSums(R*tmp1)
    tmp1 <- par$D[m]*as.vector( tmp1*outer(rep(1,par$K),tmp2) )
    gr <- gr - tmp1
  }
  return(-gr)
}

getDefaultOpts <- function(){
  #
  # A function for generating a default set of parameters.
  #
  # To run the algorithm with other values:
  #   opts <- getDefaultOpts()
  #   opts$opt.method <- "BFGS"
  #   model <- GBGFA(X,K,opts)
  
  #
  # Whether to use the rotation explained in the ICML'11 paper.
  # Using the rotation is recommended, but for some data sets
  # equally good solutions can be obtained without rotation and
  # possibly faster.
  #  - TRUE|FALSE
  #
  rotate <- TRUE
  
  #
  # Parameters for controlling how the rotation is solved
  #  - opt.method chooses the optimization method and
  #    takes values "BFGS" or "L-BFGS". The former
  #    is typically faster but takes more memory, so the latter
  #    is the default choice. For small K may use BFGS instead.
  #  - opt.iter is the maximum number of iterations
  #  - lbfgs.factr is convergence criterion for L-BFGS; smaller
  #    values increase the accuracy (10^7 or 10^10 could be tried
  #    to speed things up)
  #  - bfgs.crit is convergence criterion for BFGS; smaller
  #    values increase the accuracy (10^-7 or 10^-3 could also be used)
  #
  opt.method <- "L-BFGS"
  opt.iter <- 10^5
  lbfgs.factr <- 10^3
  bfgs.crit <- 10^-5
  
  #
  # Initial value for the noise precisions. Should be large enough
  # so that the real structure is modeled with components
  # instead of the noise parameters (see Luttinen&Ilin, 2010)
  #  Values: Positive numbers, but generally should use values well
  #          above 1
  #
  init.tau <- 10^3
  
  #
  # Parameters for controlling when the algorithm stops.
  # It stops when the relative difference in the lower bound
  # falls below iter.crit or iter.max iterations have been performed.
  #
  iter.crit <- 10^-6
  iter.max <- 10^5
  
  #
  # Additive noise level for latent variables. The latent variables
  # of inactive components (those with very large gp_gamma) occasionally
  # show some structure in the mean values, even though the distribution
  # matches very accurately the prior N(0,I). This structure disappears
  # is a tiny amount of random noise is added on top of the
  # mean estimates. Setting the value to 0 will make the predictions
  # deterministic
  #
  addednoise <- 1e-5
  
  #
  # Hyperparameters
  # - alpha_0, beta_0 for the ARD precisions
  # - alpha_0t, beta_0t for the residual noise predicions
  #
  prior.alpha_0 <- prior.beta_0 <- 1e-14
  prior.alpha_0t <- prior.beta_0t <- 1e-14
  
  #
  # Verbosity level
  #  0: Nothing
  #  1: Final cost function value for each run of GBGFAexperiment()
  #  2: Cost function values for each iteration
  #
  verbose <- 2
  
  return(list(rotate=rotate, init.tau=init.tau, iter.crit=iter.crit,
              iter.max=iter.max, opt.method=opt.method,
              lbfgs.factr=lbfgs.factr, bfgs.crit=bfgs.crit, opt.iter=opt.iter,
              addednoise=1e-6,
              prior.alpha_0=prior.alpha_0,prior.beta_0=prior.beta_0,
              prior.alpha_0t=prior.alpha_0t,prior.beta_0t=prior.beta_0t,
              verbose=verbose))
}


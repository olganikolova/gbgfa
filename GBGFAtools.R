#############################################################################
#
# This implementation is based on the CCAGFA R package:                    #
# https://cran.r-project.org/web/packages/CCAGFA/index.html                #
#
# Tools for analysing the models learned with the
# functions in GBGFA.R
#
# Contains functions:
#  GFAtrim()   : The same for GFA
#  GFApred()   : Make predictions for new data for GFA
#  GFApred()   : Make predictions for new data for GFA
#
#############################################################################

#############################################################################  
#
# Clean up the model by removing the components that were
# pushed to zero.
#
# Input:
#  model    : A model trained with GFA()
#  threshold: Components that explain less than threshold
#             of the total variation are dropped
#
# Output:
#  model    : The trimmed model (see GFA() for the contents)
#
#############################################################################  
GFAtrim <- function(model,threshold=1e-3) {

  M <- length(model$W)
  N <- dim(model$Y)[1]

  active <- matrix(1,M,model$K) 
  for(m in 1:M) {
    # The relative contribution to describing the total data
    # variation is D[m]/gp_gamma[m,k] / (datavar[m] - D[m]/tau[m])
    residual <- model$datavar[m] - model$D[m]/model$tau[m]
    if(residual < 1e-2*model$datavar[m]) {  # Noise models almost everything
      print(paste("Warning, trimming a model for which the noise explains already",
            format(model$D[m]/model$tau[m]/model$datavar[m],digits=4),"percent of the total variation in data set",m,". Might result for zero components for that data set."))
      residual <- 1e-2*model$datavar[m]
    }

    active[m,which(model$D[m]/model$gp_gamma[m,]<threshold*residual)] <- 0
  }

  keep <- which(colSums(active)>0)
  model$Y <- model$Y[,keep]
  model$covY <- model$covY[keep,keep]

  for(m in 1:M) {
    model$W[[m]] <- model$W[[m]][,keep]
    model$covW[[m]] <- model$covW[[m]][keep,keep]
    model$WW[[m]] <- model$WW[[m]][keep,keep]
  }
  model$gp_gamma <- model$gp_gamma[,keep]

  active <- active[,keep]
  model$K <- length(keep)

  model$YY <- crossprod(model$Y) + N*model$covY
  id <- rep(1,model$K)   # Vector of ones for fast matrix calculations
  for(m in 1:M) {
    for(k in 1:model$K) {
      if(active[m,k]==0) {
        model$W[[m]][,k] <- 0
      }
    }
    tmp <- 1/sqrt(model$gp_gamma[m,])
    model$covW[[m]] <- 1/model$tau[m] * outer(tmp,tmp) *
      chol2inv(chol(outer(tmp,tmp)*model$YY + diag(1/model$tau[m],model$K)))

    # An alternative way; could be tried in case of 
    #   issues with numerical stability
    #eS <- eigen( outer( tmp, id )*model$YY*outer(id,tmp) + diag(1/model$tau[m],model$K) , symmetric=TRUE)
    #model$covW[[m]] <- 1/tau[m] * outer( tmp, id ) * tcrossprod( eS$vectors*outer( id, 1/eS$values), eS$vectors ) * outer(id,tmp)

    model$WW[[m]] <- crossprod(model$W[[m]]) + model$covW[[m]]*model$D[m]
  }

  model$active <- active
  model$trimmed <- TRUE

  return(model)
}
  
#############################################################################  
#
# Function for making predictions with the model. Gives the
# mean prediction and the mean and covariance of the latent
# variables. The predictive distribution itself does not have
# a closed-form expression, so the function also allows drawing
# samples from it.
#
# Inputs:
#   pred:  Binary vector of length 2, indicating which of the
#          two data sets have been observed. (1,0) indicates
#          we observe the first data set and want to predict
#          the values for the latter, and (0,1) does the opposite.
#          Using (1,1) allows computing the latent variables
#          for new test samples where both views are observed.
#   X   :  The test data as a list of length 2, given in the
#          same format as for the function GFA(). The data
#          matrix for the missing views can be anything, e.g.
#          zeros, but it needs to exist
#   model: A model learned from training data using GFA()
#   sample: Should we sample observations from the full predictive
#           distribution?
#   nSample: How many samples to draw if sample==TRUE
#
#
# Outputs:
# A list containing:
#   X    : The mean predictions as list. Observed data sets are retained
#          as they were.
#   Y    : Mean latent variables of the test samples, given the observed
#          data; N times K matrix
#   covY : Covariance of the latent variables; K times K matrix
#   sam  : Samples drawn from the predictive distribution, only
#          returned if sample==TRUE. A list of Y, W and X.
#          Y is nSample times N times K matrix of the samples values.
#          W and X are M-element lists where only the predicted
#          views are included (to avoid storing nSample identical
#          copies of the observed data), each being a multidimensional
#          array of nSample times the size of W and X, respectively.
#
#############################################################################  
GBGFApred <- function(pred,X,model,sample=FALSE,nSample=100) {
  
  tr <- which(pred==1)  # The observed data sets
  pr <- which(pred==0)  # The data sets that need to be predicted

  N <- nrow(X[[tr[1]]])
  M <- length(model$D)
  
  # Estimate the covariance of the latent variables
  covY <- diag(1,model$K)
  for(m in tr) {
    covY <- covY + model$tau[m]*model$WW[[m]]
  }

  # Estimate the latent variables
  eS <- eigen( covY ,symmetric=TRUE)
  covY <- tcrossprod( eS$vectors*outer(rep(1,model$K),1/eS$values), eS$vectors )
  Y <- matrix(0,N,model$K)
  for(m in tr) {
    Y <- Y + X[[m]]%*%model$W[[m]]*model$tau[m]
  }
  Y <- Y%*%covY

  # Add a tiny amount of noise on top of the latent variables,
  # to supress possible artificial structure in components that
  # have effectively been turned off
  Y <- Y + model$addednoise*matrix(rnorm(N*model$K,0,1),N,model$K) %*% chol(covY)

  # The prediction
  # NOTE: The ICML'11 paper has a typo in the prediction formula
  # on page 5. The mean prediction should have W_2^T instead of W_2.
  for(m in pr) {
    X[[m]] <- tcrossprod(Y,model$W[[m]])
  }

  # Sample from the predictive distribution
  # Note that this code is fairly slow fow large nSample
  if(sample) {
    sam <- list()
    sam$Y <- array(0,c(nSample,N,model$K))
    sam$X <- vector("list",length=M)
    sam$W <- vector("list",length=M)
    cholW <- vector("list",length=M)
    for(m in pr) {
      cholW[[m]] <- chol(model$covW[[m]])
      sam$W[[m]] <- array(0,c(nSample,model$D[m],model$K))
      sam$X[[m]] <- array(0,c(nSample,N,model$D[m]))
    }

    cholZ <- chol(covY)
    for(i in 1:nSample) {
      Ztemp <- Y + matrix(rnorm(N*model$K,0,1),N,model$K) %*% cholZ
      sam$Y[i,,] <- Ztemp
      for(m in pr) {
        Wtemp <- model$W[[m]] + matrix(rnorm(model$K*model$D[[m]],0,1),model$D[m],model$K) %*% cholW[[m]]
        sam$W[[m]][i,,] <- Wtemp
        sam$X[[m]][i,,] <- tcrossprod(Ztemp,Wtemp) + matrix(rnorm(N*model$D[m],0,1/sqrt(model$tau[m])),N,model$D[m])
      }
    }
  }
  
  if(sample)
    return(list(X=X,Y=Y,covY=covY,sam=sam))
  else
    return(list(X=X,Y=Y,covY=covY))
}
#############################################################################  
#
# Generate data from a trained GFA model
#
# Inputs:
#  model : a model trained by GFA()
#  N     : the number of samples to be generated
#
# Outputs:
# A list containing the following elements
#  X : A list of M elements, N times D[m] matrices
#  Y : The latent variables used to generate the data, N times K
#
#############################################################################  
GBGFAsample <- function(model,N) {

  # Latent variables are white noise
  Y <- matrix(rnorm(N*model$K,0,1),N,model$K)
  
  X <- vector("list",length=length(model$D))
  for(view in 1:length(model$D)) {
    # The mean is given by W times Y, and the noise is diagonal
    X[[view]] <- Y %*% t(model$W[[view]]) +
      matrix(rnorm(N*model$D[view],0,1/sqrt(model$tau[view])),N,model$D[view])
  }

  return(list(X=X, Y=Y))
}


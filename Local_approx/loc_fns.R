# Additive model (1): y = f(x) + g(Z) + e  
library(MASS)

# Generates samples from f: y = f(x) + g(Z) + e
GP_gen <- function(N, y, X, newX) { 
  X <- as.matrix(X)
  d <- ncol(X)
  newX <- as.matrix(newX)
  pars <- GP.lap(y, X)  # with Lap
  rho_x <- pars["rho[1]"]^2
  Kx <- rbfkernel(as.matrix(newX[ ,1]), rho_x, as.matrix(X[ ,1]))
  Kxx <- rbfkernel(as.matrix(newX[ ,1]), rho_x)
  K <- 0
  
  for(i in 1:d) {  # additive kernel
    Ki <- rbfkernel(as.matrix(X[ ,i]), pars[i+1]^2)
    K <- Ki + K
  }
  diag(K) <- diag(K) + pars["sigma"]^2
  invK <- solve(K)
  
  Mu <- pars["mu"] + Kx %*% invK %*% y
  Sigma <- Kxx - Kx %*% invK %*% t(Kx)
  sampled_fs <- mvrnorm(N, Mu, Sigma)
  
  return(sampled_fs)
}

# Same as before but sampling from hyperparameters' posterior
GP_gen_mcmc <- function(N, y, X, newX) { 
  X <- as.matrix(X)
  d <- ncol(X)
  rhonames <- paste0("rho[", 1:d, "]")
  newX <- as.matrix(newX)
  standata <- list(N_obs = length(y), d = d, X = X, y_obs = y)
  stanfit <- sampling(GP_mod, data = standata,  # fit additive GP
                      iter = 1000+N, warmup = 1000, chains = 1, refresh = 0)
  all.samples <- rstan::extract(stanfit, pars = c("sigma", rhonames), inc_warmup = F)  # "mu",
  sampled_fs <- matrix(NA, nrow = N, ncol = nrow(newX)) 
  
  for(t in 1:N) {
    tsamples <- lapply(all.samples, `[[`, t)
    mu <- 0 # tsamples$mu
    sigma <- tsamples$sigma
    rho_x <- as.numeric(tsamples[rhonames[1]])^2
    Kx <- rbfkernel(as.matrix(newX[ ,1]), rho_x, as.matrix(X[ ,1]))
    Kxx <- rbfkernel(as.matrix(newX[ ,1]), rho_x)
    K <- 0
    
    for(i in 1:d) {
      rho_i <- as.numeric(tsamples[rhonames[i]])
      Ki <- rbfkernel(as.matrix(X[ ,i]), rho_i^2)
      K <- Ki + K
    }
    diag(K) <- diag(K) + sigma^2
    invK <- solve(K)
    
    Mu <- mu + Kx %*% invK %*% y
    Sigma <- Kxx - Kx %*% invK %*% t(Kx)
    sampled_fs[t, ] <- mvrnorm(1, Mu, Sigma)
  }
  
  return(sampled_fs)
}

rbfkernel <- function(X, sigma = 1, Y = NULL) {
  # test if X is a matrix
  if(!is.matrix(X)) {
    print("X must be a matrix containing samples in its rows")
    return()
  }
  # test if sigma is a number and > 0
  if(length(sigma) != 1 || sigma <= 0) {
    print("sigma must be a number > 0 specifying the rbf-kernel width")
    return()
  }
  
  if(!is.null(Y)) {
    # test if Y is a matrix
    if(!is.matrix(Y)) {
      print("Y must be a matrix containing samples in its rows or NULL if it should not be used")
      return()
    }
    # test if vectors in X and Y have same dimension
    if(ncol(X) != ncol(Y)) {
      print("The samples in the rows of X and Y must be of same dimension")
      return()
    }
  }
  n <- nrow(X) # number of samples in X
  
  if(is.null(Y)) {
    # calculate distance matrix
    XtX <- tcrossprod(X)
    XX <- matrix(1, n) %*% diag(XtX)
    D <- XX - 2*XtX + t(XX) # distance matrix
  } else {
    m <- nrow(Y) # number of samples in Y
    # calculate distance matrix (between vectors of X and Y)
    XX <- matrix(apply(X ^ 2, 1, sum), n, m)
    YY <- matrix(apply(Y ^ 2, 1, sum), n, m, byrow = TRUE)
    XY <- tcrossprod(X, Y)
    D <- XX - 2*XY + YY
  }
  # calculate rbf-kernel matrix
  K <- exp(-D/(2*sigma))
  
  return(K)
}

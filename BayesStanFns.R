library(rstan)
library(bridgesampling)
options(mc.cores = parallel::detectCores())

# Compile stan files
Gauss_mod <- stan_model("Gauss.stan") 
GP_mod <- stan_model("Add.stan") 

# Laplace approximation of log marginal likelihood
Laplace <- function(optim) {
  d <- length(optim$par)
  penalty <- determinant(-optim$hessian, logarithm = T)$modulus - d*log(2*pi)
  optim$value - 0.5*penalty
}

# Optimizes posterior of Gaussian and computes Laplace approximation
Gauss.lap <- function(y) {  
  Sdata <- list(N_obs = length(y), y_obs = y)
  init <- list(mu = 0, sigma = 1)
  
  opt <- optimizing(Gauss_mod, data = Sdata, init = init, hessian = T)
  lik <- setNames(Laplace(opt), "newscore")
  c(lik, opt$par)
}

# Optimizes posterior of additive GP
GP.lap <- function(y, X) {  
  X <- as.matrix(X)
  d <- ncol(X)
  Sdata <- list(N_obs = nrow(X), d = d, X = X, y_obs = y)  
  init <- list(rho = as.array(rep(1, d)), mu = 0, sigma = 1)
  
  opt <- optimizing(GP_mod, data = Sdata, init = init, hessian = T)
  lik <- setNames(Laplace(opt), "newscore")
  c(lik, opt$par)
}

# Bridgesampling for marginal likelihoods
Gauss.mcmc <- function(y) {
  N <- length(y)
  stanfit <- sampling(Gauss_mod, data = list(N_obs = N, y_obs = y), 
                      iter = 500, warmup = 200, chains = 1, refresh = 0)
  lik <- setNames(bridge_sampler(stanfit, silent = TRUE)$logml, "newscore")
  c(lik, Gauss.lap(y)[-1])
}

GP.mcmc <- function(y, X) { 
  N <- length(y)
  X <- as.matrix(X)
  d <- ncol(X)
  stanfit <- sampling(GP_mod, data = list(N_obs = N, d = d, X = X, y_obs = y),  
                      iter = 500, warmup = 200, chains = 1, refresh = 0)
  lik <- setNames(bridge_sampler(stanfit, silent = TRUE)$logml, "newscore")
  c(lik, GP.lap(y, X)[-1])
}

# Weighted generic quantile estimator
# https://aakinshin.net/posts/weighted-quantiles/
wquantile.generic <- function(x, probs, cdf.gen, weights = NA) {
  n <- length(x)
  if (any(is.na(weights)))
    weights <- rep(1 / n, n)
  nw <- sum(weights)^2 / sum(weights^2) # Kish's effective sample size
  
  indexes <- order(x)
  x <- x[indexes]
  weights <- weights[indexes]
  
  weights <- weights / sum(weights)
  cdf.probs <- cumsum(c(0, weights))
  
  sapply(probs, function(p) {
    cdf <- cdf.gen(nw, p)
    q <- cdf(cdf.probs)
    w <- tail(q, -1) - head(q, -1)
    sum(w * x)
  })
}

# Weighted Harrell-Davis quantile estimator
whdquantile <- function(x, probs, weights = NA) {
  cdf.gen <- function(n, p) return(function(cdf.probs) {
    pbeta(cdf.probs, (n + 1) * p, (n + 1) * (1 - p))
  })
  wquantile.generic(x, probs, cdf.gen, weights)
}

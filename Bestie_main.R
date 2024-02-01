# This makes Figure 3 but with Bestie
library(tidyverse)
library(BiDAG)
library(matrixStats)
library(MASS)
library(gridExtra)
library(HDInterval)

source("/Users/giudic0000/Downloads/Nonlinear_CausalFx/Fourier_fns.R")
source("/Users/giudic0000/Downloads/Nonlinear_CausalFx/BayesStanFns.R")
source("/Users/giudic0000/Downloads/Nonlinear_CausalFx/fxsampling_fns.R")
insertSource("~/Downloads/Nonlinear_CausalFx/GPscore.R", package = "BiDAG")

set.seed(99)
n <- 5      # number of nodes
N <- 50     # number of samples
lambda <- 1 # non-linearity parameter
grid_size <- 101  # steps in x axis

myDAG <- pcalg::randomDAG(n, prob = 0.4, lB = 1, uB = 2) 
trueDAG <- as(myDAG, "matrix")
truegraph <- 1*(trueDAG != 0)
knowndag <- T  # is graph known?
Fou_result <- Fou_nldata(truegraph, N, lambda = lambda, noise.sd = 0.5, standardize = T) 
data <- Fou_result$data

true.fx <- ComputeAllTrueFx(truegraph, weights = Fou_result$true.weights, 
                            scales = Fou_result$scales, mcsamples = 10000, 
                            noise.sd = 0.5, grid_size = grid_size)
testdata <- apply(data, 2, 
                  FUN = function(x) seq(-2, 2, length.out = grid_size))  # are you sure?? we are getting fx from -2 to 2

library(Bestie)
library(BiDAG)
library(MASS)
library(LaplacesDemon)

source("/Users/giudic0000/Downloads/Nonlinear_CausalFx/fxsampling_fns.R")

# This functions returns the matrix of sampled coefficients
sample_coefs <- function(adj, scoreParam) {
  n <- scoreParam$n
  res <- DAGparameters(adj, scoreParam)
  B <- adj * 0
  
  for (i in 1:n) {
    parents <- which(adj[ ,i] == 1)  # parents of i
    
    if(length(parents) > 0) {
      Mu <- res$mus[[i]]
      Sigma <- res$sigmas[[i]]
      Sigma <- (Sigma + t(Sigma))/2  # force to be symmetric...
      df <- res$dfs[[i]]
      
      B[parents,i] <- rmvt(1, Mu, Sigma, df)
    }
  }
  return(B)
}

# This function computes one sample from E(Y|do(X=0)) for different Y and a given DAG
sample_intercepts <- function(adj, x, scoreParam, beta_mat) {
  nu_prime <- scoreParam$muN
  am_prime <- scoreParam$N + scoreParam$am
  aw_prime <- scoreParam$awpN
  R <- scoreParam$TN
  
  W <- rWishart(1, aw_prime, solve(R))  # sample from parameter posterior (according to section B.1 in Gadget)
  mu <- mvrnorm(1, nu_prime, solve(am_prime*W[,,1]))
  
  top_order <- rev(BiDAG:::DAGtopartition(n, adj)$permy)  # go down order
  x_ind <- which(top_order == x)
  curr_levels <- rep(0, n)  # to store the sampled intercepts
  
  for(y in top_order[-x_ind]) {
    parents <- which(adj[, y] == 1)  # find parents of y
    
    if(length(parents) > 0) {  # one or more parents
      curr_levels[y] <- mu[y] + beta_mat[parents,y] %*% curr_levels[parents]  
    }
  }
  return(curr_levels)
}

set.seed(101)
scoreParam <- scoreparameters("bge", data)
bge.searchspace <- set.searchspace(data, dual = F, "bge")
fit <- bge.partition.mcmc(bge.searchspace, order = F, iterations = 2500)  
ndags <- length(fit$traceadd$incidence)
x_vals <- seq(-2, 2, length.out = grid_size)
testdata <- apply(data, 2, 
                  FUN = function(x) seq(-2, 2, length.out = grid_size))

lay <- matrix(NA, nrow = n, ncol = n)  # plot layout
plots <- list()  # to store all the plots 
plot.counter <- 1

# Get all linear causal fx for do(x)
for(x in 1:n) {
  all_fx <- vector(mode = "list", length = ndags)
  
  for(k in 1:ndags) {
    dag_fx <- matrix(0, n, grid_size)  # to store the samples
    dag <- as.matrix(fit$traceadd$incidence[[k]])
    
    B <- sample_coefs(dag, scoreParam)  # matrix of sampled coefficients
    A <- sample_intercepts(dag, x, scoreParam, B)  # vector of sampled intercepts
    C <- solve(diag(n) - B)  # matrix of causal fx
    
    for(i in 1:n) {
      dag_fx[i, ] <- A[i] + C[x,i] * x_vals  # E(i|do(x = xvals))
    }
    all_fx[[k]] <- dag_fx
  }
  
  # Plot bge results for different y (given x)
  for(i in 1:(n-1)) {
    y <- (1:n)[-x][i]
    xy_effect <- matrix(NA, nrow = grid_size, ncol = ndags)
    
    for(dag in 1:ndags) {
      arr <- all_fx[[dag]]
      xy_effect[ ,dag] <- arr[y, ]
    }
    E_val <- apply(xy_effect, 1, mean)  # expected value of effect
    quants <- apply(xy_effect, 1, quantile, probs = c(0.1, 0.9))
    
    fitted.vals <- data.frame(x = testdata[ ,x], E_val = E_val, 
                              lo = quants[1, ], hi = quants[2, ], xy_effect) %>%
      pivot_longer(cols = !c(x, E_val, lo, hi))
    
    true.vals <- data.frame(x = testdata[ ,x], value = true.fx[x,y, ], name = rep(NA, grid_size))
    plots[[plot.counter]] <- ggplot(data = fitted.vals, aes(x, value, group = name)) +
      geom_line(alpha = 0.2, col = '#c2c2c2') +
      geom_line(aes(x, E_val), color = "#be4173", linetype = "dashed") +
      geom_ribbon(aes(ymin = lo, ymax = hi, group = NULL, color = NULL), 
                  alpha = 0.5, fill = "#f34444") + 
      geom_line(data = true.vals, aes(color = NULL), color = '#21ba35') +
      ggtitle(as.expression(bquote("E(X"[.(y)]*"|"*"do(X"[.(x)]*"))"))) +
      theme_minimal() +
      ylim(-2.42, 2.3) + 
      theme(legend.position = "none",
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            plot.title = element_text(hjust = 0.5, size = 11))
    
    lay[y,x] <- plot.counter  # position on grid
    plot.counter <- plot.counter + 1
  }
}

grid.arrange(grobs = plots, layout_matrix = lay)


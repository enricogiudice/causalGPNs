library(Bestie)
library(BiDAG)
library(MASS)
library(tidyverse)
library(janitor)
library(gridExtra)
library(LaplacesDemon)

source("fxsampling_fns.R")

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
      Sigma <- (Sigma + t(Sigma))/2 
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
  
  W <- rWishart(1, aw_prime, solve(R))  # sample from parameter posterior 
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


# Apply to arabidopsis data
arabidopsis <- read.csv("A_thaliana_figures/arabidopsis.txt", sep="")
some_genes <- c("PPDS2","PPDS1","GPPS","IPPI1","HDR","HDS",
                "MECPS","CMK","MCT","DXR","DXPS3","DXPS2","DXPS1")

arabidopsis %>%
  filter(Genename %in% some_genes) %>%
  dplyr::select(-(1:6)[-4]) %>% 
  t() %>%
  row_to_names(1) %>%
  `class<-` ("numeric") %>%
  apply(2, log2) %>%
  scale() -> data
n <- ncol(data)

# Sample DAGs
set.seed(101)
scoreParam <- scoreparameters("bge", data)
bge.searchspace <- set.searchspace(data, dual = F, "bge")
bge.fit <- bge.partition.mcmc(bge.searchspace, order = F, iterations = 5000)  
ndags <- length(bge.fit$traceadd$incidence)

# Sample E(Y|do(X=x))
grid_size <- 101
x_vals <- seq(-2, 2, length.out = grid_size)

do <- "PPDS1"  # intervention vab
x <- which(colnames(data) == do)
all_fx <- vector(mode = "list", length = ndags)

for(k in 1:ndags) {
  dag_fx <- matrix(0, n, grid_size)  # to store the samples
  dag <- as.matrix(bge.fit$traceadd$incidence[[k]])
  
  B <- sample_coefs(dag, scoreParam)  # matrix of sampled coefficients
  A <- sample_intercepts(dag, x, scoreParam, B)  # vector of sampled intercepts
  C <- solve(diag(n) - B)  # matrix of causal fx
  
  for(i in 1:n) {
    dag_fx[i, ] <- A[i] + C[x,i] * x_vals  # E(i|do(x = xvals))
  }
  all_fx[[k]] <- dag_fx
}

# Plot bge results
testdata <- apply(data, 2, 
                  FUN = function(x) seq(-2, 2, length.out = grid_size))
plots <- list()  # to store the plots

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

  plots[[i]] <- ggplot(data = fitted.vals, aes(x, value, group = name)) +
    geom_line(alpha = 0.2, color = '#d4d4d4') +
    geom_line(aes(x, E_val), color = "#be4173", linetype = "dashed") +
    geom_ribbon(aes(ymin = lo, ymax = hi, group = NULL, color = NULL), 
                alpha = 0.5, fill = "#f34444") + 
    ggtitle(paste0("E(", colnames(data)[y], "|do(", do, "))")) +
    theme_minimal() +
    coord_cartesian(ylim = c(-4.5, 4), xlim = c(-1.8, 1.8)) +
    theme(axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          plot.title = element_text(hjust = 0.5, size = 8.5),
          plot.margin = margin(3,6,3,3)) 
}

# Visualize all plots
lay <- rbind(matrix(1:(n-1), ncol = 6))
grid.arrange(grobs = plots, layout_matrix = lay)

# size = 3.33 x 9

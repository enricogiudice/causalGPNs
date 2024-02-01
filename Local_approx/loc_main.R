library(rstan)
library(tidyverse)
library(BiDAG)
library(matrixStats)
library(gridExtra)
library(questionr)

source("/Users/giudic0000/Downloads/Nonlinear_CausalFx/Local approach/loc_fns.R")
source("/Users/giudic0000/Downloads/Nonlinear_CausalFx/Fourier_fns.R")
source("/Users/giudic0000/Downloads/Nonlinear_CausalFx/BayesStanFns.R")
source("/Users/giudic0000/Downloads/Nonlinear_CausalFx/Local approach/sampling_fns.R")
insertSource("~/Downloads/Nonlinear_CausalFx/GPscore.R", package = "BiDAG")

set.seed(99)  # was 101
n <- 5      # number of nodes
N <- 50     # number of samples
lambda <- 1 # non-linearity parameter
grid_size <- 101  # steps in x axis

myDAG <- pcalg::randomDAG(n, prob = 0.4, lB = 1, uB = 2) 
trueDAG <- as(myDAG, "matrix")
truegraph <- 1*(trueDAG != 0)
knowndag <- F  # is graph known?
Fou_result <- Fou_nldata(truegraph, N, lambda = lambda, noise.sd = 0.5, standardize = T) 
data <- Fou_result$data

true.fx <- ComputeAllTrueFx(truegraph, weights = Fou_result$true.weights, 
                            scales = Fou_result$scales, mcsamples = 10000, 
                            noise.sd = 0.5, grid_size = grid_size)

post.samples <- 500  # number of total posterior samples to generate
testdata <- apply(data, 2, 
                  FUN = function(x) seq(-2, 2, length.out = grid_size))

# Get mcmc samples
GP.searchspace <- set.searchspace(data, dual = F, method = "GP")
fit <- GP.partition.mcmc(data, GP.searchspace, order = F, burnin = 0.2, iterations = 2500)

# Make list of sampled parents for each variable
ndags <- length(fit$trace)
mcmc_mats <- fit$traceadd$incidence  # all sampled adjacency matrices 

if(knowndag) {
  mcmc_mats <- lapply(mcmc_mats, function(x) truegraph)  # sample only correct DAG to check 
}
fx_mats <- lapply(mcmc_mats,  # save non-zero causal effects for later
                  function(x) as(solve(diag(n) - x), "dgCMatrix"))

# Loop over all ordered pairs
all.pairs <- subset(expand.grid(rep(list(1:n), 2)), Var1 != Var2)
par_lengths <- rep(NA, nrow(all.pairs))  # save number of smapled parents
plots <- list()  # to store the plots
lay <- matrix(NA, nrow = n, ncol = n)  # plot layout

for(i in 1:nrow(all.pairs)) {
  x <- all.pairs[i,1]  # intervention vab
  y <- all.pairs[i,2]  # outcome vab
  
  # Estimate the empirical distribution of parent sets
  sets <- rep(NA, ndags)  # here we will store the parent sets as strings
  
  for(k in 1:ndags) {  # loop over all sampled dags

    if(fx_mats[[k]][x,y] == 0) {  # no causal effect
           pax <- "null"
    }
    
    else {
      pax <- mcmc_mats[[k]][,x]  # parents of x indicator
      pax <- paste(pax, collapse = "")  # as string
    }
    sets[k] <- pax
  }
  p.sets <- wtd.table(sets, weights = exp(fit$weights))
  n.samples <- table(sample(names(p.sets), post.samples, replace = T, prob = p.sets))  # number of samples for every parent set
  par_lengths[i] <- length(n.samples)
  
  # Get posterior for every parent set
  yD <- data[ ,y]  # data
  xy_effect <- NULL  # to store the samples from the posterior
  
  # Compute causal effects between all pairs of variables
  for(pa in names(n.samples)) {
    
    if(pa == "null") {  # no causal effect
      smus <- rt(n.samples[pa], N-1) * sd(yD) / N + mean(yD)  # Estimate for E(y)
      xy_effect <- rbind(xy_effect, replicate(grid_size, smus))
    }
    
    else {
      split.pa <- strsplit(pa, "")[[1]]  # back to numeric
      PaX <- as.numeric(split.pa)
      i.cX <- c(x, which(PaX == 1))  # x and its parents 
      cX <- data[ ,i.cX]
      testX <- testdata[ ,i.cX]
      xy_effect <- rbind(xy_effect, GP_gen_mcmc(n.samples[pa], yD, cX, testX))  # with MCMC
    }
  }
  
  # Plot the result 
  E_val <- apply(xy_effect, 2, mean)  # expected value of effect
  #SD_val <- apply(xy_effect, 2, sd)  # sd of effect
  quants <- apply(xy_effect, 2, quantile, probs = c(0.1, 0.9))
  
  fitted.vals <- data.frame(x = testdata[ ,x], t(xy_effect), E_val = E_val, 
                            lo = quants[1, ], hi = quants[2, ]) %>%
    pivot_longer(cols = !c(x, E_val, lo, hi))
  true.vals <- data.frame(x = testdata[ ,x], value = true.fx[x,y, ], name = rep(NA, grid_size))
  
  plots[[i]] <- ggplot(data = fitted.vals, aes(x, value, group = name)) +
    geom_line(alpha = 0.2, col = '#d6d6d6') +
    geom_ribbon(aes(ymin = lo, ymax = hi, group = NULL), 
                alpha = 0.5, fill = "#f34444") + 
    geom_line(data = true.vals, color = '#21ba35') +
    geom_line(aes(x, E_val), color = "#be4173", linetype = "dashed") +
    ggtitle(as.expression(bquote("E(X"[.(y)]*"|"*"do(X"[.(x)]*"))"))) +
    ylim(-2.4, 2.3) +
    theme_minimal() +
    theme(axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          plot.title = element_text(hjust = 0.5, size = 11)) 
  
  lay[y,x] <- i  # position on grid
}

# Visualize all plots
grid.arrange(grobs = plots, layout_matrix = lay)
plot_grid <- arrangeGrob(grobs = plots, layout_matrix = lay)

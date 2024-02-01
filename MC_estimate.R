library(tidyverse)
library(BiDAG)
library(matrixStats)
library(MASS)
library(gridExtra)
library(HDInterval)

source("Fourier_fns.R")
source("BayesStanFns.R")
source("fxsampling_fns.R")
insertSource("GPscore.R", package = "BiDAG")

set.seed(99)
n <- 5      # number of nodes
N <- 50     # number of samples
lambda <- 1 # non-linearity parameter
grid_size <- 101  # steps in x axis

myDAG <- pcalg::randomDAG(n, prob = 0.4, lB = 1, uB = 2) 
trueDAG <- as(myDAG, "matrix")
truegraph <- 1*(trueDAG != 0)
knowndag <- T  # known/unknwon graph
Fou_result <- Fou_nldata(truegraph, N, lambda = lambda, noise.sd = 0.5, standardize = T) 
data <- Fou_result$data

true.fx <- ComputeAllTrueFx(truegraph, weights = Fou_result$true.weights, 
                            scales = Fou_result$scales, mcsamples = 10000, 
                            noise.sd = 0.5, grid_size = grid_size)
testdata <- apply(data, 2, 
                  FUN = function(x) seq(-2, 2, length.out = grid_size))

GP.searchspace = set.searchspace(data, dual = F, "GP")
fit <- GP.mcmc.fx(data, GP.searchspace, order = F, iterations = 2500, mcsamples = 1, 
                  truedag = if(knowndag) truegraph else NULL)
weights <- exp(fit$weights)

# Plot results
all.pairs <- subset(expand.grid(rep(list(1:n), 2)), Var1 != Var2)
plots <- list()  # to store the plots
lay <- matrix(NA, nrow = n, ncol = n)  # plot layout

for(i in 1:nrow(all.pairs)) {
  x <- all.pairs[i,1]  # intervention vab
  y <- all.pairs[i,2]  # outcome vab
  
  xy_effect <- table_effect(fit$fx, x, y)  # get matrix of effects
  E_val <- apply(xy_effect, 1, weighted.mean, w = weights)  # expected value of effect
  quants <- apply(xy_effect, 1, whdquantile, probs = c(0.1, 0.9), weights = weights)
  fitted.vals <- data.frame(x = testdata[ ,x], E_val = E_val, 
                            lo = quants[1, ], hi = quants[2, ], xy_effect) %>%
    pivot_longer(cols = !c(x, E_val, lo, hi)) %>%
    mutate(weights = rep(weights, grid_size))
  true.vals <- data.frame(x = testdata[ ,x], value = true.fx[x,y, ], name = rep(NA, grid_size))
  
  plots[[i]] <- ggplot(data = fitted.vals, aes(x, value, group = name, color = weights)) +
    geom_line(alpha = 0.2) +
    scale_color_gradient(low='#d6d6d6', high='#d6d6d6') +
    {if(!knowndag) 
      scale_color_gradient2(low='#a6a6a6', high='#545454', breaks = c(0.002, 0.021))} +
    geom_ribbon(aes(ymin = lo, ymax = hi, group = NULL, color = NULL), 
                alpha = 0.5, fill = "#f34444") + 
    geom_line(data = true.vals, aes(color = NULL), color = '#21ba35') +
    geom_line(aes(x, E_val), color = "#be4173", linetype = "dashed") +
    ggtitle(as.expression(bquote("E(X"[.(y)]*"|"*"do(X"[.(x)]*"))"))) +
    theme_minimal() +
    ylim(-2.508, 2.3) + 
    theme(legend.position = "bottom",
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          plot.title = element_text(hjust = 0.5, size = 11))
  mylegend <- g_legend(plots[[i]])
  plots[[i]] <- plots[[i]] + theme(legend.position="none")
  
  lay[y,x] <- i  # position on grid
}

if(!knowndag) {
  lay <- rbind(lay, c(NA, NA, 21, NA, NA))
  plots[[21]] <- mylegend
}

# Visualize all plots
grid.arrange(grobs = plots, layout_matrix = lay)
plot_grid <- arrangeGrob(grobs = plots, layout_matrix = lay)

# size = 9 x 7 | 8.4

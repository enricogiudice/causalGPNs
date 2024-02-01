library(tidyverse)
library(janitor)
library(BiDAG)
library(MASS)
library(matrixStats)
library(questionr)
library(gridExtra)

source("BayesStanFns.R")
source("fxsampling_fns.R")
insertSource("GPscore.R", package = "BiDAG")

arabidopsis <- read.csv("A_thaliana/arabidopsis.txt", sep="")
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


# MC approach 
set.seed(101)
GP.searchspace = set.searchspace(data, dual = F, "GP")
fit <- GP.mcmc.fx(data, GP.searchspace, order = F, iterations = 5000, 
                  mcsamples = 1, truedag = NULL)
weights <- exp(fit$weights)

grid_size <- 101
testdata <- apply(data, 2, 
                  FUN = function(x) seq(-2, 2, length.out = grid_size))

do <- "PPDS1"  # intervention vab
plots <- list()  # to store the plots

for(i in 1:(n-1)) {
  x <- which(colnames(data) == do)
  y <- (1:n)[-x][i]
  xy_effect <- table_effect(fit$fx, x, y)  # get matrix of effects
  E_val <- apply(xy_effect, 1, weighted.mean, w = weights)  # expected value of effect
  quants <- apply(xy_effect, 1, whdquantile, probs = c(0.1, 0.9), weights = weights)
  
  fitted.vals <- data.frame(x = testdata[ ,x], E_val = E_val, 
                            lo = quants[1, ], hi = quants[2, ], xy_effect) %>%
    pivot_longer(cols = !c(x, E_val, lo, hi)) %>%
    mutate(weights = rep(weights, grid_size))
  
  plots[[i]] <- ggplot(data = fitted.vals, aes(x, value, group = name, color = weights)) +
    geom_line(alpha = 0.2) +
    scale_color_gradient(low='#d9d9d9', high='#545454', breaks = c(0.05, 0.15, 0.25)) +
    geom_ribbon(aes(ymin = lo, ymax = hi, group = NULL, color = NULL), 
                alpha = 0.5, fill = "#f34444") + 
    geom_line(aes(x, E_val), color = "#be4173", linetype = "dashed") +
    ggtitle(paste0("E(", colnames(data)[y], "|do(", do, "))")) +
    theme_minimal() +
    coord_cartesian(ylim = c(-4.5, 4), xlim = c(-1.8, 1.8)) +
    theme(legend.position = "bottom",
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          plot.title = element_text(hjust = 0.5, size = 8.5),
          plot.margin = margin(3,6,3,3)) 
  
  mylegend <- g_legend(plots[[i]])
  plots[[i]] <- plots[[i]] + theme(legend.position="none")
}

# Visualize all plots
plots[[n]] <- mylegend
lay <- rbind(matrix(1:(n-1), ncol = 6), c(NA, NA, n, n, NA, NA))
grid.arrange(grobs = plots, layout_matrix = lay)

# size = 3.33 | 5 x 9


# Local approximation
source("Local_approx/loc_fns.R")
source("Local_approx/loc_sampling_fns.R")

post.samples <- 800  # number of total posterior samples to generate
ndags <- length(fit$trace)
N <- nrow(data)
mcmc_mats <- fit$traceadd$incidence  # all sampled adjacency matrices 
fx_mats <- lapply(mcmc_mats,  # save non-zero causal effects for later
                  function(x) as(solve(diag(n) - x), "dgCMatrix"))

x <- 6  # intervention vab
plots <- list()  # to store the plots

for(i in 1:(n-1)) {
  y <- (1:n)[-x][i]
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
  par_lengths <- length(n.samples)

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
  E_val <- apply(xy_effect, 2, mean)
  quants <- apply(xy_effect, 2, quantile, probs = c(0.1, 0.9))
  fitted.vals <- data.frame(x = testdata[ ,x], t(xy_effect), E_val = E_val, 
                            lo = quants[1, ], hi = quants[2, ]) %>%
    pivot_longer(cols = !c(x, E_val, lo, hi))

  plots[[i]] <- ggplot(data = fitted.vals, aes(x, value, group = name)) +
    geom_line(alpha = 0.2, col = '#d4d4d4') +
    geom_line(aes(x, E_val), color = "#be4173", linetype = "dashed") +
    geom_ribbon(aes(ymin = lo, ymax = hi, group = NULL), 
                alpha = 0.5, fill = "#f34444") + 
    ggtitle(paste0("E(", colnames(data)[y], "|do(", colnames(data)[x], "))")) +
    coord_cartesian(ylim = c(-4.5, 4), xlim = c(-1.8, 1.8)) +
    theme_minimal() +
    theme(axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          plot.title = element_text(hjust = 0.5, size = 9),  # for small titles use 10
          plot.margin = margin(3,6,3,3)) 
  print(i)
}

lay <- rbind(matrix(1:(n-1), ncol = 6))
grid.arrange(grobs = plots, layout_matrix = lay)


# Joy plots
library(ggridges)

x <- 7  # intervention vab
y <- 11 # outcome vab

xy_effect <- table_effect(fit$fx, x, y)  # get matrix of effects
E_val <- apply(xy_effect, 1, weighted.mean, w = weights)  # expected value of effect
quants <- apply(xy_effect, 1, whdquantile, probs = c(0.1, 0.9), weights = weights)

joy_data <- data.frame(x = testdata[ ,x], E_val = E_val, 
                          lo = quants[1, ], hi = quants[2, ], xy_effect) %>%
  pivot_longer(cols = !c(x, E_val, lo, hi)) %>%
  mutate(weights = rep(weights, grid_size)) %>%
  filter(x %in% as.character(seq(-2, 2, length.out = 11)))

ggplot(joy_data, aes(x = value, y = fct_rev(as.factor(x)), height = after_stat(density), fill = stat(x))) +
  geom_density_ridges_gradient(aes(weight = weights), stat = "density", scale = 2.2, size = 0.1) +
  scale_fill_gradient2(low='#1196ee', mid = '#c9c9c9', high='#cf2020') +
  xlim(-3.6, 3.5) +
  xlab(paste0("E(", colnames(data)[y], "|do(", colnames(data)[x], "))")) +
  ylab(paste0("Intervention level of ", colnames(data)[x])) +
  theme_light() +
  scale_y_discrete(expand = expand_scale(add = c(0, 1))) +
  theme(legend.position = "none")

# size = 4.5 x 6.5 

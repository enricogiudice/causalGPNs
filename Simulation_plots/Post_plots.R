library(tidyverse)
library(BiDAG)
library(matrixStats)
library(MASS)
library(gridExtra)
library(gRbase)
library(questionr)
library(transport)

source("Fourier_fns.R")
source("BayesStanFns.R")
source("fxsampling_fns.R")
insertSource("GPscore.R", package = "BiDAG")

# Samples a seperate f for every sample
add.noise <- function(to_add, y, x, newX, sigma, rho) { 
  to_add <- as.matrix(to_add)
  newX <- as.matrix(newX)
  
  for(i in 1:nrow(to_add)) {
    to_add[i, ] <- to_add[i, ] + GP_gen(y, x, as.matrix(newX[i, ]), sigma, rho)
  }
  
  return(to_add)
}

# Generate DAG & data
set.seed(101)
n <- 4   # number of nodes
N <- 100  # number of samples
lambda <- 0
truegraph <- matrix(c(0,0,0,1,0,0,0,0,1,1,0,0,0,1,0,0), ncol = n, byrow = T)
Fou_result <- Fou_nldata(truegraph, N, lambda = lambda, noise.sd = 0.5, standardize = T) 
data <- Fou_result$data

# List all DAGs with n nodes
all.dags <- list()
adj <- matrix(0, nrow = n, ncol = n)
dag.counter <- 0
all.comb <- rep(list(c(0,1)), n*(n-1))
all.comb <- expand.grid(all.comb)  # all combinations outside of diagonal of adjacency matrix

for(i in 1:nrow(all.comb)) {
  adj[col(adj)!=row(adj)] <- as.numeric(all.comb[i, ])
  
  if(is.DAG(adj)) {
    dag.counter <- dag.counter + 1
    all.dags[[dag.counter]] <- adj
  }
}

# Compute posterior over DAGs (and save all scores)
true.post <- rep(NA, dag.counter)
allnames <- c("sets", "newscore", paste0("rho[", 1:(n-1), "]"), "mu", "sigma")
pars <- setNames(data.frame(matrix(NA, nrow = 0, ncol = n+3)), allnames)
parent.scores <- rep(list(pars), n)

for(k in 1:dag.counter) {
  dag <- all.dags[[k]]
  curr_score <- 0
  
  for(x in 1:n) {
    set <- parent.scores[[x]]$sets
    pax <- dag[,x]
    pax.str <- setNames(paste(pax, collapse = ""), "sets")  # parents of x
    check <- is.na(match(set, pax.str))  # check if parent set is already there
    
    if(all(check)) {  # new parent set
      
      if(sum(pax) == 0) {  # Compute local score
        loc_result <- Gauss.mcmc(data[,x])  # with MCMC
        loc_score <- loc_result[1]
      }
      else {
        loc_result <- GP.mcmc(data[ ,x], data[ ,which(pax==1)])  # with MCMC
        loc_score <- loc_result[1]
      }
      parent.scores[[x]][length(set)+1,c("sets", names(loc_result))] <- c(pax.str, loc_result)
    }
    
    else {  # fetch score from parent.scores
      ind <-  which(check == F)  # index of pax in set
      loc_score <- as.numeric(parent.scores[[x]]$newscore[ind])
    }
    curr_score <- curr_score + loc_score  # build score
  }
  true.post[k] <- curr_score
}
true.post <- true.post - logSumExp(true.post)  # normalize

# For every DAG, compute true posterior for E(Y|do(X = xval)) 
x <- 1
outcome <- 2
xval <- 0
true_fx <- vector()  # to store all the true sampled effects
totsamples <- 1e+04  # number of total samples
mcsamples <- table(sample(1:dag.counter, totsamples, replace = T, prob = exp(true.post)))

for(k in 1:length(mcsamples)) {
  adj <- all.dags[[as.numeric(names(mcsamples)[k])]]
  top_order <- rev(BiDAG:::DAGtopartition(n, adj)$permy)  # go down order
  curr_levels <- rep(list(matrix(0, ncol = 1, nrow = mcsamples[k])), n)  # list values at current nodes
  
  for(y in top_order) {
    
    if (y == x) {
      curr_levels[[y]] <- curr_levels[[y]] + xval  # set x to xval
      next
    }
    pax <- adj[, y]  # find parents of y
    pax.str <- paste(pax, collapse = "")
    parent.indexes <- which(parent.scores[[y]]$sets == pax.str) 
    hyperpars <- parent.scores[[y]][parent.indexes, ]  # find correct set in parent.scores
    hyperpars[-1] <- as.numeric(hyperpars[-1])
    
    mu <- 0  # hyperpars$mu
    sigma <- hyperpars$sigma
    
    if (sum(pax) == 1) {  # one parent
      parents <- which(pax == 1)
      rho <- hyperpars$`rho[1]`
      curr_levels[[y]] <- add.noise(curr_levels[[y]], y = data[ ,y], x = data[ ,parents], 
                                    newX = curr_levels[[parents]], sigma, rho)
    }
    
    else if (sum(pax) > 1) {  # More than one parent
      parents <- which(pax == 1)
      
      for(par in parents) {
        lscale <- paste0("rho[", match(par, parents), "]")  # rho in parents order
        rho <- as.numeric(hyperpars[lscale])
        curr_levels[[y]] <- add.noise(curr_levels[[y]], y = data[ ,y], x = data[ ,par],
                                      newX = curr_levels[[par]], sigma, rho)
      }
    }
    
    if (y == outcome) {
      true_fx <- append(true_fx, curr_levels[[y]])  # save causal fx
      break  # finish loop over top_order
    }
    
    else {
      curr_levels[[y]] <- curr_levels[[y]] + rnorm(mcsamples[k], mu, sigma)
    }
  }
}


# Loop over different samples
iters <- c(470, 550, 660, 810, 10e2, 13e2, 16e2, 21e2, 27e2, 36e2, 49e2, 67e2, 93e2, 13e3, 18e3)  
reps <- 50  # number of realizations to average over
results <- data.frame()

for(r in 1:reps) {
set.seed(101 + r)

for(i in 1:length(iters)) {
  
  # MC approach
  source("fxsampling_fns.R")
  glob_start <- Sys.time()
  GP.searchspace = set.searchspace(data, dual = F, "GP")
  fit <- GP.mcmc.fx(data, GP.searchspace, iterations = iters[i], burnin = 250/iters[i],
                    mcsamples = 1, truedag = NULL, stepsave = 2, x_levels = c(0, 1))
  glob_time <- Sys.time() - glob_start
  fit.weights <- exp(fit$weights)
  glob_effect <- table_effect(fit$fx, x, outcome)[1, ]

  wa_dis <- wasserstein1d(glob_effect, true_fx, wa = fit.weights)
  results <- rbind(results, data.frame(wass = wa_dis, iter = length(fit.weights),
                                       group = "MC approach",
                                       time = as.numeric(glob_time, units = "secs")))

  # Local approximation
  source("Local_approx/loc_fns.R")
  source("Local_approx/loc_sampling_fns.R")
  post.samples <- (iters[i] - 250)/2
  loc_start <- Sys.time()

  # Make list of sampled parents for each variable
  ndags <- length(fit$trace)
  mcmc_mats <- fit$traceadd$incidence  # all sampled adjacency matrices
  fx_mats <- lapply(mcmc_mats,  # save non-zero causal effects for later
                    function(x) as(solve(diag(n) - x), "dgCMatrix"))

  sets <- rep(NA, ndags)  # here we will store the parent sets as strings

  for(k in 1:ndags) {  # loop over all sampled dags

    if(fx_mats[[k]][x,outcome] == 0) {  # no causal effect
      pax <- "null"
    }

    else {
      pax <- mcmc_mats[[k]][,x]  # parents of x indicator
      pax <- paste(pax, collapse = "")  # as string
    }
    sets[k] <- pax
  }
  p.sets <- wtd.table(sets, weights = fit.weights)
  n.samples <- table(sample(names(p.sets), post.samples, replace = T, prob = p.sets))  # number of samples per parent set
  yD <- data[ ,outcome]  # data
  xy_effect <- NULL  # to store the samples from the posterior
  testdata <- apply(data, 2,
                    FUN = function(x) c(0, 1))

  for(pa in names(n.samples)) {

    if(pa == "null") {  # no causal effect
      smus <- rt(n.samples[pa], N-1) * sd(yD) / N + mean(yD)  # Estimate for E(y)
      xy_effect <- rbind(xy_effect, replicate(2, smus))
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
  loc_time <- Sys.time() - loc_start + fit$fit_time

  wal_dis <- wasserstein1d(xy_effect[ ,1], true_fx)
  results <- rbind(results, data.frame(wass = wal_dis, iter = nrow(xy_effect),
                                       group = "Local approximation",
                                       time = as.numeric(loc_time, units = "secs")))
  
  # True posteriors
  source("fxsampling_fns.R")
  post_fx <- vector()  # to store all the true sampled effects
  totsamples <- (iters[i] - 250)/2  # number of total samples
  mcsamples <- table(sample(1:dag.counter, totsamples, replace = T, prob = exp(true.post)))
  
  for(k in 1:length(mcsamples)) {
    adj <- all.dags[[as.numeric(names(mcsamples)[k])]]
    top_order <- rev(BiDAG:::DAGtopartition(n, adj)$permy)  # go down order
    curr_levels <- rep(list(matrix(0, ncol = 1, nrow = mcsamples[k])), n)  # list values at current nodes
    
    for(y in top_order) {
      
      if (y == x) {
        curr_levels[[y]] <- curr_levels[[y]] + xval  # set x to xval
        next
      }
      pax <- adj[, y]  # find parents of y
      pax.str <- paste(pax, collapse = "")
      parent.indexes <- which(parent.scores[[y]]$sets == pax.str) 
      hyperpars <- parent.scores[[y]][parent.indexes, ]  # find correct set in parent.scores
      hyperpars[-1] <- as.numeric(hyperpars[-1])
      
      mu <- 0  # hyperpars$mu
      sigma <- hyperpars$sigma
      
      if (sum(pax) == 1) {  # one parent
        parents <- which(pax == 1)
        rho <- hyperpars$`rho[1]`
        curr_levels[[y]] <- add.noise(curr_levels[[y]], y = data[ ,y], x = data[ ,parents], 
                                      newX = curr_levels[[parents]], sigma, rho)
      }
      
      else if (sum(pax) > 1) {  # More than one parent
        parents <- which(pax == 1)
        
        for(par in parents) {
          lscale <- paste0("rho[", match(par, parents), "]")  # rho in parents order
          rho <- as.numeric(hyperpars[lscale])
          curr_levels[[y]] <- add.noise(curr_levels[[y]], y = data[ ,y], x = data[ ,par],
                                        newX = curr_levels[[par]], sigma, rho)
        }
      }
      
      if (y == outcome) {
        post_fx <- append(post_fx, curr_levels[[y]])  # save causal fx
        break  # finish loop over top_order
      }
      
      else {
        curr_levels[[y]] <- curr_levels[[y]] + rnorm(mcsamples[k], mu, sigma)
      }
    }
  }
  
  wap_dis <- wasserstein1d(post_fx, true_fx)
  results <- rbind(results, data.frame(wass = wap_dis, iter = length(post_fx), 
                                          group = "True posterior",
                                          time = 0))
  }
}


results %>%  # average reps
  group_by(iter, group) %>%
  summarise(wass = mean(wass),
            time = mean(time)) %>%
  mutate(group = factor(group, 
                        levels = c("MC approach", "Local approximation", "True posterior"), 
                        ordered = T)) -> avg_results

# Plot of distances
color3 <- c('#db0000','#00acc7',"darkgreen")
ggplot(avg_results, aes(x = iter, y = wass, group = group)) +
  geom_line(aes(linetype = group)) +
  geom_point(aes(shape = group, color = group)) +
  scale_shape_manual(values = c(4, 2, 3)) +
  scale_color_manual(values = color3) +
  scale_linetype_manual(values = c("dashed","dotdash","dotted")) +
  scale_x_continuous(trans = 'log10', labels = function(x) format(x, scientific = TRUE)) +
  ylim(0.0122, 0.25821) +
  xlab("Samples") + ylab("Wasserstein distance") + 
  theme_light() +
  theme(legend.title = element_blank(), legend.position = "bottom",
        plot.margin = margin(5,14,5,10)) 
  
# size = 3.5 x 4.3

# Plot of times
avg_results %>%
  filter(group %in% c("MC approach", "Local approximation")) %>%
  ggplot(aes(x = time, y = wass, group = group)) +
    geom_point(aes(shape = group, color = group)) +
    scale_shape_manual(values = c(2, 4)) +
    scale_color_manual(values = color3[-3]) +
    scale_linetype_manual(values = c("dashed","dotdash")) +
    scale_x_continuous(trans = 'log10') +
    ylim(0.0122, 0.25821) +
    xlab("Time (s)") + ylab("Wasserstein distance") +
    theme_light() +
    theme(legend.title = element_blank(), legend.position = "bottom",
          plot.margin = margin(5,14,5,10)) +
    annotate(geom = "text", x = 43.23611, y = 0.25505229, label = "110", hjust = "center", color = "#858585", size = 3.1) +
    annotate(geom = "text", x = 170.84060, y = 0.22954721, label = "110", hjust = "center", color = "#858585", size = 3.1) +
    annotate(geom = "text", x = 124.13234, y = 0.20578470, label = "6375", hjust = "center", color = "#858585", size = 3.1) +
    annotate(geom = "text", x = 421.88783, y = 0.05969357, label = "6375", hjust = "center", color = "#858585", size = 3.1) 

# Density plots of last (18e3) estimates
coll_fx <- rbind(data.frame(fx = glob_effect, Method = "MC approach", weights = fit.weights),
                 data.frame(fx = xy_effect[ ,1], Method = "Local approximation", weights = 1/length(xy_effect[ ,1])),
                 data.frame(fx = true_fx, Method = "True posterior", weights = 1/length(true_fx)))
coll_fx$Method <- factor(coll_fx$Method, 
                         levels = c("MC approach", "Local approximation", "True posterior"), ordered = T)

bw <- bw.nrd0(coll_fx$fx[coll_fx$Method == "True posterior"])  # use same bandwidth
ggplot(coll_fx) +
  geom_density(aes(x = fx, fill = Method, color = Method, weight = weights), alpha = 0.11, bw = bw, adjust = 5) +
  scale_y_sqrt() +
  scale_fill_manual(values = c("red", "blue", "darkgreen")) +
  scale_color_manual(values = c("red", "blue", "darkgreen")) +
  theme_light() +
  theme(legend.title = element_blank(), legend.position = "bottom") +
  xlab(as.expression(bquote("E(X"[.(outcome)]*"|"*"do(X"[.(x)]*"=0"*")"*")"))) +
  ylab("Density") 

# size = 3.6 x 4.7

# Plots campring three densities
ggplot(coll_fx) +
  geom_density(aes(x = fx, fill = Method, color = Method, weight = weights), alpha = 0.11, bw = bw, adjust = 6) +
  scale_y_sqrt(limits = c(0, 1.2)) +
  scale_fill_manual(values = c("red", "blue", "darkgreen")) +
  scale_color_manual(values = c("red", "blue", "darkgreen")) +
  theme_light() +
  theme(legend.title = element_blank(), legend.position = "bottom",
        plot.title = element_text(hjust = 0.5)) +
  xlab(as.expression(bquote("E(X"[.(outcome)]*"|"*"do(X"[.(x)]*"=0"*")"*")"))) +
  ylab("Density") +
  ggtitle(paste0(length(post_fx), " Samples"))

ggarrange(plots[[1]], plots[[2]], plots[[3]], plots[[4]], ncol = 2, nrow = 2, common.legend = TRUE, legend = "bottom")

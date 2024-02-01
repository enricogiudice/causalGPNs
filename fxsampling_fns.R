set.searchspace <- function(data, dual, method, par = 1) {
  startspace <- NULL
  
  if(dual) {
    cor_mat <- cor(data)
    startspace <- dual_pc(cor_mat, nrow(data), alpha = 0.05, skeleton = T)
  }
  
  if(method == "GP") {
    n <- ncol(data)
    score <- scoreparameters("usr", data, 
                             usrpar = list(pctesttype = "bge", edgepmat = matrix(par, n, n)))
  }
  
  if(method == "bge") {
    score <- scoreparameters("bge", data, bgepar = list(am = par))
  }
  searchspace <- iterativeMCMC(scorepar = score, startspace = startspace, hardlimit = 14, 
                               verbose = F, scoreout = TRUE, alphainit = 0.01)
  
  list(score = score, scoretable = searchspace$scoretable, DAG = searchspace$DAG,
       maxorder = searchspace$maxorder, endspace = searchspace$endspace)
}


# Main function: For now hyperparameters are maximized; use model_gpcovf$sample directly for mcmc approach
GP.mcmc.fx <- function(data, searchspace, alpha = 0.05, order = FALSE, 
                       burnin = 0.2, iterations = 1000, mcsamples = 1, 
                       truedag = NULL, stepsave = 4, x_levels = seq(-2, 2, length.out = 101)) {
  n <- ncol(data)
  myScore <- searchspace$score
  if(is.null(truedag)) {
    start <- Sys.time()
  
    if(order) {
      fit <- orderMCMC(myScore, MAP = FALSE, chainout = TRUE, alpha = alpha, 
                      startorder = searchspace$maxorder, scoretable = searchspace$scoretable, 
                      startspace = searchspace$endspace, iterations = iterations, stepsave = stepsave)
    }
    
    else {
      fit <- partitionMCMC(myScore, alpha = alpha, startDAG = searchspace$DAG,
                          scoretable = searchspace$scoretable, iterations = iterations, stepsave = stepsave)
    }
    # Make list of sampled parents for each variable
    weights <- NULL  # store importance weights for all sampled DAGs
    allnames <- c("sets", "newscore", paste0("rho[", 1:(n-1), "]"), "mu", "sigma")
    pars <- setNames(data.frame(matrix(NA, nrow = 0, ncol = n+3)), allnames)
    parent.scores <- rep(list(pars), n)  # also store estimated (MAP) parameters of each variable given its parents
    parent.indx <- rep(NA, n)  # to store index of current parent set for each variable
    
    toburn <- round(burnin * fit$info$samplesteps)
    fit$traceadd$incidence <- fit$traceadd$incidence[-(1:toburn)]
    
    fit$trace <- fit$trace[-(1:toburn)]
    ndags <- length(fit$trace)
    fit$fx <- vector(mode = "list", length = ndags)
    fit$fit_time <- Sys.time() - start
  }
  
  else {  # known dag
    ndags <- iterations*0.8/stepsave
    fit <- list(traceadd = list(incidence = lapply(seq_len(ndags), function(x) truedag)),
                fx = vector(mode = "list", length = ndags),
                trace = rep(0, ndags))
    allnames <- c("sets", "newscore", paste0("rho[", 1:(n-1), "]"), "mu", "sigma")
    pars <- setNames(data.frame(matrix(NA, nrow = 0, ncol = n+3)), allnames)
    parent.scores <- rep(list(pars), n)  # also store estimated (MAP) parameters of each variable given its parents
    parent.indx <- rep(NA, n)  # to store index of current parent set for each variable
    weights <- NULL
  }
  
  for(k in 1:ndags) {
    dag <- fit$traceadd$incidence[[k]]
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
        parent.indx[x] <- length(set)+1  # save position of set for faster causal effect computation
      }
      
      else {  # fetch score from parent.scores
        ind <-  which(check == F)  # index of pax in set
        parent.indx[x] <- ind
        loc_score <- as.numeric(parent.scores[[x]]$newscore[ind])
      }
      curr_score <- curr_score + loc_score  # build score
    }
    # here now compute all causal fx
    fit$fx[[k]] <- EstimateAllFx(as.matrix(dag), parent.indx, parent.scores, data, mcsamples, x_levels)
    weights[k] <- curr_score - fit$trace[k]   # weight for current DAG 
  }
  
  if(is.null(truedag)) {
    fit$weights <- weights - logSumExp(weights)  # normalize weights
  }
  else {
    fit$weights <- rep(log(1/ndags), ndags)
  }
  return(fit)
}


bge.partition.mcmc <- function(searchspace, alpha = 0.05, 
                               order = FALSE, burnin = 0.2, iterations = 1000) {
  BGEScore <- searchspace$score
  
  if(order) {
    bge.fit <- orderMCMC(BGEScore, MAP = FALSE, chainout = TRUE, alpha = alpha, 
                         startorder = searchspace$maxorder, scoretable = searchspace$scoretable,
                         iterations = iterations, stepsave = 4)
  }
  
  else {
    bge.fit <- partitionMCMC(BGEScore, alpha = alpha, startDAG = searchspace$DAG, 
                             scoretable = searchspace$scoretable, 
                             iterations = iterations, stepsave = 4)
  }
  toburn <- round(burnin * bge.fit$info$samplesteps)
  bge.fit$traceadd$incidence <- bge.fit$traceadd$incidence[-(1:toburn)]
  
  return(bge.fit)
}

# New version
EstimateAllFx <- function(adj, parent.indexes, parent.scores, data, mcsamples, x_levels) {
  n <- ncol(adj)
  fxmat <- solve(diag(n) - adj)
  grid_size <- length(x_levels)
  true_fx <- array(0, dim = c(n, n, grid_size))  # to store the true effects
  top_order <- rev(BiDAG:::DAGtopartition(n, adj)$permy)  # go down order
  
  for(i in 1:n) {
    x <- top_order[i]
    curr_levels <- rep(list(matrix(0, nrow = grid_size, ncol = mcsamples)), n)  # list values at current nodes
    curr_levels[[x]] <- replicate(mcsamples, x_levels)
    
    for(y in top_order[-i]) {
      parents <- which(adj[, y] == 1)  # find parents of y
      lp <- length(parents)  # number of parents
      hyperpars <- parent.scores[[y]][parent.indexes[y], ]  # find correct set in parent.scores
      hyperpars[-1] <- as.numeric(hyperpars[-1])
      
      mu <- 0 # hyperpars$mu
      sigma <- hyperpars$sigma
      
      if (lp == 1) {  # one parent
        rho <- hyperpars$`rho[1]`
        curr_levels[[y]] <- curr_levels[[y]] + GP_gen(y = data[ ,y], x = data[ ,parents],
                                                      newX = curr_levels[[parents]], sigma, rho)
      }
      
      else if (lp > 1) {  # More than one parent

        for(par in parents) {
          lscale <- paste0("rho[", match(par, parents), "]")   # rho in parents order
          rho <- as.numeric(hyperpars[lscale])
          trans.par <- GP_gen(y = data[ ,y], x = data[ ,par],  # sample one function from posterior
                              newX = curr_levels[[par]], sigma, rho)
          curr_levels[[y]] <- curr_levels[[y]] + trans.par
        }
      }
      
      if(fxmat[x,y] != 0) {
        true_fx[x,y, ] <- rowMeans(curr_levels[[y]])  # save avg causal effect
      }
      curr_levels[[y]] <- apply(curr_levels[[y]], 2,  # add noise
                                function(f) f + rnorm(1, mu, sigma))  # mu zero works, but why?
    }
  }
  
  return(true_fx)
}


GP_gen <- function(y, x, newX, sigma, rho, mu = 0) {  # different f for every column of newX, faster
  y <- as.matrix(y)
  x <- as.matrix(x)
  #newX <- as.matrix(newX)
  K <- rbfkernel(x, rho^2)
  sampled_fs <- newX * 0  
  
  for(i in 1:ncol(newX)) {
    newx <- as.matrix(newX[ ,i])
    Ks <- rbfkernel(x, rho^2, newx)
    Kss <- rbfkernel(newx, rho^2)
    invK <- solve(K + sigma^2 * diag(nrow(x)))
    Mu <- mu + t(Ks) %*% invK %*% y
    Sigma <- Kss - t(Ks) %*% invK %*% Ks
    sampled_fs[ ,i] <- mvrnorm(1, Mu, Sigma)
  }
  
  return(sampled_fs)
}

# Turns fit$fx into a table for plotting the effect of x on y
table_effect <- function(fx, x, y) {
  ndags <- length(fx)
  grid_size <- dim(fx[[1]])[3]
  xy_effect <- matrix(NA, nrow = grid_size, ncol = ndags)
  
  for(dag in 1:length(fx)) {
    arr <- fx[[dag]]
    xy_effect[ ,dag] <- arr[x,y, ]
  }
  
  return(xy_effect)
}

# RBF kernel
`rbfkernel` <- function(X, sigma = 1, Y = NULL) {
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

# function to compute weighted variance
weighted.var <- function(x, w, na.rm = FALSE) {
  if (na.rm) {
    w <- w[i <- !is.na(x)]
    x <- x[i]
  }
  sum.w <- sum(w)
  sum.w2 <- sum(w^2)
  mean.w <- sum(x * w) / sum(w)
  (sum.w / (sum.w^2 - sum.w2)) * sum(w * (x - mean.w)^2, na.rm =
                                       na.rm)
}

# function for common legend
g_legend<-function(a.gplot){
  tmp <- ggplot_gtable(ggplot_build(a.gplot))
  leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
  legend <- tmp$grobs[[leg]]
  return(legend)}

# These are the functions to generate non-linear data with the Fourier transform approach
# Scale does not currently work for the true.fx? 
# Also the data should be scaled, otherwise we know the GP has trouble estimating. Or change the priors, max and min values
# Right now we recompute all the causal effects every time, even with two identical dags. Is there a faster way? It's very slo
# Mu, add it to the equations for the GP!
# Rho should be squared?? Check cov_exp_quad vs rbfkernel!

# Fourier transform data with known or random weights
Fou_trans <- function(x, lambda = NULL, sampled_weights = NULL, n_Four = 6) {
  
  n_Four <- min(n_Four, 10) # check we don't have too many
  n_Four <- max(n_Four, 2) # or too few!
  
  # store the multiples of x
  x_mults <- matrix(rep(x, each = n_Four)*(1:n_Four), nrow = n_Four)
  
  if(is.null(sampled_weights)){
    # compute the weights for the frequencies (we have most twice)
    if(lambda == 0) {
      freq_weights <- c(1, rep(0, n_Four, each = 2))
    }
    
    else {
      freq_weights <- exp(-rep(0:n_Four, each = 2)[-1] / lambda)
      freq_weights <- freq_weights/sum(freq_weights)
    }
    
    # Dirichlet sampler
    sampled_weights <- rgamma(length(freq_weights), freq_weights)
    sampled_weights <- sampled_weights/sum(sampled_weights)
    # chose random signs
    sampled_weights <- sampled_weights*sample(c(-1, 1), length(sampled_weights), replace = TRUE)
    # add beta
    sampled_weights <- runif(1, 0.5, 2)*sampled_weights
  }
  
  # transform x by taking cosine terms, weighted and ignoring intercept
  x_transform_cos <- colSums(sampled_weights[2*1:n_Four]*cos(x_mults))
  # transform x by taking sine terms, weighted and ignoring intercept
  x_transform_sin <- colSums(sampled_weights[2*1:n_Four + 1]*sin(x_mults))  
  # include linear term  
  x_transform <- sampled_weights[1]*x + x_transform_cos + x_transform_sin
  
  list("data" = x_transform, "sampled_weights" = sampled_weights)  
}


# Generate data with random weights
Fou_nldata <- function(dag, N, lambda, noise.sd, standardize = FALSE) {
  trueDAG <- 1*(dag != 0)  # the edge presence in the DAG
  n <- ncol(trueDAG)  # number of variables
  data <- matrix(0, nrow = N, ncol = n)  # to store the simulated data
  true.weights <- array(0, dim = c(n, n, 13))  # to store the true parameters
  scales <- rbind(rep(0, n), rep(1, n))  # to store scaling factors
  top_order <- rev(BiDAG:::DAGtopartition(n, trueDAG)$permy)  # go down order
  
  for (jj in top_order) {
    parents <- which(trueDAG[, jj] == 1)  # find parents
    lp <- length(parents)  # number of parents
    
    if (lp == 0) {  # zero parents
      data[ ,jj] <- rnorm(N, 0, 1)
    }
    
    else if (lp == 1) {  # one parent
      trans.pa <- Fou_trans(data[ ,parents], lambda)
      true.weights[parents, jj, ] <- trans.pa$sampled_weights
      data[ ,jj] <- trans.pa$data + rnorm(N, 0, noise.sd)
    }
    
    else {  # More than one parent
      trans.pa <- Fou_trans(data[ ,parents], lambda)
      true.weights[parents, jj, ] <- t(replicate(lp, trans.pa$sampled_weights))  
      data[, jj] <- rowSums(trans.pa$data) + rnorm(N, 0, noise.sd)
    }
    
    if(standardize) { 
      scd <- scale(data[, jj])
      data[, jj] <- scd
      scales[ ,jj] <- c(attr(scd,"scaled:center"), attr(scd,"scaled:scale"))
    }
  }
  list("data" = data, "true.weights" = true.weights, "scales" = scales)
}


# Compute causal effect between all variables
ComputeAllTrueFx <- function(adj, weights, scales, mcsamples, noise.sd, grid_size = 401) {
  n <- ncol(adj)
  fxmat <- solve(diag(n) - adj)
  true.fx <- array(0, dim = c(n, n, grid_size))  # to store the true effects
  top_order <- rev(BiDAG:::DAGtopartition(n, adj)$permy)  # go down order
  
  for(i in 1:n) {
    x <- top_order[i]
    x_levels <- seq(-2, 2, length.out = grid_size)
    curr_levels <- rep(list(matrix(0, nrow = grid_size, ncol = mcsamples)), n)  # list values at current nodes
    curr_levels[[x]] <- replicate(mcsamples, x_levels)
    
    for(y in top_order[-i]) {
      parents <- which(adj[, y] == 1)  # find parents of y
      lp <- length(parents)  # number of parents
      
      if (lp == 0) {  # zero parents
        curr_levels[[y]] <- curr_levels[[y]] + 
          replicate(mcsamples, rep(rnorm(1, 0, noise.sd), grid_size))
      }
      
      else if (lp == 1) {  # one parent
        curr_levels[[y]] <- curr_levels[[y]] + 
          Fou_trans(curr_levels[[parents]], sampled_weights = weights[parents,y, ])$data
      }
      
      else {  # More than one parent
        
        for(par in parents) {
          trans.par <- Fou_trans(curr_levels[[par]], sampled_weights = weights[par,y, ])$data
          curr_levels[[y]] <- curr_levels[[y]] + trans.par
        }
      }
      curr_levels[[y]] <- (curr_levels[[y]] - scales[1,y]) / scales[2,y] # scale
      
      if(fxmat[x,y] != 0) {
        true.fx[x,y, ] <- rowMeans(curr_levels[[y]])  # save avg causal effect
      }
      curr_levels[[y]] <- apply(curr_levels[[y]], 2,  # add noise
                                   function(f) f + rnorm(1, 0, noise.sd))
    }
  }
  
  return(true.fx)
}









# Deprecated
# ComputeAllTrueFx <- function(adj, weights, mcsamples, noise.sd) {  
#   n <- ncol(adj)
#   fxmat <- solve(diag(n) - adj)
#   true.fx <- array(0, dim = c(n, n, 401))  # to store the true effects
#   top_order <- rev(BiDAG:::DAGtopartition(n, adj)$permy)  # go down order
#   
#   for(i in 1:n) {
#     x <- top_order[i]
#     x_levels <- seq(-2, 2, length.out = 401)
#     
#     for(y in top_order[-(1:i)]) {
#       pars <- x
#       curr_levels <- vector(mode = "list", length = n)  # list values at current nodes
#       curr_levels[[x]] <- replicate(mcsamples, x_levels)  # of mcsamples x 401 matrices
#       reach_end <- FALSE
# 
#       while(reach_end == FALSE) {
#         
#         for(par in pars) {
#           chln <- which(adj[par, ] == 1, arr.ind = TRUE)  # find children of parents 
#           act_chln <- chln[fxmat[chln,y] != 0]  # children that have an effect on y 
#         
#           for(chld in act_chln) {
#             curr_levels[[chld]] <- curr_levels[[chld]] + Fou_trans(curr_levels[[par]], 
#                                           sampled_weights = weights[par,chld, ])$data
#           }
#             if(chld != y) {  # add noise
#               curr_levels[[chld]] <- apply(curr_levels[[chld]], 2, 
#                                            function(f) f + rnorm(1, 0, noise.sd))
#             }
#             
#             else {  # add path to the total effect
#               true.fx[x,y, ] <- true.fx[x,y, ] + rowMeans(curr_levels[[chld]])
#             }
#           }
#         finished <- act_chln == y  # which paths have ended
#         reach_end <- all(finished)  # check if all paths have ended
#         
#         if(sum(finished) > 0) {  # move to next step and delete paths that ended
#           pars <- act_chln[-which(finished)]
#         }
#         
#         else {
#           pars <- act_chln
#         }
#       }
#     }
#   }
#   
#   return(true.fx)
# }



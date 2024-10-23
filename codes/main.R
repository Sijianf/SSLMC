################################# Logs
# 2024.10.04: Added initialization based on SVD
# 2024.10.04: Developed from the adjusted version 4 for the BiSSLMC: 
#             SSL + DRIMC + confidence level
#             Updated using gradient  
#             Updated column-wise




mysplit = function(vec){
  
  checking_list = cumsum(vec!=0)
  counts = table(checking_list)
  non_zeros_cuts = as.integer(names(counts[counts>1]))
  non_zeros_cuts = non_zeros_cuts[non_zeros_cuts>0]
  out_list = c()
  
  c1=vec[checking_list>=1 & checking_list<non_zeros_cuts[1]]
  c1=c(c1,vec[checking_list==non_zeros_cuts[1]][1])
  out_list[[1]] = c1
  
  if(length(non_zeros_cuts)>1){
    for(i in 2:length(non_zeros_cuts)){
      this_cluster=vec[checking_list>non_zeros_cuts[i-1] & checking_list<non_zeros_cuts[i]]
      this_cluster=c(this_cluster,vec[checking_list==non_zeros_cuts[i]][1])
      out_list[[i]] = this_cluster
    }
  }
  
  if(length(non_zeros_cuts)==0){
    out_list = list(0)
  }
  
  return(out_list)
  
}



################################# Get the logit link
get_logit <- function(x){
  if(is.list(x)){
    lapply(x, function(z) 1/(1+exp(-z)))
  }else{
    1/(1+exp(-x))
  }
}



################################# Get the Laplace density
get_Laplace <- function(x, lambda){

  lambda/2*exp(-lambda*abs(x))

}



################################# Get the lambdastar
get_lambdastar <- function(x, thetas, lambdas) {

  if(exists("get_lambdastar_matrix", mode = "function") == FALSE){ 
    # this is the original function when Rcpp is not loaded. 
    pstar0 = (1-thetas)*lambdas[1]/2*exp(-lambdas[1]*abs(x))
    pstar1 = thetas*lambdas[2]/2*exp(-lambdas[2]*abs(x))
    pstar  = pstar1/(pstar0+pstar1)
    lambdastar = lambdas[1]*(1-pstar)+lambdas[2]*pstar
    return(lambdastar)

  } else { 
    # this is the Rcpp function with available for either double or matrix as input. 
    if(is.vector(x)) x = t(as.matrix(x))
    if (is.numeric(x) && length(x) == 1) {
      return(get_lambdastar_double(x, thetas, lambdas))
    } else if (is.matrix(x) | is.vector(x)) {
      if (!is.numeric(thetas) || length(thetas) != ncol(x)) {
        stop("For matrix input, thetas must be a numeric vector with length equal to the number of columns of the matrix.")
      }
      return(get_lambdastar_matrix(as.matrix(x), thetas, lambdas))
    } else {
      stop("Input x must be either a single numeric value or a matrix.")
    }
  }
  
}



################################# Get the soft thresholding
# soft_thresholding <- function(x, lambdastar, eta){
# 
#   if(x > lambdastar*eta){
#     return(x - lambdastar*eta)
#   }
#   else{
#     if(x < -lambdastar*eta){
#       return(x + lambdastar*eta)
#     }
#     else{
#       return(0)
#     }
#   }
# 
# }



################################# Get the proximal
# get_proximal <- function(Y, xi, A, B, U, V, vec_thetas, lambdas, eta){
# 
#   I = nrow(Y)
#   J = ncol(Y)
#   K = ncol(A)
#   proximal_A = A
#   proximal_B = B 
# 
#   H = (1+xi*Y-Y)/(1+exp(-U%*%A%*%t(B)%*%t(V)))
# 
#   dA = - xi*t(U)%*%Y%*%V%*%B + t(U)%*%H%*%V%*%B
#   dB = - xi*t(V)%*%t(Y)%*%U%*%A + t(V)%*%t(H)%*%U%*%A
#   
#   for(k in 1:K){
# 
#     for(i in 1:I){
#       lambdastar = get_lambdastar(x = A[i,k], thetas = vec_thetas[k], lambdas = lambdas)
#       proximal_A[i,k] = soft_thresholding(x = A[i,k] - eta*dA[i,k], lambdastar, eta)
#     }
# 
#     for(j in 1:J){
#       lambdastar = get_lambdastar(x = B[j,k], thetas = vec_thetas[k], lambdas = lambdas)
#       proximal_B[j,k] = soft_thresholding(x = B[j,k] - eta*dB[j,k], lambdastar, eta)
#     }
# 
#   }
# 
#   return(list(proximal_A = proximal_A,
#               proximal_B = proximal_B))
# 
# }



################################# Update the stick-breaking fractions
get_thetas <- function(mat, a, b, tol = 0){

  K = ncol(mat)
  N = nrow(mat)
  thetas = rep(0,K)
  counts = rep(0,K)

  for (k in 1:K){

    q = sum(abs(mat[,k]) > tol, na.rm = TRUE)
    counts[k] <- q
    thetas[k] <- (a+q) / (a+b+N)
    
  }

  return(list(thetas = thetas, counts = counts))

}



################################# get the log likelihood
get_logLikelihood <- function(Y, xi, A, B, tilde_thetas, thetas, tilde_lambdas, lambdas, tilde_alpha, tilde_beta, alpha, beta){

  I = nrow(Y)
  J = ncol(Y)

  M = A%*%t(B)

  out = sum(xi*Y * M) - sum((1 + xi*Y - Y) * log(1+exp(M))) - 
    sum(get_lambdastar(x = A, thetas = tilde_thetas, lambdas = tilde_lambdas) * A) - 
    sum(get_lambdastar(x = B, thetas = thetas, lambdas = lambdas) * B) + 
    sum((tilde_alpha - 1) * log(tilde_thetas) + (tilde_beta - 1) * log(1 - tilde_thetas)) + 
    sum((alpha - 1) * log(thetas) + (beta - 1) * log(1 - thetas))

  return(out)

}



################################# Update AB
SSLMC <- function(Y, A = NULL, B = NULL, xi = 1, 
                  tilde_lambda_0 = 5, tilde_lambda_1 = 1,
                  lambda_0 = 5, lambda_1 = 1, method = "svd",
                  tilde_alpha = 1, tilde_beta = 1, alpha = 1, beta = 1, 
                  K_init = 50, thisSeed = 123, eta = 0.001, IBP = 1, 
                  tol = 1e-10, max_iter = 500, show_plot = FALSE) {
  
  I <- nrow(Y)
  J <- ncol(Y)
  K <- K_init

  tilde_lambda_0 <- tilde_lambda_0       # Matrix A's spike parameter
  tilde_lambda_1 <- tilde_lambda_1       # Matrix A's slab parameter
  lambda_0 <- lambda_0                   # Matrix B's spike parameter
  lambda_1 <- lambda_1                   # Matrix B's slab parameter
  tilde_lambdas <- c(tilde_lambda_0, tilde_lambda_1)
  lambdas <- c(lambda_0, lambda_1)

  A_in <- A
  B_in <- B

  if(method == "svd") {
    Y_svd = svd(Y)
    if(length(Y_svd$d)<K_init){
      A <- Y_svd$u%*%diag(sqrt(Y_svd$d))
      B <- Y_svd$v%*%diag(sqrt(Y_svd$d))
    } else{
      A <- Y_svd$u[,1:K_init]%*%diag(sqrt(Y_svd$d[1:K_init]))
      B <- Y_svd$v[,1:K_init]%*%diag(sqrt(Y_svd$d[1:K_init]))
    }
    if(!is.null(A_in)) A <- A_in
    if(!is.null(B_in)) B <- B_in
  } else {
    A <- matrix(rexp(I * K_init, rate = 1), nrow = I, ncol = K_init)
    B <- matrix(rexp(J * K_init, rate = 1), nrow = J, ncol = K_init)
    if(!is.null(A_in)) A <- A_in
    if(!is.null(B_in)) B <- B_in
  }

  K = K_init = min(K_init, ncol(A))
  tilde_nus <- sort(rbeta(K_init, tilde_alpha, tilde_beta), decreasing = TRUE)   # K x 1 vector
  tilde_thetas <- rep(0.5,K_init)                                                # K x 1 vector
  nus <- sort(rbeta(K_init, alpha, beta), decreasing = TRUE)                     # K x 1 vector
  thetas <- rep(0.5,K_init)                                                      # K x 1 vector
  
  # monitor the log likelihood
  logLikelihood_old <- get_logLikelihood(
    Y = Y, 
    xi = xi, 
    A = A, 
    B = B, 
    tilde_thetas = tilde_thetas, 
    thetas = thetas, 
    tilde_lambdas = tilde_thetas, 
    lambdas = lambdas, 
    tilde_alpha = tilde_alpha, 
    tilde_beta = tilde_beta, 
    alpha = alpha, 
    beta = beta)

  delta_logLikelihood_old <- 1000
  delta_logLikelihood_save = c()
  logLikelihood_save = c()

  # main loop
  A_lag <- A
  B_lag <- B

  for (i in 1:max_iter) {
    cat("Iteration: ", i, "\n")

    if(i == 1){
      A <- A
      B <- B
    } else {
      A_momentum <- A+(i-2)/(i+1)*(A-A_lag)
      B_momentum <- B+(i-2)/(i+1)*(B-B_lag)
      A_lag <- A
      B_lag <- B
      
      proximals <- get_proximal(
          Y = Y,
          xi = xi,
          eta = eta,
          A = A_momentum,
          B = B_momentum,
          tilde_lambdas = tilde_lambdas,
          lambdas = lambdas,
          tilde_thetas = tilde_thetas,
          thetas = thetas
          )
      # proximal_A
      A <- proximals$proximal_A      

      # proximal_B
      B <- proximals$proximal_B

    }
    
    # Update sparsity
    A_sparse <- get_thetas(A, tilde_alpha, tilde_beta, tol)
    B_sparse <- get_thetas(B, alpha, beta, tol)
    counts = A_sparse$counts+B_sparse$counts
    tilde_thetas <- A_sparse$thetas
    thetas <- B_sparse$thetas
    cat("Thetas update finished.","\n")
    
    # Get the change of probability with one cluster leave out
    overall_prob = get_logit(A%*%t(B))
    bicluster_prob = unlist(
      lapply(1:K, FUN = function(k){
        max(abs(get_logit(A[,-k]%*%t(B[,-k]))-overall_prob))
      }))
    
    # Re-order the columns: 
    if(IBP == 1){
      re_order = order(counts, decreasing = TRUE)
      A        = A[    ,re_order]           
      B        = B[    ,re_order]
      A_lag    = A_lag[,re_order]           
      B_lag    = B_lag[,re_order]
      tilde_thetas = tilde_thetas[re_order]
      thetas = thetas[re_order]
      cat("AB update finished.","\n")
      
    } else {
      re_order = order(bicluster_prob, decreasing = TRUE)
      A        = A[    ,re_order]           
      B        = B[    ,re_order]
      A_lag    = A_lag[,re_order]           
      B_lag    = B_lag[,re_order]
      tilde_thetas = tilde_thetas[re_order]
      thetas = thetas[re_order]
      cat("AB update finished.","\n")
    }
    
    # remove all zero columns 
    # keep = bicluster_prob >= 0.001*(median(bicluster_prob))
    # keep = counts <= (I+J)*1
    keep = colSums(A==0) < I & colSums(B==0) < J
    # cat(sort(bicluster_prob,decreasing = TRUE), "\n")
    K = sum(keep)
    cat("K =", K, "\n")
    if(K <= 2){
      break
    }
    A            = A[    ,keep]
    B            = B[    ,keep]
    A_lag        = A_lag[,keep]
    B_lag        = B_lag[,keep]
    tilde_thetas = tilde_thetas[keep]
    thetas       = thetas[keep]

    if(show_plot == TRUE){
      library(cowplot)
      plt_A = plot_matrix(A, main = glue("Estimated A, Columns = {K}"), legend = FALSE)
      plt_B = plot_matrix(B, main = glue("Estimated B, Columns = {K}"), legend = FALSE)
      plt_AB = plot_matrix(get_logit(tcrossprod(A,B))>0.5, main = glue("Estimated Probability"), legend = FALSE)
      plt_Y = plot_matrix(Y, main = glue("Observed Matrix"), legend = TRUE)
      legend = get_legend(plt_Y)
      gridExtra::grid.arrange(
        arrangeGrob(plt_A, plt_B, plt_AB, plt_Y + theme(legend.position = "none"), ncol = 2),
        legend = legend,
        ncol = 2,
        widths = c(3, 0.5)  # Adjust widths to allocate space for the legend)
      )
    }

    # stop earlier by log likelihood 
    logLikelihood <- get_logLikelihood(
      Y = Y, 
      xi = xi, 
      A = A, 
      B = B, 
      tilde_thetas = tilde_thetas, 
      thetas = thetas, 
      tilde_lambdas = tilde_thetas, 
      lambdas = lambdas, 
      tilde_alpha = tilde_alpha, 
      tilde_beta = tilde_beta, 
      alpha = alpha, 
      beta = beta)
    
    if(is.finite(abs(logLikelihood - logLikelihood_old) / abs(logLikelihood_old))){
      delta_logLikelihood <- abs(logLikelihood - logLikelihood_old) / abs(logLikelihood_old)
      
      if (delta_logLikelihood < tol) {
        break
      }
      # if ((i > 100) & (delta_logLikelihood > delta_logLikelihood_old)) {
      #   break
      # }
    }
    
    logLikelihood_old <- logLikelihood
    delta_logLikelihood_old <- delta_logLikelihood
    
    delta_logLikelihood_save = c(delta_logLikelihood_save, delta_logLikelihood_old)
    logLikelihood_save = c(logLikelihood_save, logLikelihood_old)
    
  }
  
  out <- list(A = A, B = B, tilde_thetas = tilde_thetas, thetas = thetas, counts = sum(counts))
  
  if(show_plot == TRUE) plot(logLikelihood_save)

  return(out)
}



SSLMC_ladder <- function(Y, U = NULL, V = NULL, xi = 1, 
                         tilde_lambda0s, tilde_lambda1,
                         lambda0s, lambda1,
                         tilde_alpha = 1, tilde_beta = 1, alpha = 1, beta = 1, 
                         K_init = 50, thisSeed = 123, eta = 0.001, IBP = 1, 
                         tol = 1e-10, max_iter = 100, show_plot = FALSE){
  
  L = length(lambda0s)
  update_tilde_lambda0 = 1
  counts = rep(0,L)
  tilde_lambda0 = tilde_lambda0s[1]
  lambda0 = lambda0s[1]
  
  for (l in 1:L) {
    print(glue::glue('Ladder {l}: 
    tilde_lambdas = ({tilde_lambda0},{tilde_lambda1}) 
    lambdas = ({lambda0},{lambda1}) \n'))
    
    if(l == 1) {
      A = NULL
      B = NULL 
      K_init = K_init
      this_SSLMC <- SSLMC(A = A, B = B, 
                          Y = Y, U = U, V = V, xi = xi, 
                          tilde_lambdas = c(tilde_lambda0, tilde_lambda1),
                          lambdas = c(lambda0, lambda1),
                          tilde_alpha = tilde_alpha, tilde_beta = tilde_beta, alpha = alpha, beta = beta, 
                          K_init = K_init,
                          thisSeed = thisSeed, eta = eta, tol = tol, IBP = IBP, 
                          max_iter = max_iter, show_plot = show_plot)
      counts[1] = this_SSLMC$counts
      
    } else {
      A = this_SSLMC$A
      B = this_SSLMC$B
      K_init = ncol(A)
      if(update_tilde_lambda0){
        tilde_lambda0 = tilde_lambda0s[l]
      }
      lambda0 = lambda0s[l]
      
      this_SSLMC <- SSLMC(A = A, B = B, 
                          Y = Y, U = U, V = V, xi = xi, 
                          tilde_lambdas = c(tilde_lambda0, tilde_lambda1),
                          lambdas = c(lambda0, lambda1),
                          tilde_alpha = tilde_alpha, tilde_beta = tilde_beta, alpha = alpha, beta = beta, 
                          K_init = K_init,
                          thisSeed = thisSeed, eta = eta, tol = tol, IBP = IBP, 
                          max_iter = max_iter, show_plot = show_plot)
      counts[l] = this_SSLMC$counts
      if(counts[l]>counts[l-1]){
        tilde_lambda0 = tilde_lambda0s[l-1]
        update_tilde_lambda0 = 0
      }
    }
  }
  
  out <- list(A = this_SSLMC$A, B = this_SSLMC$B, tilde_thetas = this_SSLMC$tilde_thetas, thetas = this_SSLMC$thetas)
  return(out)
}



# This is a simulation setting adjusted from Lee_Huang_2014_A_biclustering algorithm for binary matrices based on penalized Bernoulli likelihood.
my_levelplot <- function(mat, 
                         main="Binary Matrix Plot", 
                         xlab = "X", 
                         ylab = "Y", 
                         col.regions=c("white", "black"), 
                         colorkey = TRUE, 
                         aspect="fill", 
                         cuts = 1, 
                         at=c(0, 0.5, 1)){

  levelplot(t(apply(mat, 2, rev)),
      main = main,
      xlab = xlab,
      ylab = ylab,
      col.regions = col.regions,
      colorkey = colorkey,
      aspect=aspect,
      cuts=cuts,
      at=at)

}


sim_data <- function(I=120,J=300,p=0.4,row1=0.3,col1=0.4,p1=0.95,row0=0.4,col0=0.5,p0=0.05,seed,...){
  require(lattice)
  set.seed(seed)
  p = p 
  I = I
  J = J
  mat = matrix(rbinom(I*J,size=1,prob=p), nrow = I, ncol = J)
  p1 = p1
  I1 = as.integer(I*row1)
  J1 = as.integer(J*col1)
  mat1 = matrix(rbinom(I1*J1,size=1,prob=p1), nrow = I1, ncol = J1)
  p0 = p0
  I0 = as.integer(I*row0)
  J0 = as.integer(J*col0)
  mat0 = matrix(rbinom(I0*J0,size=1,prob=p0), nrow = I0, ncol = J0)
  
  mat[(I-I1+1):I,1:J1] = mat1
  mat[1:I0,(J-J0+1):J] = mat0

  my_levelplot(mat,...)
  return(mat)
  
}


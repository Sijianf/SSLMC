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
get_logLikelihood <- function(Q, A, B, tilde_thetas, thetas, tilde_lambdas, lambdas, tilde_alpha, tilde_beta, alpha, beta){

  I = nrow(Q)
  J = ncol(Q)

  M = A%*%t(B)

  out = sum(log(get_logit(Q*M))) - 
    sum(get_lambdastar(x = A, thetas = tilde_thetas, lambdas = tilde_lambdas) * A) - 
    sum(get_lambdastar(x = B, thetas = thetas, lambdas = lambdas) * B) + 
    sum((tilde_alpha - 1) * log(tilde_thetas) + (tilde_beta - 1) * log(1 - tilde_thetas)) + 
    sum((alpha - 1) * log(thetas) + (beta - 1) * log(1 - thetas))

  out = -2*sum(log(get_logit(Q*M)))+(log(I*J)*(1+sum(A!=0)+sum(B!=0)))
  
  return(out)
}



################################# Update AB
SSLMC <- function(Y, A = NULL, B = NULL, xi = 2, 
                  tilde_lambda_0 = 5, tilde_lambda_1 = 1,
                  lambda_0 = 5, lambda_1 = 1, method = "svd",
                  tilde_alpha = 1, tilde_beta = 1, alpha = 1, beta = 1, 
                  K_init = 50, thisSeed = 123, eta = 0.001, IBP = 1, 
                  tol = 1e-10, max_iter = 500, show_plot = FALSE) {
  
  I <- nrow(Y)
  J <- ncol(Y)
  K <- K_init

  Q <- as.matrix(xi*Y-1)
  Qc <- Q - sum(Q)/(I*J)

  tilde_lambda_0 <- tilde_lambda_0       # Matrix A's spike parameter
  tilde_lambda_1 <- tilde_lambda_1       # Matrix A's slab parameter
  lambda_0 <- lambda_0                   # Matrix B's spike parameter
  lambda_1 <- lambda_1                   # Matrix B's slab parameter
  tilde_lambdas <- c(tilde_lambda_0, tilde_lambda_1)
  lambdas <- c(lambda_0, lambda_1)

  A_in <- A
  B_in <- B

  if(method == "svd") {
    Q_svd = svd(Qc)
    if(length(Q_svd$d)<K_init){
      A <- Q_svd$u%*%diag(sqrt(Q_svd$d))
      B <- Q_svd$v%*%diag(sqrt(Q_svd$d))
      # A <- Q_svd$u
      # B <- Q_svd$v
    } else {
      A <- Q_svd$u[,1:K_init]%*%diag(sqrt(Q_svd$d[1:K_init]))
      B <- Q_svd$v[,1:K_init]%*%diag(sqrt(Q_svd$d[1:K_init]))
      # A <- Q_svd$u[,1:K_init]
      # B <- Q_svd$v[,1:K_init]
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
    Q = Q, 
    A = A, 
    B = B, 
    tilde_thetas = tilde_thetas, 
    thetas = thetas, 
    tilde_lambdas = tilde_lambdas, 
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
    # cat("Iteration: ", i, "\n")

    if(i == 1){
      A <- A
      B <- B
    } else {
      A_momentum <- A+(i-2)/(i+1)*(A-A_lag)
      B_momentum <- B+(i-2)/(i+1)*(B-B_lag)
      A_lag <- A
      B_lag <- B
      
      proximals <- get_proximal(
          Q = Q,
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
    # cat("Thetas update finished.","\n")
    
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
      # cat("AB update finished.","\n")
      
    } else {
      re_order = order(bicluster_prob, decreasing = TRUE)
      A        = A[    ,re_order]           
      B        = B[    ,re_order]
      A_lag    = A_lag[,re_order]           
      B_lag    = B_lag[,re_order]
      tilde_thetas = tilde_thetas[re_order]
      thetas = thetas[re_order]
      # cat("AB update finished.","\n")
    }
    
    # remove all zero columns 
    # keep = bicluster_prob >= 0.001*(median(bicluster_prob))
    # keep = counts <= (I+J)*1
    keep = colSums(A==0) < I & colSums(B==0) < J
    # cat(sort(bicluster_prob,decreasing = TRUE), "\n")
    K = sum(keep)
    # cat("K =", K, "\n")
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

    # Stop earlier by log likelihood 
    if(i > 1){
      
      logLikelihood <- get_logLikelihood(
        Q = Q, 
        A = A, 
        B = B, 
        tilde_thetas = tilde_thetas, 
        thetas = thetas, 
        tilde_lambdas = tilde_lambdas, 
        lambdas = lambdas, 
        tilde_alpha = tilde_alpha, 
        tilde_beta = tilde_beta, 
        alpha = alpha, 
        beta = beta)

      if(is.finite(abs(logLikelihood - logLikelihood_old) / abs(logLikelihood_old))){
        delta_logLikelihood <- abs(logLikelihood - logLikelihood_old) / abs(logLikelihood_old)

        logLikelihood_old <- logLikelihood
        delta_logLikelihood_old <- delta_logLikelihood
      
        delta_logLikelihood_save = c(delta_logLikelihood_save, delta_logLikelihood_old)
        logLikelihood_save = c(logLikelihood_save, logLikelihood_old)

        if (delta_logLikelihood < tol) break
      }
      
    } # End the earlier termination
    
  } # End of the main iterations 
  
  out <- list(A = A, B = B, tilde_thetas = tilde_thetas, thetas = thetas, counts = sum(counts), BIC = min(logLikelihood_save))
  
  plot(1:length(logLikelihood_save),logLikelihood_save,type = ifelse(length(logLikelihood_save)==1,"o","l"),
        # xlim = c(0,max_iter),
        xlab = "Iterations",ylab = "BIC")

  return(out)
}


SSLMC_ladder <- function(Y, xi = 2, A = NULL, B = NULL, 
                         tilde_lambda_0s, tilde_lambda_1,
                         lambda_0s, lambda_1, method = "svd",
                         tilde_alpha = 1, tilde_beta = 1, alpha = 1, beta = 1, 
                         K_init = 30, thisSeed = 123, eta = 0.001, IBP = 1, 
                         tol = 1e-10, max_iter = 100, show_plot = FALSE){
  
  L = length(lambda_0s)
  counts = rep(0,L)
  BIC = rep(0,L)
  update_tilde_lambda_0 = 1
  update_lambda_0 = 1
  tilde_lambda_0 = tilde_lambda_0s[1]
  lambda_0 = lambda_0s[1]
  
  for (l in 1:L) {
    print(glue::glue('Ladder {l}: 
    tilde_lambdas = ({tilde_lambda_0},{tilde_lambda_1}) 
    lambdas = ({lambda_0},{lambda_1}) \n'))
    
    if(l == 1) {
      A = A
      B = B
      K_init = K_init
      this_SSLMC <- SSLMC(A = A, B = B, 
                          Y = Y, xi = xi, 
                          tilde_lambda_0 = tilde_lambda_0, tilde_lambda_1 = tilde_lambda_1,
                          lambda_0 = lambda_0, lambda_1 = lambda_1,
                          tilde_alpha = tilde_alpha, tilde_beta = tilde_beta, alpha = alpha, beta = beta, 
                          K_init = K_init, method = method,
                          thisSeed = thisSeed, eta = eta, tol = tol, IBP = IBP, 
                          max_iter = max_iter, show_plot = show_plot)
      counts[1] = this_SSLMC$counts
      BIC[1] = this_SSLMC$BIC
      print(glue::glue('K = {ncol(this_SSLMC$A)} and counts = {this_SSLMC$counts} \n'))
      
    } else {
      A = this_SSLMC$A
      B = this_SSLMC$B
      K_init = ncol(A)
      # A = A
      # B = B
      # K_init = K_init
      if(update_tilde_lambda_0 == 1){
        tilde_lambda_0 = tilde_lambda_0s[l]
      }
      if(update_lambda_0 == 1){
        lambda_0 = lambda_0s[l]
      }
      
      this_SSLMC <- SSLMC(A = A, B = B, 
                          Y = Y, xi = xi, 
                          tilde_lambda_0 = tilde_lambda_0, tilde_lambda_1 = tilde_lambda_1,
                          lambda_0 = lambda_0, lambda_1 = lambda_1,
                          tilde_alpha = tilde_alpha, tilde_beta = tilde_beta, alpha = alpha, beta = beta, 
                          K_init = K_init, method = method,
                          thisSeed = thisSeed, eta = eta, tol = tol, IBP = IBP, 
                          max_iter = max_iter, show_plot = show_plot)
      counts[l] = this_SSLMC$counts
      BIC[l] = this_SSLMC$BIC
      if(counts[l]>counts[l-1]){
        tilde_lambda_0 = tilde_lambda_0s[l-1]
        update_tilde_lambda_0 = 0
        lambda_0 = lambda_0s[l-1]
        update_lambda_0 = 0
      }

      print(glue::glue('K = {ncol(this_SSLMC$A)} and counts = {this_SSLMC$counts} \n'))

    }
  }
  
  out <- list(A = this_SSLMC$A, B = this_SSLMC$B, tilde_thetas = this_SSLMC$tilde_thetas, thetas = this_SSLMC$thetas, BIC = BIC)
  return(out)
}



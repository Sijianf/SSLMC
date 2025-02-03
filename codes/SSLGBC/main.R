################################# Logs
# 2025.01.18: Changed as SSLGBC only for biclusterings. 
# 2024.12.30: Changed from transferred Q to original Y matrix.
# 2024.12.07: Changed to Rcpp with more functions
#             Corrected the mis-match of the centering step
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
get_lambdastarR <- function(x, rhos, lambdas) {

  if(exists("get_lambdastar_matrix", mode = "function") == FALSE){ 
    # this is the original function when Rcpp is not loaded. 
    pstar0 = (1-rhos)*lambdas[1]/2*exp(-lambdas[1]*abs(x))
    pstar1 = rhos*lambdas[2]/2*exp(-lambdas[2]*abs(x))
    pstar  = pstar1/(pstar0+pstar1)
    lambdastar = lambdas[1]*(1-pstar)+lambdas[2]*pstar
    return(lambdastar)

  } else { 
    # this is the Rcpp function with available for either double or matrix as input. 
    if(is.vector(x)) x = t(as.matrix(x))
    if (is.numeric(x) && length(x) == 1) {
      return(get_lambdastar_double(x, rhos, lambdas))
    } else if (is.matrix(x) | is.vector(x)) {
      if (!is.numeric(rhos) || length(rhos) != ncol(x)) {
        stop("For matrix input, rhos must be a numeric vector with length equal to the number of columns of the matrix.")
      }
      return(get_lambdastar_matrix(as.matrix(x), rhos, lambdas))
    } else {
      stop("Input x must be either a single numeric value or a matrix.")
    }
  }
  
}


################################# Update the stick-breaking fractions
get_rhosR <- function(mat, a, b, tol = 0){

  K = ncol(mat)
  N = nrow(mat)
  rhos = rep(0,K)
  counts = rep(0,K)

  for (k in 1:K){

    q = sum(abs(mat[,k]) > tol, na.rm = TRUE)
    counts[k] <- q
    rhos[k] <- (a+q) / (a+b+N)
    
  }

  return(list(rhos = rhos, counts = counts))

}


################################# get the log likelihood
get_logLikelihoodR <- function(Y, A, B, tilde_rhos, rhos, tilde_lambdas, lambdas, tilde_alpha, tilde_beta, alpha, beta){

  I = nrow(Y)
  J = ncol(Y)

  M = A%*%t(B)

  out = sum(log(get_logit(Y*M))) - 
    sum(get_lambdastar(x = A, rhos = tilde_rhos, lambdas = tilde_lambdas) * A) - 
    sum(get_lambdastar(x = B, rhos = rhos, lambdas = lambdas) * B) + 
    sum((tilde_alpha - 1) * log(tilde_rhos) + (tilde_beta - 1) * log(1 - tilde_rhos)) + 
    sum((alpha - 1) * log(rhos) + (beta - 1) * log(1 - rhos))

  out = -2*sum(log(get_logit(Y*M)))+(log(I*J)*(1+sum(A!=0)+sum(B!=0)))
  
  return(out)

}


################################# Main function
SSLGBC <- function(Y, A = NULL, B = NULL, mu = NULL, 
                   tilde_lambda_0 = 5, tilde_lambda_1 = 1,
                   lambda_0 = 5, lambda_1 = 1, method = "svd", centering = FALSE, 
                   tilde_alpha = 1, tilde_beta = 1, alpha = 1, beta = 1, 
                   K_init = 50, thisSeed = 123, eta = 0.001, IBP = 1, 
                   tol = 1e-10, max_iter = 500, show_plot = FALSE) {
  
  set.seed(thisSeed)  # Ensure reproducibility
  
  I <- nrow(Y)
  J <- ncol(Y)
  K <- K_init

  Y <- as.matrix(Y)
  if(centering){
    Yc <- Y - sum(Y) / (I * J)
  } else {
    Yc <- Y
  }

  tilde_lambda_0 <- tilde_lambda_0  # Matrix A's spike parameter
  tilde_lambda_1 <- tilde_lambda_1  # Matrix A's slab parameter
  lambda_0 <- lambda_0              # Matrix B's spike parameter
  lambda_1 <- lambda_1              # Matrix B's slab parameter
  tilde_lambdas <- c(tilde_lambda_0, tilde_lambda_1)
  lambdas <- c(lambda_0, lambda_1)

  A_in <- A
  B_in <- B

  # Initialize A and B
  if (method == "svd") {
    Y_svd <- tryCatch(
      svd(Yc), 
      error = function(e) {
        message("SVD failed. Using random initialization.")
        NULL
      }
    )
    
    if (!is.null(Y_svd)) {
      if (length(Y_svd$d) < K_init) {
        A <- Y_svd$u %*% diag(sqrt(Y_svd$d))
        B <- Y_svd$v %*% diag(sqrt(Y_svd$d))
      } else {
        A <- Y_svd$u[, 1:K_init] %*% diag(sqrt(Y_svd$d[1:K_init]))
        B <- Y_svd$v[, 1:K_init] %*% diag(sqrt(Y_svd$d[1:K_init]))
      }
    } else {
      # Fallback to random initialization
      A <- matrix(rexp(I * K_init, rate = 1), nrow = I, ncol = K_init)
      B <- matrix(rexp(J * K_init, rate = 1), nrow = J, ncol = K_init)
    }
  } else {
    # Random initialization if method is not 'svd'
    A <- matrix(rexp(I * K_init, rate = 1), nrow = I, ncol = K_init)
    B <- matrix(rexp(J * K_init, rate = 1), nrow = J, ncol = K_init)
  }
  if (!is.null(A_in)) A <- A_in
  if (!is.null(B_in)) B <- B_in

  # Initialize other variables
  if(is.null(mu)) mu <- rep(0,I)

  K <- K_init <- min(K_init, ncol(A))
  tilde_nus <- sort(rbeta(K_init, tilde_alpha, tilde_beta), decreasing = TRUE)
  tilde_rhos <- rep(0.5, K_init)
  nus <- sort(rbeta(K_init, alpha, beta), decreasing = TRUE)
  rhos <- rep(0.5, K_init)

  # The main loop: Use the Rcpp implementation
  result <- main_iterations(
    A = A,
    B = B,
    mu = mu,
    Y = Yc,
    tilde_lambdas = tilde_lambdas,
    lambdas = lambdas,
    tilde_rhos = tilde_rhos,
    rhos = rhos,
    tilde_alpha = tilde_alpha,
    tilde_beta = tilde_beta,
    alpha = alpha,
    beta = beta,
    tol = tol,
    max_iter = max_iter,
    IBP = IBP,
    I = I,
    J = J
  )
  
  # Extract results
  A <- result$A
  B <- result$B
  mu <- result$mu
  tilde_rhos <- result$tilde_rhos
  rhos <- result$rhos
  logLikelihood_save <- result$logLikelihood_save
  delta_logLikelihood_save <- result$delta_logLikelihood_save

  # Prepare output
  out <- list(
    A = A,
    B = B,
    mu = mu,
    tilde_rhos = tilde_rhos,
    rhos = rhos,
    counts = sum(colSums(A != 0) + colSums(B != 0)),
    BIC = min(logLikelihood_save)
  )

  # Optionally show the plot
  if (show_plot) {
    plot(1:length(logLikelihood_save), logLikelihood_save, type = ifelse(length(logLikelihood_save) == 1, "o", "l"),
         xlab = "Iterations", ylab = "Log Likelihood")
  }

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
  
  out <- list(A = this_SSLMC$A, B = this_SSLMC$B, tilde_rhos = this_SSLMC$tilde_rhos, rhos = this_SSLMC$rhos, BIC = BIC)
  return(out)
}



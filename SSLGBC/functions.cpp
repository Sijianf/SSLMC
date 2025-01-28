// This file contains Pcpp functions for SSLGBC
// Author: Gengar
// Date: 2024-10-04
// Note:
// (2025-01-19): Changed the algorithm named as SSLGBC.  
// (2024-12-30): Changed the algorithm from using the transferred Q to the origianl Y.  
// (2024-12-07): Changed 'get_proximal()' to passing function: use shadow links to the target A and B matrix
// (2024-12-07): Added more Rcpp functions. 


// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace arma;

// [[Rcpp::export]]
arma::mat get_logit_matrix(const arma::mat& X) {
  
  return 1 / (1 + exp(-X));
  
}


// lambdastar function for double input
// [[Rcpp::export]]
double get_lambdastar_double(double x, double rho, const arma::vec &lambdas) {
  
  double pstar0 = (1 - rho) * lambdas[0] / 2 * exp(-lambdas[0] * std::abs(x));
  double pstar1 = rho * lambdas[1] / 2 * exp(-lambdas[1] * std::abs(x));
  double pstar = pstar1 / (pstar0 + pstar1);
  double lambdastar = lambdas[0] * (1 - pstar) + lambdas[1] * pstar;
  return lambdastar;
  
}


// lambdastar function for vector input
// [[Rcpp::export]]
arma::vec get_lambdastar_vector(const arma::vec& x, double rho, const arma::vec& lambdas) {
  
  // Compute pstar0 and pstar1 for the entire vector
  arma::vec abs_x = arma::abs(x);
  arma::vec pstar0 = (1 - rho) * lambdas[0] / 2 * exp(-lambdas[0] * abs_x);
  arma::vec pstar1 = rho * lambdas[1] / 2 * exp(-lambdas[1] * abs_x);
  
  // Compute pstar
  arma::vec pstar = pstar1 / (pstar0 + pstar1);
  
  // Compute lambdastar for the vector
  arma::vec lambdastar = lambdas[0] * (1 - pstar) + lambdas[1] * pstar;
  
  return lambdastar;
  
}


// lambdastar function for matrix input
// [[Rcpp::export]]
arma::mat get_lambdastar_matrix(const arma::mat &X, const arma::vec &rhos, const arma::vec &lambdas) {
  
  arma::mat lambdastar_matrix(X.n_rows, X.n_cols);
  
  for (arma::uword i = 0; i < X.n_rows; ++i) {
    for (arma::uword k = 0; k < X.n_cols; ++k) {
      double value = X(i, k);
      double rho = rhos[k];
      double pstar0 = (1 - rho) * lambdas[0] / 2 * exp(-lambdas[0] * std::abs(value));
      double pstar1 = rho * lambdas[1] / 2 * exp(-lambdas[1] * std::abs(value));
      double pstar = pstar1 / (pstar0 + pstar1);
      lambdastar_matrix(i, k) = lambdas[0] * (1 - pstar) + lambdas[1] * pstar;
    }
  }
  
  return lambdastar_matrix;
  
}


// [[Rcpp::export]]
double g(double x, double rho, const arma::vec &lambdas, double n){
  
  double lambdastar=get_lambdastar_double(x,rho,lambdas); 
  double pstar0 = (1 - rho) * lambdas[0] / 2 * exp(-lambdas[0] * std::abs(x));
  double pstar1 = rho * lambdas[1] / 2 * exp(-lambdas[1] * std::abs(x));
  double pstar = pstar1 / (pstar0 + pstar1);
  return pow((lambdastar-lambdas[1]),2)+2*n/4*log(pstar);
  
}


// [[Rcpp::export]]
double get_delta(double rho, const arma::vec &lambdas, double n){
  
  if (lambdas[0] == lambdas[1]){
    return 4*lambdas[1];
  } else {
    if (g(0, rho, lambdas, n) > 0){
      double pstar0 = (1 - rho) * lambdas[0] / 2 * exp(-lambdas[0] * 0);
      double pstar1 = rho * lambdas[1] / 2 * exp(-lambdas[1] * 0);
      double pstar = pstar1 / (pstar0 + pstar1);
      return sqrt(8 * n * log(1/pstar)) + 4*lambdas[1];
    }  else {
      return 4*get_lambdastar_double(0,rho,lambdas);
    }
  }
  
}


// [[Rcpp::export]]
arma::vec soft_thresholding_vector(const arma::vec& z, const arma::vec& lambdastar, double n, double delta) {
  // Compute absolute values and signs
  arma::vec abs_z = arma::abs(z);
  arma::vec sign_z = arma::sign(z);
  
  // Apply the condition |z| > Delta
  arma::uvec mask_uvec = (abs_z > delta); // Indicator as arma::uvec
  arma::vec mask = arma::conv_to<arma::vec>::from(mask_uvec); // Convert to arma::vec
  
  // Compute (|z| - lambdastar)_+ * sign(z)
  arma::vec thresholded = arma::max(abs_z - lambdastar, arma::zeros<arma::vec>(sign_z.n_elem)) % sign_z;
  
  // Apply the mask and scale by 1/n
  arma::vec result = (thresholded % mask) / n;
  
  return result;
}


// [[Rcpp::export]]
Rcpp::List get_rhos(const arma::mat& mat, double alpha, double beta, double tol = 0) {
  int K = mat.n_cols;
  int N = mat.n_rows;
  arma::vec rhos(K, fill::zeros);
  arma::vec counts(K, fill::zeros);

  for (int k = 0; k < K; ++k) {
    int q = sum(abs(mat.col(k)) > tol);
    counts[k] = q;
    rhos[k] = (alpha + q) / (alpha + beta + N);
  }

  return Rcpp::List::create(Rcpp::Named("rhos") = rhos, Rcpp::Named("counts") = counts);
}


// [[Rcpp::export]]
double get_logLikelihood(const arma::mat& Y, const arma::vec mu, 
                         const arma::mat& A, const arma::mat& B,
                         const arma::vec& tilde_rhos, const arma::vec& rhos,
                         const arma::vec& tilde_lambdas, const arma::vec& lambdas,
                         double tilde_alpha, double tilde_beta, 
                         double alpha, double beta) {
  // int I = Y.n_rows;
  int J = Y.n_cols;
  
  // Compute Theta
  arma::vec ones = arma::ones(J);
  arma::mat Theta = mu * ones.t() + A * B.t();

  // Compute log_Pi
  arma::mat Pi = 1 / (1 + exp(-Theta));
  arma::mat log_Pi = arma::log(Pi);
  
  // Compute lambdastar
  arma::mat lambdastar_A = get_lambdastar_matrix(A, tilde_rhos, tilde_lambdas) * abs(A);
  arma::mat lambdastar_B = get_lambdastar_matrix(B, rhos, lambdas) * abs(B);
  
  // Compute the sums:
  double log_likelihood_data = arma::accu(log_Pi - (1 - Y) % Theta);
  double pen_A = -arma::accu(lambdastar_A);
  double pen_B = -arma::accu(lambdastar_B);
  
  // Combine all components
  double log_likelihood = log_likelihood_data + pen_A + pen_B;
  
  // Additional adjustment term
  // double BIC = -2 * log_likelihood + (std::log(I*J) * (1 + accu(A != 0) + accu(B != 0)));

  return log_likelihood;
}


// [[Rcpp::export]]
arma::mat get_X(const arma::mat& Theta, const arma::mat& Y) {
  
  arma::mat Pi = 1 / (1 + exp(-Theta));
  arma::mat X = Theta + 4 * (Y - Pi);
  
  return X;
}


// [[Rcpp::export]]
void update_mu(arma::vec& mu, const arma::mat& X, const arma::mat& A, const arma::mat& B) {
  
  // Compute the sum of rows
  arma::vec rowsums = arma::sum(X - A * B.t(), 1);
  
  // Update mu in place
  mu = rowsums / X.n_cols;
}


// [[Rcpp::export]]
void update_A_B(arma::mat& A, arma::mat& B, 
                const arma::mat& A_momentum, const arma::mat& B_momentum, 
                const arma::mat& X, const arma::vec& mu, 
                const arma::vec& tilde_lambdas, const arma::vec& lambdas,
                const arma::vec& tilde_rhos, const arma::vec& rhos) {
  
  int J = X.n_cols; // Columns in X
  int K = A.n_cols; // Number of columns in A (and B)
  arma::vec ones(J, arma::fill::ones);
  
  for (int k = 0; k < K; ++k) {
    // Get X_k excluding the k-th component
    arma::mat Xstar_k = X - mu * ones.t() - A_momentum * B_momentum.t() + A_momentum.col(k) * B_momentum.col(k).t();
    
    // Update a_k
    arma::vec b_k = B_momentum.col(k); // Current b_k
    double n_a = arma::dot(b_k, b_k);
    
    if (n_a > 0) {
      arma::vec xstar_k = Xstar_k * b_k; 
      arma::vec lambdastar_a = get_lambdastar_vector(A_momentum.col(k), tilde_rhos(k), tilde_lambdas); // Get lambdastar for a_k
      double delta_a = get_delta(tilde_rhos(k), tilde_lambdas, n_a);
      
      // Apply soft-thresholding to update a_k
      A.col(k) = soft_thresholding_vector(xstar_k, lambdastar_a, n_a, delta_a);
    } else {
      A.col(k).zeros(); // Set a_k to zero if b_k_norm is zero
    }
    
    // Update b_k
    arma::vec a_k = A_momentum.col(k); // Current a_k
    double n_b = arma::dot(a_k, a_k);
    
    if (n_b > 0) {
      arma::vec xdagger_k = Xstar_k.t() * a_k; 
      arma::vec lambdastar_b = get_lambdastar_vector(B_momentum.col(k), rhos(k), lambdas); // Get lambdastar for b_k
      double delta_b = get_delta(rhos(k), lambdas, n_b);
      
      // Apply soft-thresholding to update b_k
      B.col(k) = soft_thresholding_vector(xdagger_k, lambdastar_b, n_b, delta_b);
    } else {
      B.col(k).zeros(); // Set b_k to zero if a_k_norm is zero
    }
  }

}


// [[Rcpp::export]]
Rcpp::List main_iterations(arma::mat A, 
                           arma::mat B, 
                           arma::vec mu, 
                           arma::mat Y,
                           arma::vec tilde_lambdas,
                           arma::vec lambdas,
                           arma::vec tilde_rhos,
                           arma::vec rhos,
                           double tilde_alpha,
                           double tilde_beta,
                           double alpha,
                           double beta,
                           double tol,
                           int max_iter, 
                           int IBP,
                           int I,
                           int J) {
  
  // Prepare work
  arma::vec ones = arma::ones(J);
  arma::mat Theta = mu * ones.t() + A * B.t();
  arma::mat X = get_X(Theta, Y);
  arma::mat A_lag = A;
  arma::mat B_lag = B;
  
  double logLikelihood_old = get_logLikelihood(Y, mu, A, B, tilde_rhos, rhos, tilde_lambdas, lambdas, tilde_alpha, tilde_beta, alpha, beta);
  double delta_logLikelihood_old = 1.0;
  arma::vec delta_logLikelihood_save;
  arma::vec logLikelihood_save;
  //Rcpp::Rcout << "logLikelihood_old: " << logLikelihood_old << std::endl;
  
  int K = rhos.size();
  
  for (int i = 1; i <= max_iter; ++i) {
    
    if (i > 1) {
      arma::mat A_momentum = A + (i - 2.0) / (i + 1.0) * (A - A_lag); 
      arma::mat B_momentum = B + (i - 2.0) / (i + 1.0) * (B - B_lag); 
      
      A_lag = A;  // Store the previous state of A
      B_lag = B;  // Store the previous state of B
      
      // Call the function to update A and B directly (since they are linked in the function itself)
      update_A_B(A, B, A_momentum, B_momentum, X, mu, tilde_lambdas, lambdas, tilde_rhos, rhos);
      Rcpp::Rcout << "AB updated: " << A << std::endl;
      
      // Update mu
      update_mu(mu, X, A, B);
      Rcpp::Rcout << "mu updated: " << mu.t() << std::endl;
    }
    
    // Update sparsity
    Rcpp::List A_sparse = get_rhos(A, tilde_alpha, tilde_beta, tol);
    Rcpp::List B_sparse = get_rhos(B, alpha, beta, tol);
     
    // Extract and sum counts
    arma::vec A_counts = Rcpp::as<arma::vec>(A_sparse["counts"]);
    arma::vec B_counts = Rcpp::as<arma::vec>(B_sparse["counts"]);
    arma::vec counts = A_counts + B_counts;
    
    // Extract rhos
    tilde_rhos = Rcpp::as<arma::vec>(A_sparse["rhos"]);
    rhos = Rcpp::as<arma::vec>(B_sparse["rhos"]);
    
    // Re-order the columns
    arma::uvec re_order;
    if (IBP == 1) {
      re_order = sort_index(counts, "descend");
    } else { 
      re_order = sort_index(counts, "descend");
    } 
    
    A = A.cols(re_order);
    B = B.cols(re_order);
    A_lag = A_lag.cols(re_order);
    B_lag = B_lag.cols(re_order);
    tilde_rhos = tilde_rhos(re_order);
    rhos = rhos(re_order);
     
    // Remove all zero columns
    arma::uvec keep = find((sum(A == 0, 0) < I) % (sum(B == 0, 0) < J));
    K = keep.n_elem;
     
    if (K <= 2) {
      break;
    } 
    
    A = A.cols(keep);
    B = B.cols(keep);
    A_lag = A_lag.cols(keep);
    B_lag = B_lag.cols(keep);
    tilde_rhos = tilde_rhos(keep);
    rhos = rhos(keep);
    
    // Update Theta and X
    Theta = mu * ones.t() + A * B.t();
    X = get_X(Theta, Y);
     
    // Stop earlier by log likelihood
    if (i > 1) {
      double logLikelihood = get_logLikelihood(Y, mu, A, B, tilde_rhos, rhos, tilde_lambdas, lambdas, tilde_alpha, tilde_beta, alpha, beta);
      Rcpp::Rcout << "logLikelihood: " << logLikelihood << std::endl;
      
      if (std::isfinite(abs(logLikelihood - logLikelihood_old) / abs(logLikelihood_old))) {
        double delta_logLikelihood = abs(logLikelihood - logLikelihood_old) / abs(logLikelihood_old);
         
        logLikelihood_old = logLikelihood;
        delta_logLikelihood_old = delta_logLikelihood;
         
        delta_logLikelihood_save = arma::join_vert(delta_logLikelihood_save, arma::vec({delta_logLikelihood_old}));
        logLikelihood_save = arma::join_vert(logLikelihood_save, arma::vec({logLikelihood_old}));
         
        if (delta_logLikelihood < tol) {
          break;
        }
      } 
    }
  }
  
  return Rcpp::List::create(
    Rcpp::Named("A") = A,
    Rcpp::Named("B") = B,
    Rcpp::Named("mu") = mu,
    Rcpp::Named("tilde_rhos") = tilde_rhos,
    Rcpp::Named("rhos") = rhos,
    Rcpp::Named("delta_logLikelihood_save") = delta_logLikelihood_save,
    Rcpp::Named("logLikelihood_save") = logLikelihood_save
  );
} 


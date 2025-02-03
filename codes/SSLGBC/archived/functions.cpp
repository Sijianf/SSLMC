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
arma::mat get_logit_matrix(const arma::mat& x) {
  return 1 / (1 + exp(-x));
}


// Function for double input
// [[Rcpp::export]]
double get_lambdastar_double(double x, double thetas, const arma::vec &lambdas) {
  double pstar0 = (1 - thetas) * lambdas[0] / 2 * exp(-lambdas[0] * std::abs(x));
  double pstar1 = thetas * lambdas[1] / 2 * exp(-lambdas[1] * std::abs(x));
  double pstar = pstar1 / (pstar0 + pstar1);
  double lambdastar = lambdas[0] * (1 - pstar) + lambdas[1] * pstar;
  return lambdastar;
}


// Function for matrix input
// [[Rcpp::export]]
arma::mat get_lambdastar_matrix(const arma::mat &x, const arma::vec &thetas, const arma::vec &lambdas) {
  arma::mat lambdastar_matrix(x.n_rows, x.n_cols);
  
  for (arma::uword i = 0; i < x.n_rows; ++i) {
    for (arma::uword k = 0; k < x.n_cols; ++k) {
      double value = x(i, k);
      double theta = thetas[k];
      double pstar0 = (1 - theta) * lambdas[0] / 2 * exp(-lambdas[0] * std::abs(value));
      double pstar1 = theta * lambdas[1] / 2 * exp(-lambdas[1] * std::abs(value));
      double pstar = pstar1 / (pstar0 + pstar1);
      lambdastar_matrix(i, k) = lambdas[0] * (1 - pstar) + lambdas[1] * pstar;
    }
  }
  
  return lambdastar_matrix;
}


// [[Rcpp::export]]
double g(double x, double theta, double eta, const arma::vec &lambdas){
  double value=get_lambdastar_double(x,theta,lambdas); 
  double pstar0 = (1 - theta) * lambdas[0] / 2 * exp(-lambdas[0] * std::abs(x));
  double pstar1 = theta * lambdas[1] / 2 * exp(-lambdas[1] * std::abs(x));
  double pstar = pstar1 / (pstar0 + pstar1);
  return pow((value-lambdas[1]),2)+2/eta*log(pstar);
}


// [[Rcpp::export]]
double get_delta(double theta, double eta, const arma::vec &lambdas){
  if (lambdas[0] == lambdas[1]){
    return eta * lambdas[1];
  } else {
    if (g(0, theta, eta, lambdas) > 0){
      double pstar0 = (1 - theta) * lambdas[0] / 2 * exp(-lambdas[0] * 0);
      double pstar1 = theta * lambdas[1] / 2 * exp(-lambdas[1] * 0);
      double pstar = pstar1 / (pstar0 + pstar1);
      return sqrt(2 * eta * log(1/pstar)) + eta * lambdas[1];
    }  else {
      return eta * get_lambdastar_double(0,theta,lambdas);
    }
  }
}


// [[Rcpp::export]]
double soft_thresholding(double x, double lambdastar, double eta, double delta) {
  double s = 0;
  if (x > 0) s = 1;
  else if (x < 0) s = -1;
  if (fabs(x) <= delta) {
    return(0);
  } else { 
    double temp;
    temp = fabs(x) - eta * lambdastar;
    if (temp > 0) {
      return(temp * s);
    } else {
      return(0);  
    }
  }
}


// [[Rcpp::export]]
void get_proximal(const arma::mat& Y, double eta, double xi, const arma::vec& mu,
                  const arma::mat& A_momentum, const arma::mat& B_momentum, 
                  arma::mat& A, arma::mat& B,
                  const arma::vec& tilde_lambdas, const arma::vec& lambdas,
                  const arma::vec& tilde_thetas, const arma::vec& thetas) {
    int I = Y.n_rows;
    int J = Y.n_cols;
    int K = A_momentum.n_cols;
    arma::vec one_vec = arma::ones(J);
    
    // Compute shared terms using A_momentum and B_momentum
    arma::mat W = (1 + xi * Y - Y) / (1 + arma::exp(- mu * one_vec.t() - A_momentum * B_momentum.t()));
    arma::mat dA = - xi * Y * B_momentum + W * B_momentum;
    arma::mat dB = - xi * Y.t() * A_momentum + W.t() * A_momentum;

    // Update A using A_momentum
    for (int k = 0; k < K; ++k) {
        for (int i = 0; i < I; ++i) {
            double lambdastar = get_lambdastar_double(A_momentum(i, k), tilde_thetas(k), tilde_lambdas);
            double delta = get_delta(tilde_thetas(k), eta, tilde_lambdas);
            A(i, k) = soft_thresholding(A_momentum(i, k) - eta * dA(i, k), lambdastar, eta, delta);
        }
    }

    // Update B using B_momentum
    for (int k = 0; k < K; ++k) {
        for (int j = 0; j < J; ++j) {
            double lambdastar = get_lambdastar_double(B_momentum(j, k), thetas(k), lambdas);
            double delta = get_delta(thetas(k), eta, lambdas);
            B(j, k) = soft_thresholding(B_momentum(j, k) - eta * dB(j, k), lambdastar, eta, delta);
        }
    }
}


// [[Rcpp::export]]
Rcpp::List rescale_A_B(arma::mat A, arma::mat B) {
  int K = A.n_cols;
  int N = A.n_rows;
  int G = B.n_rows;
  double A_norm = 0; 
  double B_norm = 0;
  rowvec d = ones<rowvec>(K);
  
  for(int k = 0; k < K; k++) {
    for (int i = 0; i < N; i++) {
      A_norm += fabs(A(i, k));
    }
    for (int j = 0; j < G; j++) {
      B_norm += fabs(B(j, k));
    }
    if ((A_norm > 0) && (B_norm > 0)) {
      d(k) = pow(A_norm / B_norm, 0.5);
    }
    A_norm = 0;
    B_norm = 0;
  }
  A.each_row() /= d;
  B.each_row() %= d;
  
  return Rcpp::List::create(Rcpp::Named("A") = A,
                            Rcpp::Named("B") = B);
  
}


// [[Rcpp::export]]
Rcpp::List get_thetas(const arma::mat& mat, double a, double b, double tol = 0) {
  int K = mat.n_cols;
  int N = mat.n_rows;
  arma::vec thetas(K, fill::zeros);
  arma::vec counts(K, fill::zeros);

  for (int k = 0; k < K; ++k) {
    int q = sum(abs(mat.col(k)) > tol);
    counts[k] = q;
    thetas[k] = (a + q) / (a + b + N);
  }

  return Rcpp::List::create(Rcpp::Named("thetas") = thetas, Rcpp::Named("counts") = counts);
}


// [[Rcpp::export]]
double get_logLikelihood(const arma::mat& Y, double xi, const arma::vec mu, 
                         const arma::mat& A, const arma::mat& B,
                         const arma::vec& tilde_thetas, const arma::vec& thetas,
                         const arma::vec& tilde_lambdas, const arma::vec& lambdas,
                         double tilde_alpha, double tilde_beta, double alpha, double beta) {
  int I = Y.n_rows;
  int J = Y.n_cols;
  arma::vec one_vec = arma::ones(J);
  
  // Compute the latent matrix M
  arma::mat M = mu * one_vec.t() + A * B.t();

  // Likelihood component
  arma::mat log_likelihood_matrix = xi * Y % M - (xi * Y + 1 - Y) % arma::log(1 + arma::exp(M));
  double log_likelihood_data = arma::accu(log_likelihood_matrix);

  // Initialize accumulators
  double lambdastar_sum_A = 0.0;
  double lambdastar_sum_B = 0.0;

  // Compute lambdastar sums for A
  for (int i = 0; i < I; ++i) {
    for (int k = 0; k < A.n_cols; ++k) {
      lambdastar_sum_A += get_lambdastar_double(A(i, k), tilde_thetas(k), tilde_lambdas) * abs(A(i, k));
    }
  }

  // Compute lambdastar sums for B
  for (int j = 0; j < J; ++j) {
    for (int k = 0; k < B.n_cols; ++k) {
      lambdastar_sum_B += get_lambdastar_double(B(j, k), thetas(k), lambdas) * abs(B(j, k));
    }
  }
  
  // Combine all components
  double total_log_likelihood = log_likelihood_data - lambdastar_sum_A - lambdastar_sum_B;
  
  // Additional adjustment term
  double BIC = -2 * total_log_likelihood + (2 * (1 + accu(A != 0) + accu(B != 0)));

  return BIC;
}


// [[Rcpp::export]]
Rcpp::List main_iterations(int max_iter, 
                           arma::mat A, 
                           arma::mat B, 
                           arma::vec mu, 
                           arma::mat Y,
                           double eta, 
                           double xi, 
                           arma::vec tilde_lambdas,
                           arma::vec lambdas,
                           arma::vec tilde_thetas,
                           arma::vec thetas,
                           double tilde_alpha,
                           double tilde_beta,
                           double alpha,
                           double beta,
                           double tol,
                           int IBP,
                           int I,
                           int J) {
  
  arma::mat A_lag = A;
  arma::mat B_lag = B;
  double logLikelihood_old = get_logLikelihood(Y, xi, mu, A, B, tilde_thetas, thetas, tilde_lambdas, lambdas, tilde_alpha, tilde_beta, alpha, beta);
  double delta_logLikelihood_old = 1000;
  arma::vec delta_logLikelihood_save;
  arma::vec logLikelihood_save;

  int K = thetas.size();

  for (int i = 1; i <= max_iter; ++i) {

  if (i > 1) {
      arma::mat A_momentum = A + (i - 2.0) / (i + 1.0) * (A - A_lag);
      arma::mat B_momentum = B + (i - 2.0) / (i + 1.0) * (B - B_lag);

      A_lag = A;  // Store the previous state of A
      B_lag = B;  // Store the previous state of B

      // Call the function to update A and B directly (since they are linked in the function itself)
      get_proximal(Y, eta, xi, mu, A_momentum, B_momentum, A, B, tilde_lambdas, lambdas, tilde_thetas, thetas);
  }
  
    // Update mu
    arma::vec one_vec = arma::ones(J);
    arma::mat p = 1 / (1 + arma::exp(-(mu * one_vec.t() + A * B.t())));  // Compute p_ij
    for (int i = 0; i < I; ++i) {
      double sum_xi_yij = 0.0;
      double sum_adjusted_pij = 0.0;
      double sum_xi_yij_plus_other = 0.0;
      
      for (int j = 0; j < J; ++j) {
        sum_xi_yij += xi * Y(i, j);
        sum_adjusted_pij += (xi * Y(i, j) + 1 - Y(i, j)) * p(i, j);
        sum_xi_yij_plus_other += xi * Y(i, j) + 1 - Y(i, j);
      }
      
      // Update mu_i
      mu[i] += 4 * (1 / sum_xi_yij_plus_other) * (sum_xi_yij - sum_adjusted_pij);
    }

    // Update sparsity
    Rcpp::List A_sparse = get_thetas(A, tilde_alpha, tilde_beta, tol);
    Rcpp::List B_sparse = get_thetas(B, alpha, beta, tol);

    // Extract and sum counts
    arma::vec A_counts = Rcpp::as<arma::vec>(A_sparse["counts"]);
    arma::vec B_counts = Rcpp::as<arma::vec>(B_sparse["counts"]);
    arma::vec counts = A_counts + B_counts;

    // Extract thetas
    tilde_thetas = Rcpp::as<arma::vec>(A_sparse["thetas"]);
    thetas = Rcpp::as<arma::vec>(B_sparse["thetas"]);

    // Get the change of probability with one cluster leave out
    // arma::mat overall_prob = 1 / (1 + exp(-A * B.t()));
    // arma::vec bicluster_prob(K);

    // for (int k = 0; k < K; ++k) {
    //   arma::mat A_k = A; A_k.col(k).zeros();
    //   arma::mat B_k = B; B_k.col(k).zeros();
    //   arma::mat prob_diff = 1 / (1 + exp(-A_k * B_k.t())) - overall_prob; // Element-wise difference
    //   bicluster_prob[k] = arma::abs(prob_diff).max();                      // Maximum absolute value
    // }

    // Re-order the columns
    arma::uvec re_order;
    if (IBP == 1) {
      re_order = sort_index(counts, "descend");
    } else {
      re_order = sort_index(counts, "descend");
      // re_order = sort_index(bicluster_prob, "descend");
    }

    A = A.cols(re_order);
    B = B.cols(re_order);
    A_lag = A_lag.cols(re_order);
    B_lag = B_lag.cols(re_order);
    tilde_thetas = tilde_thetas(re_order);
    thetas = thetas(re_order);

    // Remove all zero columns
    arma::uvec keep = find(sum(A == 0, 0) < I && sum(B == 0, 0) < J);
    K = keep.n_elem;

    if (K <= 2) {
      break;
    }

    A = A.cols(keep);
    B = B.cols(keep);
    A_lag = A_lag.cols(keep);
    B_lag = B_lag.cols(keep);
    tilde_thetas = tilde_thetas(keep);
    thetas = thetas(keep);

    // Stop earlier by log likelihood
    if (i > 1) {
      double logLikelihood = get_logLikelihood(Y, xi, mu, A, B, tilde_thetas, thetas, tilde_lambdas, lambdas, tilde_alpha, tilde_beta, alpha, beta);
      
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
    Rcpp::Named("tilde_thetas") = tilde_thetas,
    Rcpp::Named("thetas") = thetas,
    Rcpp::Named("delta_logLikelihood_save") = delta_logLikelihood_save,
    Rcpp::Named("logLikelihood_save") = logLikelihood_save
  );
}



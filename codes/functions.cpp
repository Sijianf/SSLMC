// This file contains Pcpp functions for SSLMC
// Author: Gengar
// Date: 2024-10-04

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace arma;


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
Rcpp::List get_proximal(const arma::mat Y, double xi, double eta, 
                        arma::mat A, arma::mat B, 
                        const arma::vec tilde_lambdas, const arma::vec lambdas,
                        const arma::vec tilde_thetas, const arma::vec thetas) {

  int I = Y.n_rows;
  int J = Y.n_cols;
  int K = A.n_cols;
  arma::mat proximal_A = A;
  arma::mat proximal_B = B;

  arma::mat W = (1 + xi * Y - Y) /
                    (1 + exp(- A * B.t()));

  arma::mat dA = - xi * Y * B + W * B;
  arma::mat dB = - xi * Y.t() * A + W.t() * A;

  for (int k = 0; k < K; ++k) {
    for (int i = 0; i < I; ++i) {
      double lambdastar = get_lambdastar_double(A(i, k), tilde_thetas(k), tilde_lambdas);
      double delta = get_delta(tilde_thetas(k), eta, tilde_lambdas);
      proximal_A(i, k) = soft_thresholding(A(i, k) - eta * dA(i, k), lambdastar, eta, delta);
    }

    for (int j = 0; j < J; ++j) {
      double lambdastar = get_lambdastar_double(B(j, k), thetas(k), lambdas);
      double delta = get_delta(thetas(k), eta, lambdas);
      proximal_B(j, k) = soft_thresholding(B(j, k) - eta * dB(j, k), lambdastar, eta, delta);
    }
  }

  return Rcpp::List::create(Rcpp::Named("proximal_A") = proximal_A,
                            Rcpp::Named("proximal_B") = proximal_B);
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








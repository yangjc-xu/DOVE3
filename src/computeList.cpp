#include <iostream>
#include <armadillo>
#include <RcppArmadillo.h>
using namespace Rcpp;
// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
arma::field<arma::vec> computeList(arma::field<arma::vec> object, arma::field<arma::vec> add, arma::vec index) {
  int len = index.n_elem;
  arma::field<arma::vec> object_new = object;
  for(int i = 0; i < len; ++i){
    object_new(index(i)) = arma::join_vert(object(index(i)), add(i));
  }
  return(object_new);
}

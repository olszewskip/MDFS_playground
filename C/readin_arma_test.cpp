// Copyright: me
// g++ arma_1.cpp -DARMA_DONT_USE_WRAPPER -lopenblas -llapack

#include <armadillo>
#include <iostream>
#include <vector>

int main() {
  int n_cols = 2;
  // int cols[] = {200, 500};
  // std::vector<int> cols_vector(cols, cols + n_cols);
  arma::uvec cols_vec = {200, 500};

  arma::Mat<double> M;

  M.load("madelon.csv");

  arma::Mat<double> column_bunch = M.cols(cols_vec);
  
}

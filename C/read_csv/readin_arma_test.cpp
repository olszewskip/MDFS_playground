// Copyright: me
// g++ arma_1.cpp -DARMA_DONT_USE_WRAPPER -lopenblas -llapack

#include <armadillo>
#include <iostream>
#include <vector>

int main() {
  arma::uvec cols_vec = {9, 99, 199};

  arma::Mat<double> M;
  arma::Mat<double> column_bunch;

  int n_loops = 100;

  for (int i = 0; i < n_loops; i++) {
    M.load("madelonX16.csv");
    column_bunch = M.cols(cols_vec);
  }
}

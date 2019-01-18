// Copyright: me
// g++ arma_1.cpp -DARMA_DONT_USE_WRAPPER -lopenblas -llapack

#include <armadillo>
#include <iostream>
#include <vector>

int main(int argc, char *argv[]) {
  arma::uvec cols_vec = {9, 99, 199};

  arma::Mat<double> M;
  arma::Mat<double> column_bunch;

  M.load(argv[1]);
  int runsB = 100000;

  for (int i = 0; i < runsB; i++) {
    column_bunch = M.cols(cols_vec);
  }
}

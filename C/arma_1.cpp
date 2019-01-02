// Copyright: me
// g++ arma_1.cpp -DARMA_DONT_USE_WRAPPER -lopenblas -llapack

#include <armadillo>
#include <iostream>

int main(int argc, const char **argv) {
  arma::Mat<double> M;

  M.load("madelon.csv");

  std::cout << "first element: " << M(0, 0) << std::endl;
  std::cout << "number of rows: " << M.n_rows << std::endl;
  std::cout << "number of cols: " << M.n_cols << std::endl;
  std::cout << "number of elems: " << M.n_elem << std::endl;
  std::cout << "dimensions together: " << arma::size(M) << std::endl;

  std::cout
      << "number of elements in the 1st columns that are greater than 480: "
      << arma::sum(M.col(0) > 480) << std::endl;

  arma::Mat<double> M2 = M.cols(1, 2);
  M2(0, 0) = 123;
  std::cout << "first element again: " << M(0, 0) << std::endl;

  return 0;
}

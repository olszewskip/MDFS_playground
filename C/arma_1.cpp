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

  arma::Mat<double> M2 = M.cols(0, 2);
  std::cout << "number of columns in the new copy: " << M2.n_cols << std::endl;
  M2(0, 0) = 123;
  std::cout << "first element of the copy: " << M2(0, 0) << std::endl;
  std::cout << "first element of the original: " << M(0, 0) << std::endl;

  // overflowing copy: runtime error
  arma::Mat<double> M3 = M.cols(M.n_cols, M.n_cols + 1);

  return 0;
}

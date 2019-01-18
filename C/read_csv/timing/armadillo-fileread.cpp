// Copyright: me
// g++ arma_1.cpp -DARMA_DONT_USE_WRAPPER -lopenblas -llapack

#include <armadillo>
#include <iostream>
#include <vector>

int main(int argc, char *argv[]) {
  arma::Mat<double> M;

  int runsA = atoi(argv[2]);

  for (int i = 0; i < runsA; i++) {
    M.load(argv[1]);
  }

  std::cout << runsA << " runs in " << argv[1] << "\n";
}

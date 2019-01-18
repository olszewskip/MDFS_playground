// Copyright: me
// g++ arma_1.cpp -DARMA_DONT_USE_WRAPPER -lopenblas -llapack

#include <armadillo>
#include <chrono>
#include <iostream>
#include <vector>

using Clock = std::chrono::steady_clock;
using std::chrono::milliseconds;
using std::chrono::duration_cast;

int main(int argc, char *argv[]) {
  arma::uvec cols_vec = {9, 99, 199};

  arma::Mat<double> M;
  arma::Mat<double> column_bunch;

  M.load(argv[1]);
  int runsB = atoi(argv[2]);

  auto time_0 = Clock::now();
  for (int i = 0; i < runsB; i++) {
    column_bunch = M.cols(cols_vec);
  }
  auto time_1 = Clock::now();
  milliseconds diff = duration_cast<milliseconds>(time_1 - time_0);

  std::cout << runsB << " runs in " << argv[1] << ": " << diff.count() << " millisec.\n";
}

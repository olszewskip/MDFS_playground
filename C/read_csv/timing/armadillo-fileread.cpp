// Copyright: me
// g++ arma_1.cpp -DARMA_DONT_USE_WRAPPER -lopenblas -llapack

#include <armadillo>
#include <chrono>
#include <iostream>
#include <vector>

using Clock = std::chrono::steady_clock;
using std::chrono::duration_cast;
using std::chrono::milliseconds;

int main(int argc, char *argv[]) {
  arma::Mat<double> M;

  int runsA = atoi(argv[2]);

  auto time_0 = Clock::now();
  for (int i = 0; i < runsA; i++) {
    M.load(argv[1]);
  }
  auto time_1 = Clock::now();
  milliseconds diff = duration_cast<milliseconds>(time_1 - time_0);

  std::cout << runsA << " runs in " << argv[1] << ": " << diff.count()
            << " millisec.\n";
}

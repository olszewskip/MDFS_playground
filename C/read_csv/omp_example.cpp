// Copyright: me
#include <omp.h>
#include <iostream>

int main() {
  int n = 1E8;
  double step = 1.0 / n;
  double s = 0;

#pragma omp parallel
  {
    double x;
#pragma omp for reduction(+ : s)
    for (int i = 0; i < n; i++) {
      x = (i + 0.5) * step;
      s += 4.0 / (1 + x * x);
    }
  }
  double result = step * s;
  std::cout << result << std::endl;
}

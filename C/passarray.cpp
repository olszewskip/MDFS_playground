// Copyright: me
#include <algorithm>
#include <fstream>
#include <iostream>
#include <iterator>
#include <sstream>
#include <string>
#include <vector>

double avg(double* arr, int size);

int main() {
  double A[3] = {1., 2., 3.};
  std::cout << avg(A, 3) << std::endl;
}

double avg(double* arr, int size) {
  double sum = 0;
  for (int i = 0; i < size; i++) {
    sum += arr[i];
  }
  return sum / size;
}

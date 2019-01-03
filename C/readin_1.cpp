// Copyright: internet
#include <algorithm>
#include <fstream>
#include <iostream>
#include <iterator>
#include <sstream>
#include <string>
#include <vector>

int main() {
  // using namespace std;

  // open file, create vector of vectors M
  std::ifstream in("simple.csv");
  std::vector<std::vector<double>> M;
  // stream from file to M
  if (in) {
    std::string line;

    while (getline(in, line)) {
      std::stringstream nums(line);
      std::string num;
      M.emplace_back(std::vector<double>());
      while (getline(nums, num, ',')) {
        M.back().emplace_back(stod(num));
      }
    }
  }
  in.close();

  // print M to screen
  for (auto row : M) {
    for (auto num : row) {
      std::cout << num << " ";
    }
    std::cout << '\n';
  }

  std::cout << M.size() << std::endl;  // 5 ??
  std::cout << M[0].size() << std::endl;

  // copy(!) to another vector-vector
  // std::vector<std::vector<double>> M2(M);
  // which is equivalent to
  // std::vector<std::vector<double>> M2 = M;
  // I guess; see also:
  std::vector<std::vector<double>> M2;
  // followed by
  // std::copy(M.begin(), M.end(), std::back_inserter(M2));
  // or
  M2.assign(M.begin(), M.end());

  M[0][0] = 123.;
  std::cout << "M[0][0]=" << M[0][0] << " ,  M2[0][0]=" << M2[0][0]
            << std::endl;

  
}

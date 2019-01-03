// Copyright: internet
#include <algorithm>
#include <fstream>
#include <iostream>
#include <iterator>
#include <sstream>
#include <string>
#include <vector>

void printM(std::vector<std::vector<double>> *M) {
  for (auto row : *M) {
    for (auto num : row) {
      std::cout << num << " ";
    }
    std::cout << '\n';
  }
}

int main() {
  // using namespace std;

  // open file, create vector of vectors M
  std::ifstream in("simple.csv");
  std::vector<std::vector<double>> M;
  // stream from file to M
  if (in) {
    std::string line;

    while (std::getline(in, line)) {
      std::stringstream nums(line);
      std::string num;
      M.emplace_back(std::vector<double>());
      while (std::getline(nums, num, ',')) {
        // std::cout << stod(num) << '\n';
        M.back().emplace_back(stod(num));
      }
    }
  }
  in.close();  // should I do this?

  ///
  std::cout << M.size() << " " << M[0].size() << std::endl;
  printM(&M);

  // open file, create vector of vectors Mt
  std::ifstream in2("simple.csv");
  std::vector<std::vector<double>> Mt;
  // stream from file to Mt in a transposed order

  // prepare global variables (I guess... ?)
  std::string line;
  std::string num;

  // fill the first column of Mt with first line from the csv
  std::getline(in2, line);
  std::stringstream firstnums(line);
  while (std::getline(firstnums, num, ',')) {
    std::vector<double> row = {stod(num)};
    Mt.emplace_back(row);
  }

  // printM(&Mt);

  // fill other columns with proceeding lines from the csv
  while (std::getline(in2, line)) {
    std::stringstream nums(line);
    int i = 0;
    for (auto &row : Mt) {
      /* get num */ std::getline(nums, num, ',');
      std::cout << i << " " << stod(num) << std::endl;
      i++;
      // row.emplace_back(0.0);
    }
  }

  printM(&Mt);

  // if (in) {
  //   std::string line;

  //   while (getline(in, line)) {
  //     std::stringstream nums(line);
  //     std::string num;
  //     M.emplace_back(std::vector<double>());
  //     while (getline(nums, num, ',')) {
  //       M.back().emplace_back(stod(num));
  //     }
  //   }
  // }
  //

  // print M to screen
  printM(&M);

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

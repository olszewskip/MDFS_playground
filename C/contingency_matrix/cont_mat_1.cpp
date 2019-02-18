// copyright: me
#include <math.h>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

void print_matrix(const std::vector<std::vector<int>> &M);

void populate_from_file_transposing(std::vector<std::vector<int>> *M,
                                    std::string file_name,
                                    const char separator = ',');

std::tuple<double, double> work(int dim1, int divisions,
                                const std::vector<int> &X0,
                                const std::vector<int> &X1, int n_classes,
                                const std::vector<int> &y);

int main(int argc, char *argv[]) {
  std::string file_name = argv[1];
  std::vector<std::vector<int>> data;
  populate_from_file_transposing(&data, file_name);

  std::vector<int> X_0, X_1, y;
  X_0 = data[0];
  X_1 = data[1];
  y = data[2];

  print_matrix(data);
  // std::cout << y.data() << '\n';
  // std::cout << y.size() << '\n';
  // std::cout << "-----" << '\n';

  auto [IG_0, IG_1] = work(10, 1, X_0, X_1, 2, y);
  std::cout << "IG_0=" << IG_0 << ", " << "IG_1=" << IG_1 << std::endl;
}

// print matrix to screen
void print_matrix(const std::vector<std::vector<int>> &M) {
  for (auto row : M) {
    for (auto num : row) {
      std::cout << num << " ";
    }
    std::cout << '\n';
  }
  std::cout << "-----\n";
}

// fill the matrix M with values from a csv-file in a way
// that transposition of M equals the data from file
void populate_from_file_transposing(std::vector<std::vector<int>> *M,
                                    std::string file_name,
                                    const char separator) {
  std::ifstream csv_file(file_name);
  if (csv_file) {
    std::string csv_line;
    // get first line from the csv
    // and turn it into the first column in M
    std::getline(csv_file, csv_line);
    std::stringstream text_nums(csv_line);
    std::string text_num;
    while (std::getline(text_nums, text_num, separator)) {
      int num = std::stoi(text_num);
      // for each num in the csv_line create new row in M
      std::vector<int> row = {num};
      (*M).emplace_back(row);
    }
    // get proceeding lines from the csv
    while (std::getline(csv_file, csv_line)) {
      std::stringstream text_nums(csv_line);
      for (auto &row : *M) {
        std::string text_num;
        std::getline(text_nums, text_num, separator);
        int num = std::stoi(text_num);
        row.emplace_back(num);
      }
    }
  }
  csv_file.close();
}

std::tuple<double, double> work(const int kData_dim, const int kDivisions,
                                const std::vector<int> &X0,
                                const std::vector<int> &X1,
                                const int kN_classes,
                                const std::vector<int> &y) {
  const int kC_Xdim = kDivisions + 1;
  const int kC_ydim = kN_classes;

  // zero-initialized
  int contingency_m[kC_Xdim + 1][kC_Xdim + 1][kC_ydim + 1] = {};

  for (int data_idx = 0; data_idx < kData_dim; data_idx++) {
    // std::cout << "data_idx=" << data_idx << std::endl;
    // std::cout << "X0[data_idx]=" << X0[data_idx] << std::endl;
    // std::cout << "X1[data_idx]=" << X1[data_idx] << std::endl;
    // std::cout << "y[data_idx]=" << y[data_idx] << std::endl;
    contingency_m[X0[data_idx]][X1[data_idx]][y[data_idx]] += 1;
    contingency_m[X0[data_idx]][X1[data_idx]][kC_ydim] += 1;
    contingency_m[X0[data_idx]][kC_Xdim][y[data_idx]] += 1;
    contingency_m[kC_Xdim][X1[data_idx]][y[data_idx]] += 1;
    contingency_m[X0[data_idx]][kC_Xdim][kC_ydim] += 1;
    contingency_m[kC_Xdim][X1[data_idx]][kC_ydim] += 1;
  }

  for (int C_yidx = 0; C_yidx < kC_ydim + 1; C_yidx++) {
    for (int C_Xidx_0 = 0; C_Xidx_0 < kC_Xdim + 1; C_Xidx_0++) {
      for (int C_Xidx_1 = 0; C_Xidx_1 < kC_Xdim + 1; C_Xidx_1++) {
        std::cout << contingency_m[C_Xidx_0][C_Xidx_1][C_yidx] << " ";
      }
      std::cout << std::endl;
    }
    std::cout << "--" << std::endl;
  }

  double neg_H = 0, neg_H_y = 0, neg_H_X0 = 0, neg_H_X1 = 0;

  for (int C_Xidx_0 = 0; C_Xidx_0 < kC_Xdim; C_Xidx_0++) {
    for (int C_Xidx_1 = 0; C_Xidx_1 < kC_Xdim; C_Xidx_1++) {
      for (int C_yidx = 0; C_yidx < kC_ydim; C_yidx++) {
        int &count = contingency_m[C_Xidx_0][C_Xidx_1][C_yidx];
        neg_H += count * log2(count);
      }
      int &count_y = contingency_m[C_Xidx_0][C_Xidx_1][kC_ydim];
      neg_H -= count_y * log2(count_y);
    }
  }
  // std::cout << neg_H << std::endl;

  for (int C_Xidx = 0; C_Xidx < kC_Xdim; C_Xidx++) {
    for (int C_yidx = 0; C_yidx < kC_ydim; C_yidx++) {
      int &count_X0 = contingency_m[kC_Xdim][C_Xidx][C_yidx];
      neg_H_X0 += count_X0 * log2(count_X0);
      int &count_X1 = contingency_m[C_Xidx][kC_Xdim][C_yidx];
      neg_H_X1 += count_X1 * log2(count_X1);
    }
    int &count_X0_y = contingency_m[kC_Xdim][C_Xidx][kC_ydim];
    neg_H_X0 -= count_X0_y * log2(count_X0_y);
    int &count_X1_y = contingency_m[C_Xidx][kC_Xdim][kC_ydim];
    neg_H_X1 -= count_X1_y * log2(count_X1_y);
  }

  // std::cout << neg_H_X0 << std::endl;
  return {neg_H - neg_H_X0, neg_H - neg_H_X1};
}

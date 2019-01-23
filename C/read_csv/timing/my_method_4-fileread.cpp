// copyright: me
#include <chrono>
#include <fstream>
#include <iostream>
// #include <iterator>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

void print_matrix(const int dim0, const int dim1, double M[]);

std::tuple<int, int> data_dims(const std::string &file_name,
                               const char separator = ',');

void populate_from_file_transposing(double M[], const int dim0, const int dim1,
                                    const std::string &file_name,
                                    const char separator = ',');

double *return_rows(double source[], int n_rows, int rows[], int dim1);

int main() {
  const std::string file_name = "simple.csv";
  const auto [dim0, dim1] = data_dims(file_name);

  double *M = new double[dim0 * dim1];
  populate_from_file_transposing(M, dim0, dim1, "simple.csv");
  print_matrix(dim1, dim0, M);

  int n_rows = 3;
  int rows[] = {0, 1, 3};
  double *column_bunch = return_rows(M, n_rows, rows, dim0);
  print_matrix(n_rows, dim0, column_bunch);
}

// print matrix to screen
void print_matrix(const int dim0, const int dim1, double M[]) {
  for (int i = 0; i < dim0; i++) {
    for (int j = 0; j < dim1; j++) {
      std::cout << M[i * dim1 + j] << " ";
    }
    std::cout << '\n';
  }
  std::cout << "-----\n";
}

std::tuple<int, int> data_dims(const std::string &file_name,
                               const char separator) {
  int dim0 = 0;
  int dim1 = 0;

  std::ifstream csv_file(file_name);
  std::string first_line;
  std::getline(csv_file, first_line);
  dim0++;

  std::stringstream text_nums(first_line);
  std::string text_placeholder;
  while (std::getline(text_nums, text_placeholder, separator)) {
    dim1++;
  }
  while (std::getline(csv_file, text_placeholder)) {
    dim0++;
  }
  csv_file.close();

  return {dim0, dim1};
}

void populate_from_file_transposing(double M[], const int dim0, const int dim1,
                                    const std::string &file_name,
                                    const char separator) {
  std::ifstream csv_file(file_name);
  std::string csv_line, text_num;

  for (int i = 0; i < dim0; i++) {
    std::getline(csv_file, csv_line);
    std::stringstream text_nums(csv_line);
    for (int j = 0; j < dim1; j++) {
      std::getline(text_nums, text_num, separator);
      double num = std::stod(text_num);
      M[j * dim0 + i] = num;
    }
  }
}

double *return_rows(double source[], int n_rows, int rows[], int dim1) {
  double *column_bunch = new double[n_rows * dim1];
  for (int i = 0; i < n_rows; i++) {
    double *ptr = &source[rows[i] * dim1];
    for (int j = 0; j < dim1; j++) {
      column_bunch[i * dim1 + j] = ptr[j];
    }
  }
  return column_bunch;
}

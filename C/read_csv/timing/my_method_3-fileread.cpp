// copyright: me
#include <chrono>
#include <fstream>
#include <iostream>
// #include <iterator>
#include <sstream>
#include <string>
#include <vector>

using Clock = std::chrono::steady_clock;
using std::chrono::duration_cast;
using std::chrono::milliseconds;

void print_matrix(const int dim1, const std::vector<double *> &M);

int populate_from_file(std::vector<double *> *M, const std::string &file_name,
                       const char separator = ',');

std::vector<double *> return_transposed_columns(
    const int dim0, const std::vector<double *> &source, int n_cols,
    int cols[]);

void clear_matrix(std::vector<double *> *M);

int main(int argc, char *argv[]) {
  // matrix M to polulate
  std::vector<double *> M;
  int dim1;

  int n_cols = 3;
  int cols[] = {0, 1, 2};

  int runsA = atoi(argv[2]);

  auto time_0 = Clock::now();
  // method 3
  for (int i = 0; i < runsA; i++) {
    dim1 = populate_from_file(&M, argv[1]);
    clear_matrix(&M);
  }

  auto time_1 = Clock::now();
  milliseconds diff = duration_cast<milliseconds>(time_1 - time_0);

  std::cout << runsA << " runs in " << argv[1] << ": " << diff.count()
            << " millisec.\n";
}

// print matrix to screen
void print_matrix(const int dim1, const std::vector<double *> &M) {
  for (auto row : M) {
    for (int j = 0; j < dim1; j++) {
      std::cout << row[j] << " ";
    }
    std::cout << '\n';
  }
  std::cout << "-----\n";
}

// fill the matrix M with values from a csv-file
// return its second dimension, the length of rows
int populate_from_file(std::vector<double *> *M, const std::string &file_name,
                       const char separator) {
  // the file to stream to M from
  std::ifstream csv_file(file_name);
  std::string csv_line;
  std::getline(csv_file, csv_line);
  std::stringstream text_nums(csv_line);

  std::vector<double> first_row_vec;
  std::string text_num;
  while (std::getline(text_nums, text_num, ',')) {
    double num = std::stod(text_num);
    // std::cout << num << " \n";
    first_row_vec.emplace_back(num);
  }

  const int dim1 = first_row_vec.size();

  double *first_row = new double[dim1];

  for (int i = 0; i < dim1; i++) {
    first_row[i] = first_row_vec[i];
  }

  (*M).emplace_back(first_row);

  while (std::getline(csv_file, csv_line)) {
    double *row = new double[dim1];
    std::stringstream text_nums(csv_line);
    std::string text_num;
    for (int i = 0; i < dim1; i++) {
      std::getline(text_nums, text_num, ',');
      double num = std::stod(text_num);
      row[i] = num;
    }
    (*M).emplace_back(row);
  }
  csv_file.close();

  return dim1;
}

// get selected columns from a source-matrix
// return them transposed
std::vector<double *> return_transposed_columns(
    const int dim0, const std::vector<double *> &source, int n_cols,
    int cols[]) {
  std::vector<double *> column_bunch;
  column_bunch.reserve(n_cols);

  for (int i = 0; i < n_cols; i++) {
    column_bunch.emplace_back(new double[dim0]);
    for (int j = 0; j < dim0; j++) {
      column_bunch[i][j] = source[j][cols[i]];
    }
  }
  return column_bunch;
}

// free memory allocoted for the matrix
void clear_matrix(std::vector<double *> *M) {
  for (auto row : *M) {
    delete[] row;
  }
  M->clear();
}

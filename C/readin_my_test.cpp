// Copyright: me
// #include <algorithm>
#include <fstream>
#include <iostream>
// #include <iterator>
#include <sstream>
#include <string>
#include <vector>

void print_matrix(std::vector<std::vector<double>> *M);

// there must be a better way of returning M without copying
void populate_from_file(std::vector<std::vector<double>> *M,
                        std::string file_name, const char separator = ',');

std::vector<std::vector<double>> return_columns(
    std::vector<std::vector<double>> *source, int n_cols, int cols[]);

void populate_from_file_transposing(std::vector<std::vector<double>> *M,
                                    std::string file_name,
                                    const char separator = ',');

std::vector<std::vector<double>> return_rows(
    std::vector<std::vector<double>> *source, int n_rows, int rows[]);

int main() {
  // matrix M1 to populate
  std::vector<std::vector<double>> M;
  int n_cols = 2;
  int cols[] = {200, 500};

  // method 1
  // populate_from_file(&M, "madelon.csv");
  // std::vector<std::vector<double>> column_bunch_A =
  //     return_columns(&M, n_cols, cols);

  // method 2
  populate_from_file_transposing(&M, "madelon.csv");
  std::vector<std::vector<double>> column_bunch_A =
      return_rows(&M, n_cols, cols);
}

// print matrix to screen
void print_matrix(std::vector<std::vector<double>> *M) {
  for (auto row : *M) {
    for (auto num : row) {
      std::cout << num << " ";
    }
    std::cout << '\n';
  }
  std::cout << "-----\n";
}

// fill the matrix M with values from a csv-file
void populate_from_file(std::vector<std::vector<double>> *M,
                        std::string file_name, const char separator) {
  // the file to stream to M from
  std::ifstream csv_file(file_name);
  if (csv_file) {
    // parse csv_file with defaults of getline
    std::string csv_line;
    while (std::getline(csv_file, csv_line)) {
      // for each csv_line create a row in the matrix M
      (*M).emplace_back(std::vector<double>());
      // parse csv_line with getline with the separator
      std::stringstream text_nums(csv_line);
      // above: std::stringstream text_nums = std::stringstream(csv_line);
      std::string text_num;
      while (std::getline(text_nums, text_num, separator)) {
        double num = std::stod(text_num);
        // fill the newly created row in M
        (*M).back().emplace_back(num);
      }
    }
  }
  csv_file.close();
}

// get selected columns from a sourse-matrix
// return them transposed
std::vector<std::vector<double>> return_columns(
    std::vector<std::vector<double>> *source, int n_cols, int cols[]) {
  // number of rows in the source
  int n_rows = (*source).size();
  // initialize the returned matrix
  std::vector<std::vector<double>> column_bunch(n_cols);
  // fill in i-th row of column_bunch with cols[i]-th column of source
  for (int i = 0; i < n_cols; i++) {
    int col_idx = cols[i];
    for (int j = 0; j < n_rows; j++) {
      column_bunch[i].emplace_back((*source)[j][col_idx]);
    }
  }
  return column_bunch;
}

// fill the matrix M with values from a csv-file in a way
// that transposition of M equals the data from file
void populate_from_file_transposing(std::vector<std::vector<double>> *M,
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
      double num = std::stod(text_num);
      // for each num in the csv_line create new row in M
      std::vector<double> row = {num};
      (*M).emplace_back(row);
    }
    // get proceeding line from the csv
    while (std::getline(csv_file, csv_line)) {
      std::stringstream text_nums(csv_line);
      for (auto &row : *M) {
        std::string text_num;
        std::getline(text_nums, text_num, separator);
        double num = std::stod(text_num);
        row.emplace_back(num);
      }
    }
  }
  csv_file.close();
}

std::vector<std::vector<double>> return_rows(
    std::vector<std::vector<double>> *source, int n_rows, int rows[]) {
  // initialize the returned matrix
  std::vector<std::vector<double>> row_bunch(n_rows);
  // fill the i-th row of row_bunch with rows[i]-th row of source
  for (int i = 0; i < n_rows; i++) {
    row_bunch[i] = (*source)[rows[i]];
  }
  return row_bunch;
}

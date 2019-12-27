#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include "matrix_class.h"

void parse_csv_dims(const std::string& csv_file_name,
                    const char separator,
                    uint& csv_row_count,
                    uint& csv_column_count) {

  csv_row_count = 0;
  csv_column_count = 0;

  std::ifstream csv_file(csv_file_name);
  std::string first_line;
  std::getline(csv_file, first_line);
  csv_row_count++;

  std::stringstream text_nums(first_line);
  std::string text_placeholder;
  while (std::getline(text_nums, text_placeholder, separator)) {
    csv_column_count++;
  }
  while (std::getline(csv_file, text_placeholder)) {
    csv_row_count++;
  }
  csv_file.close();
}


void populate_from_csv_transposing(Matrix<double>& matrix,
                                   const std::string csv_file_name,
                                   const char separator,
                                   const uint csv_column_count,
                                   const uint csv_row_count){
  std::ifstream csv_file(csv_file_name);
  std::string csv_line, text_num;

  for (int i = 0; i < csv_row_count; i++) {
    std::getline(csv_file, csv_line);
    std::stringstream text_nums(csv_line);
    for (int j = 0; j < csv_column_count; j++) {
      std::getline(text_nums, text_num, separator);
      double num = std::stod(text_num);
      matrix(j, i) = num;
    }
  }
  csv_file.close();
}

int main(int argc, char *argv[]) {
    std::string csv_file_name = argv[1];
    uint csv_row_count, csv_column_count;
    parse_csv_dims(csv_file_name,
                   ',',
                   csv_row_count,
                   csv_column_count);
    // std::cout << csv_row_count << " " << csv_column_count << std::endl;
    Matrix<double> data(csv_column_count, csv_row_count);
    populate_from_csv_transposing(data,
                                  csv_file_name,
                                  ','   ,
                                  csv_column_count, 
                                  csv_row_count);
    data.print();
//    print_matrix(matrix,
//                 csv_column_count,
//                 csv_row_count);
}


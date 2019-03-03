// Copyright: https://isocpp.org/wiki/faq/operator-overloading#matrix-subscript-op
#include <stdexcept>
#include <iostream>

template<typename T>
class Matrix {
public:
  Matrix(unsigned rows, unsigned cols);
  T& operator() (unsigned row, unsigned col);        // Subscript operators often come in pairs
  T  operator() (unsigned row, unsigned col) const;  // Subscript operators often come in pairs
  // ...
 ~Matrix();                              // Destructor
  Matrix(const Matrix& m);               // Copy constructor
  Matrix& operator= (const Matrix& m);   // Assignment operator
  void print();
  // ...
private:
  unsigned rows_, cols_;
  T* data_;
};
template<typename T>
Matrix<T>::Matrix(unsigned rows, unsigned cols)
  : rows_ (rows)
  , cols_ (cols)
//, data_ ← initialized below after the if...throw statement
{
  if (rows == 0 || cols == 0)
    throw std::invalid_argument("Matrix constructor was given 0 size");
  data_ = new T[rows * cols]();
}
template<typename T>
Matrix<T>::~Matrix()
{
  delete[] data_;
}
template<typename T>
T& Matrix<T>::operator() (unsigned row, unsigned col)
{
  if (row >= rows_ || col >= cols_)
    throw std::out_of_range("Matrix subscript out of bounds");
  return data_[cols_*row + col];
}
template<typename T>
T Matrix<T>::operator() (unsigned row, unsigned col) const
{
  if (row >= rows_ || col >= cols_)
    throw std::out_of_range("const Matrix subscript out of bounds");
  return data_[cols_*row + col];
}
template<typename T>
void Matrix<T>::print() {
    std::cout << "shape: " << rows_ << ", " << cols_ << "\n";
    for (unsigned int i = 0; i < rows_; ++i) {
        for (unsigned int j = 0; j < cols_; ++j) {
            std::cout << (*this)(i, j) << " ";
        }
        std::cout << "\n";
    }
    std::cout << "---\n";
}

int main()
{
  Matrix<double> m(10,10);
  m(5,8) = 106.15;
  m.print();
}

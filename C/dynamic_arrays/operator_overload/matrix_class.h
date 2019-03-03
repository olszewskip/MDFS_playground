#pragma once

#include <iostream>
#include <stdexcept>

template <typename T>
class Matrix {
 public:

  Matrix(unsigned dim0, unsigned dim1);
  T& operator()(unsigned i, unsigned j);
  T operator()(unsigned i, unsigned j) const;

  Matrix(unsigned dim0, unsigned dim1, unsigned dim2);
  T& operator()(unsigned i, unsigned j, unsigned k);
  T operator()(unsigned i, unsigned j, unsigned k) const;

  Matrix(unsigned dim0, unsigned dim1, unsigned dim2, unsigned dim3);
  T& operator()(unsigned i, unsigned j, unsigned k, unsigned l);
  T operator()(unsigned i, unsigned j, unsigned k, unsigned l) const;

  Matrix(unsigned dim0, unsigned dim1, unsigned dim2, unsigned dim3,
         unsigned dim4);
  T& operator()(unsigned i, unsigned j, unsigned k, unsigned l, unsigned m);
  T operator()(unsigned i, unsigned j, unsigned k, unsigned l,
               unsigned m) const;

  Matrix(unsigned dim0, unsigned dim1, unsigned dim2, unsigned dim3,
         unsigned dim4, unsigned dim5);
  T& operator()(unsigned i, unsigned j, unsigned k, unsigned l, unsigned m,
                unsigned n);
  T operator()(unsigned i, unsigned j, unsigned k, unsigned l, unsigned m,
               unsigned n) const;

  ~Matrix();                           // Destructor
  Matrix(const Matrix& m);             // Copy constructor
  Matrix& operator=(const Matrix& m);  // Assignment operator
  void print();

 private:
  unsigned num_dims_;
  unsigned dim0_, dim1_, dim2_, dim3_, dim4_, dim5_;
  T* data_;
};

template <typename T>
Matrix<T>::Matrix(unsigned dim0, unsigned dim1)
    : num_dims_(2),
      dim0_(dim0),
      dim1_(dim1)
//, data_ ← initialized below after the if...throw statement
{
  if (dim0 == 0 || dim1 == 0)
    throw std::invalid_argument("Matrix constructor was given 0 size");
  data_ = new T[dim0 * dim1]();
}

template <typename T>
Matrix<T>::Matrix(unsigned dim0, unsigned dim1, unsigned dim2)
    : num_dims_(3),
      dim0_(dim0),
      dim1_(dim1),
      dim2_(dim2)
//, data_ ← initialized below after the if...throw statement
{
  if (dim0 == 0 || dim1 == 0 || dim2 == 0)
    throw std::invalid_argument("Matrix constructor was given 0 size");
  data_ = new T[dim0 * dim1 * dim2]();
}

template <typename T>
Matrix<T>::Matrix(unsigned dim0, unsigned dim1, unsigned dim2, unsigned dim3)
    : num_dims_(4),
      dim0_(dim0),
      dim1_(dim1),
      dim2_(dim2),
      dim3_(dim3)
//, data_ ← initialized below after the if...throw statement
{
  if (dim0 == 0 || dim1 == 0 || dim2 == 0 || dim3 == 0)
    throw std::invalid_argument("Matrix constructor was given 0 size");
  data_ = new T[dim0 * dim1 * dim2 * dim3]();
}

template <typename T>
Matrix<T>::Matrix(unsigned dim0, unsigned dim1, unsigned dim2, unsigned dim3,
                  unsigned dim4)
    : num_dims_(5),
      dim0_(dim0),
      dim1_(dim1),
      dim2_(dim2),
      dim3_(dim3),
      dim4_(dim4)
//, data_ ← initialized below after the if...throw statement
{
  if (dim0 == 0 || dim1 == 0 || dim2 == 0 || dim3 == 0 || dim4 == 0)
    throw std::invalid_argument("Matrix constructor was given 0 size");
  data_ = new T[dim0 * dim1 * dim2 * dim3 * dim4]();
}

template <typename T>
Matrix<T>::Matrix(unsigned dim0, unsigned dim1, unsigned dim2, unsigned dim3,
                  unsigned dim4, unsigned dim5)
    : num_dims_(6),
      dim0_(dim0),
      dim1_(dim1),
      dim2_(dim2),
      dim3_(dim3),
      dim4_(dim4),
      dim5_(dim5)
//, data_ ← initialized below after the if...throw statement
{
  if (dim0 == 0 || dim1 == 0 || dim2 == 0 || dim3 == 0 || dim4 == 0)
    throw std::invalid_argument("Matrix constructor was given 0 size");
  data_ = new T[dim0 * dim1 * dim2 * dim3 * dim4 * dim5]();
}


template <typename T>
Matrix<T>::~Matrix() {
  delete[] data_;
}

template <typename T>
T& Matrix<T>::operator()(unsigned i, unsigned j) {
  if (num_dims_ != 2)
    throw std::invalid_argument(
        "Number of indeces doesn't match the matrix dimension");
  if (i >= dim0_ || j >= dim1_)
    throw std::out_of_range("const Matrix subscript out of bounds");
  return data_[dim1_ * i + j];
}

template <typename T>
T Matrix<T>::operator()(unsigned i, unsigned j) const {
  if (num_dims_ != 2)
    throw std::invalid_argument(
        "Number of indeces doesn't match the matrix dimension");
  if (i >= dim0_ || j >= dim1_)
    throw std::out_of_range("const Matrix subscript out of bounds");
  return data_[dim1_ * i + j];
}

template <typename T>
T& Matrix<T>::operator()(unsigned i, unsigned j, unsigned k) {
  if (num_dims_ != 3)
    throw std::invalid_argument(
        "Number of indeces doesn't match the matrix dimension");
  if (i >= dim0_ || j >= dim1_ || k >= dim2_)
    throw std::out_of_range("const Matrix subscript out of bounds");
  return data_[dim1_ * dim2_ * i + dim2_ * j + k];
}

template <typename T>
T Matrix<T>::operator()(unsigned i, unsigned j, unsigned k) const {
  if (num_dims_ != 3)
    throw std::invalid_argument(
        "Number of indeces doesn't match the matrix dimension");
  if (i >= dim0_ || j >= dim1_ || k >= dim2_)
    throw std::out_of_range("const Matrix subscript out of bounds");
  return data_[dim1_ * dim2_ * i + dim2_ * j + k];
}

template <typename T>
T& Matrix<T>::operator()(unsigned i, unsigned j, unsigned k, unsigned l) {
  if (num_dims_ != 4)
    throw std::invalid_argument(
        "Number of indeces doesn't match the matrix dimension");
  if (i >= dim0_ || j >= dim1_ || k >= dim2_ || l >= dim3_)
    throw std::out_of_range("const Matrix subscript out of bounds");
  return data_[dim1_ * dim2_ * dim3_ * i + dim2_ * dim3_ * j + dim3_ * k + l];
}

template <typename T>
T Matrix<T>::operator()(unsigned i, unsigned j, unsigned k, unsigned l) const {
  if (num_dims_ != 4)
    throw std::invalid_argument(
        "Number of indeces doesn't match the matrix dimension");
  if (i >= dim0_ || j >= dim1_ || k >= dim2_ || l >= dim3_)
    throw std::out_of_range("const Matrix subscript out of bounds");
  return data_[dim1_ * dim2_ * dim3_ * i + dim2_ * dim3_ * j + dim3_ * k + l];
}

template <typename T>
T& Matrix<T>::operator()(unsigned i, unsigned j, unsigned k, unsigned l, unsigned m) {
  if (num_dims_ != 5)
    throw std::invalid_argument(
        "Number of indeces doesn't match the matrix dimension");
  if (i >= dim0_ || j >= dim1_ || k >= dim2_ || l >= dim3_ || m >= dim4_)
    throw std::out_of_range("const Matrix subscript out of bounds");
  return data_[dim1_ * dim2_ * dim3_ * dim4_* i + dim2_ * dim3_ * dim4_ * j + dim3_ * dim4_ * k + dim4_ * l];
}

template <typename T>
T Matrix<T>::operator()(unsigned i, unsigned j, unsigned k, unsigned l, unsigned m) const {
  if (num_dims_ != 5)
    throw std::invalid_argument(
        "Number of indeces doesn't match the matrix dimension");
  if (i >= dim0_ || j >= dim1_ || k >= dim2_ || l >= dim3_)
    throw std::out_of_range("const Matrix subscript out of bounds");
  return data_[dim1_ * dim2_ * dim3_ * dim4_* i + dim2_ * dim3_ * dim4_ * j + dim3_ * dim4_ * k + dim4_ * l];
}

template <typename T>
void Matrix<T>::print() {
  switch (num_dims_) {
    case 2:
      std::cout << "|" << __PRETTY_FUNCTION__ << "\n";
      std::cout << "|dimensions: " << num_dims_ << "; shape: " << dim0_ << ", "
                << dim1_ << "\n";
      for (unsigned int i = 0; i < dim0_; ++i) {
        for (unsigned int j = 0; j < dim1_; ++j) {
          std::cout << (*this)(i, j) << " ";
        }
        std::cout << "\n";
      }
      std::cout << "|------\n";
      break;

    case 3:
      std::cout << "|" << __PRETTY_FUNCTION__ << "\n";
      std::cout << "|dimensions: " << num_dims_ << "; shape: " << dim0_ << ", "
                << dim1_ << ", " << dim2_ << "\n";
      for (unsigned int i = 0; i < dim0_; ++i) {
        for (unsigned int j = 0; j < dim1_; ++j) {
          for (unsigned int k = 0; k < dim2_; ++k) {
            std::cout << (*this)(i, j, k) << " ";
          }
          std::cout << "\n";
        }
        std::cout << "\n";
      }
      std::cout << "|------\n";
      break;

    default:
      std::cout << "|" << __PRETTY_FUNCTION__ << "\n";
      std::cout << "|dimensions: " << num_dims_ << "\n";
      std::cout << "NOT IMPLEMENTED"
                << "\n";
  }
}

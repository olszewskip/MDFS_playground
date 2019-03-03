#include "matrix_class.h"

int main()
{
  Matrix<double> m1(2, 3, 4);
  m1(1, 2, 3) = 12.134;
  m1.print();

  Matrix<int> m2(2, 3, 4);
  m2(1, 2, 3) = 12134;
  m2.print();

  Matrix<int> m3(2, 3, 4, 5);
  m3(1, 2, 3, 4) = 121345;
  m3.print();

}

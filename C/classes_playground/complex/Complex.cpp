#include <iostream>
#include <cmath>
#include <math.h>

#include "Complex.h"

Complex Complex::operator-() const {
   return Complex(-m_re, -m_im);
}
Complex Complex::operator+(const Complex &z) {
   return Complex(m_re + z.m_re, m_im + z.m_im);
}
Complex operator+(const double &x, const Complex &z) {
   return Complex(x) + z;
}
Complex Complex::operator-(const Complex &z) {
   return Complex(m_re - z.m_re, m_im - z.m_im);
}
Complex operator-(const double &x, const Complex &z) {
   return -Complex(x) + z;
}
Complex Complex::operator*(const Complex &z) {
   return Complex(m_re * z.m_re - m_im * z.m_im, m_re * z.m_im + m_im * z.m_re);
}
Complex operator*(const double &x, const Complex &z) {
   return Complex(x) * z;
}
Complex Complex::operator/(const Complex &z) {
   double denominator = z.m_re * z.m_re + z.m_im * z.m_im;
   return Complex((m_re*z.m_re + m_im*z.m_im)/denominator, (-m_re*z.m_im + m_im*z.m_re)/denominator);
}
Complex operator/(const double &x, const Complex &z) {
   return Complex(x) / z;
}

//Complex imag(const double &y) {
//   return Complex(0, y);
//}

double Complex::radius() {
   std::sqrt(m_re * m_re + m_im * m_im);
}
double Complex::phase(double offset) {
   double phi = atan(m_im / m_re);
   if (m_re < 0.0)
      phi += M_PI;
   else if (m_im < 0.0)
      phi += 2*M_PI;
   return phi + offset;
}
std::ostream& operator<<(std::ostream &out, const Complex &z) {
   return out << z.m_re << "+" << z.m_im <<"i";
}

   

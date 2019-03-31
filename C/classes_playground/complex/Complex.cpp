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

Complex imag(const double &y) {
   return Complex(0, y);
}

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

int main(){
   Complex z;
   std::cout << z << std::endl;
   Complex z1 {1.0};
   std::cout << z1 << std::endl;
   z1 = 2.0;
   std::cout << z1 << std::endl;
   z1 = {3.0, 4.0};
   std::cout << z1 << std::endl;
   z1.setRe(5.0);
   z1.setIm(6.0);
   std::cout << z1 << std::endl;

   Complex z2(7.0, 8.0);
   std::cout << z2 << std::endl;
   std::cout << z2 + 1 << std::endl;
   std::cout << z2 + imag(1) << std::endl;
   
   std::cout << (1 + imag(2)) * Complex(0,3) / imag(0.5) << std::endl;
   std::cout << Complex(-4,-3).radius() << " " << imag(-1).phase() << " " << -imag(1).phase() << std::endl;
   return 0;
}



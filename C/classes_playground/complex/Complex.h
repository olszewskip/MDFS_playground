#ifndef COMPLEX
#define COMPLEX

#include <iostream>

class Complex {
   private:
      double m_re = 0.0;
      double m_im = 0.0;
   public:
      Complex() {}
      Complex(double re) : m_re {re} {}
      Complex(double re, double im) : m_re {re}, m_im {im} {}
      double getRe() { return m_re; }
      double getIm() { return m_im; }
      void setRe(double re) { m_re = re; }
      void setIm(double im) { m_im = im; }
      
      Complex operator-() const;
      Complex operator+(const Complex &z);
      Complex operator-(const Complex &z);
      Complex operator*(const Complex &z);
      Complex operator/(const Complex &z);
      friend Complex operator*(const double &x, const Complex &z);
      friend Complex operator-(const double &x, const Complex &z);
      friend Complex operator+(const double &x, const Complex &z);
      friend Complex operator/(const double &x, const Complex &z);
      
      //friend Complex imag(const double &y);
      
      double radius();
      double phase(double offset = 0.0);
      friend std::ostream& operator<<(std::ostream &out, const Complex &z);
};

#endif


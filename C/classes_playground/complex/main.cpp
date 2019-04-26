#include <iostream>
#include <cmath>
#include <math.h>

#include "Complex.h"

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
   //std::cout << z2 + imag(1) << std::endl;
   
   //std::cout << (1 + imag(2)) * Complex(0,3) / imag(0.5) << std::endl;
   std::cout << Complex(-4,-3).radius() << " ";// << imag(-1).phase() << " " << -imag(1).phase() << std::endl;
   return 0;
}



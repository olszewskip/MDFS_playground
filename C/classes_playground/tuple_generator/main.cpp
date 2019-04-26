// uncomment to disable assert()
// #define NDEBUG
#include <cassert>
#include <stdlib.h>
#include <iostream>

#define kDim 4 


class tuple{
   /* Simple data structure that holds a sequence
      of kDim strictly increasing indeces.
      It starts with (0, 1, ..., kDim-1),
      and it can be taken one step forward
      with the method `up`, eg.
      (0, 1, 2, 3) -> (0, 1, 2, 4)
      (0, 1, 2, 4) -> (0, 1, 3, 4),
      assuming that kDim=4 and that
      tuple(n) was called with n > 4.
      The field `exhausted` states, whether
      the last possible sequence, i.e
      (n-kDim, ..., n-2, n-1) was reached.
   */
   int m_n;
   bool m_with_diagonals;
   public:
   int indeces[kDim];
   bool exhausted = false;
   tuple(int n, bool with_diagonals) : m_n(n), m_with_diagonals(with_diagonals) {
      assert(kDim <= n);
      for(int i = 0; i < kDim; i++)
         indeces[i] = 0 + !with_diagonals * i;
   }
   void up() {
      assert(!exhausted);
      int step = 1;
      for(int i = 0; i < kDim - 1; i++){
         div_t division = div(indeces[i] + step, indeces[i+1] + m_with_diagonals);
         step = division.quot;
         indeces[i] = division.rem + !m_with_diagonals * step * i;
      }
      div_t division = div(indeces[kDim - 1] + step, m_n + m_with_diagonals);
      indeces[kDim - 1] = division.rem;
      exhausted = division.quot;
   }
   friend std::ostream& operator<<(std::ostream &out, const tuple &t) {
      out << "(";
      for (int i = 0; i < kDim-1; i++)
         out << t.indeces[i] << ", ";
      out << t.indeces[kDim-1] << ")";
      return out;
   }
};



int main() {
   
   tuple my_tuple = tuple(5, false);
   while (!my_tuple.exhausted) {
      std::cout << my_tuple << std::endl;
      my_tuple.up();
   }
   
   return 0;
}


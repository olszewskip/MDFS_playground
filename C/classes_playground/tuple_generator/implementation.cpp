#include <stdlib.h>
//#include "interface.h"
#define kDim 3

class tuple{
   int m_n;
   public:
      int indeces[kDim];
      
      tuple(int n) : m_n(n) {
         for(int i = 0; i < kDim - 1; i++)
            indeces[i] = i;
      }
      void up();
      bool is_outside() { return indeces[kDim - 1] == m_n; }
};


void tuple::up(){
   indeces[0] + 1;
   for(int i = 0; i < kDim - 1; i++){
      div_t dv = div(indeces[i], indeces[i+1]);
      indeces[i] = div.rem;
      indeces[i+1] += div.quot; 
   }
}


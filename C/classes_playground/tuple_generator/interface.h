#define kDim 3

class tuple{
   int n;
   public:
      int indeces[kDim];
      
      tuple(int n);
      void up();
      bool is_outside() { return indeces[kDim - 1] == n; }
};


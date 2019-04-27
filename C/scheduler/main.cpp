// uncomment to disable assert()
// #define NDEBUG
#include <cassert>
#include <stdlib.h>
#include <iostream>
#include <sstream>

const int kDim = 5;


template <typename T>
class Matrix2D {
   // Selfexplanatory
   private:
      T* data;
   public:
      unsigned dim_0;
      unsigned dim_1;
      Matrix2D(unsigned dim_0_arg, unsigned dim_1_arg) :
         dim_0(dim_0_arg),
         dim_1(dim_1_arg) {
         assert(dim_0 > 0 & dim_1 > 1);
         // zero-initialized
         data = new T[dim_0_arg * dim_1_arg]();
      }
      ~Matrix2D() {
         delete[] data;
      }
      Matrix2D(const Matrix2D&);
      Matrix2D& operator=(const Matrix2D&);
      T& operator()(unsigned i, unsigned j) {
         assert(i < dim_0 & j < dim_1);
         return data[dim_1 * i + j];
      }
      T operator()(unsigned i, unsigned j) const {
         assert(i < dim_0 & j < dim_1);
         return data[dim_1 * i + j];
      }
      void print() {
         std::cout << "|" << __PRETTY_FUNCTION__ << "\n";
         std::cout << "|shape: " << dim_0 << ", " << dim_1 << "\n";
            for (int i = 0; i < dim_0; ++i) {
               for (int j = 0; j < dim_1; ++j)
                  std::cout << (*this)(i, j) << " ";
               std::cout << "\n";
            }
      }
};


class Indices {
   // Abstraction for a finite sequence of indices,
   // where each 'index' is a multi-index i.e. it 
   // is itself a fixed-sized sequence of integers.
   // The latter are by assumption stored in a fixed
   // buffer, and this class serves to modify them.
   public:
      // By the sequence being exhausted we mean
      // that it has been taken beyond its last element.
      // The flag is not yet set when the indices
      // represent the last element of the sequence.
      virtual bool get_exhausted() = 0;
      virtual void set_exhausted(bool) = 0;
      // pointer to a buffer with actual integers
      int* index = NULL;
      // connect pointer to indices
      virtual void use_buff(int* const) = 0;
      // number of indices (lenght of buffer)
      int k;
      // move the index one step forward in the sequence
      virtual void up() = 0; 
      // change the index back to the first element of the sequence
      virtual void reset() = 0;
      virtual std::string to_str() = 0;
      void print(unsigned limit = -1) {
         assert(index);
         unsigned count = 0;
         while(!get_exhausted()) {
            if (count == limit)
               break;
            count++;
            std::cout << to_str() << std::endl;
            up();
         }
         reset();
      }
};


class Indices_triangle : public Indices {
   // Each index in the sequence is an increasing sequence,
   // e.g. (3,4,6,8) for k = 4 and n > 8
   protected:
      // is the whole sequence exhauseted
      bool exhausted;
      // strictly upper bound of values of all indices
      int n;
      // is the index-sequence not strictly increasing
      bool with_diag;
   public:
      // constructor
      Indices_triangle (int n_arg, int k_arg, bool with_diag_arg) :
         n(n_arg),
         with_diag(with_diag_arg) {
         assert(k_arg > 0); 
         k = k_arg;
         if (!with_diag_arg & n_arg < k_arg) {
            k = 0;
            exhausted = true;
         }
      }
      void use_buff(int* const buff) {
         index = buff;
         reset();
      }
      bool get_exhausted() { return exhausted; }
      void set_exhausted(bool exhausted_var) { exhausted = exhausted_var; }
      void up() {
         assert(k > 0);
         bool increment = 1;
         for (int i = 0; i < k - 1; i++) {
            div_t division = div(index[i] + increment, index[i+1] + with_diag);
            increment = division.quot;
            index[i] = division.rem + !with_diag * increment * i;
         }
         div_t division = div(index[k - 1] + increment, n);
         index[k - 1] = division.rem;
         exhausted = division.quot;
      }
      void reset() {
         if (!with_diag & n < k) {
            k = 0;
            exhausted = true;
            return;
         }
         for(int i = 0; i < k; i++) {
            index[i] = 0 + !with_diag * i;
            exhausted = false;
         }
      }
      std::string to_str() {
         std::stringstream ss;
         ss << "(";
         for (int i = 0; i < k - 1; i++)
            ss << index[i] << ", ";
         ss << index[k - 1] << ")";
         return ss.str();
      }
};


class Indices_product : public Indices {   
   // Cartesian product of two sequences.
   // Again results in a sequence (the order
   // is "row-wise", i.e. second index changes
   // faster).
   protected:
      Indices& ind_L;
      Indices& ind_R;
   public:
      int semicolon_index;
      Indices_product(Indices& ind_L_arg, Indices& ind_R_arg) :
         ind_L(ind_L_arg),
         ind_R(ind_R_arg),
         semicolon_index(ind_L_arg.k) {
            k = ind_L_arg.k + ind_R_arg.k;
      }
      void use_buff(int* const buff) {
         index = buff;
         ind_L.use_buff(buff);
         ind_R.use_buff(buff + semicolon_index);
      }
      bool get_exhausted() {
         return ind_L.get_exhausted();
      }
      void set_exhausted(bool exhausted_var) {
         ind_L.set_exhausted(exhausted_var);
      }
      void up() {
         ind_R.up();
         if (ind_R.get_exhausted()) {
            ind_R.reset();
            ind_L.up();
         }
      }
      void reset() {
         ind_L.reset();
         ind_R.reset();
      }
      std::string to_str() {
         std::stringstream ss;
         ss << ind_L.to_str() << ind_R.to_str();
         return ss.str();
      }
};



class Indices_product_with_semicolon : public Indices_product {
   // It uses a buffer that is assumed to be one integer
   // longer than for the Indices_product class. Value
   // stored in this additional memory is the length of
   // first factor's buffer. 
   public:
      Indices_product_with_semicolon(Indices& ind_L_arg, Indices& ind_R_arg) :
         Indices_product(ind_L_arg, ind_R_arg) {}
      void use_buff(int* const buff) {
         index = buff;
         ind_L.use_buff(buff);
         ind_R.use_buff(buff + ind_L.k);
         // Append another value to the buffer, namely
         // the border-position within the buffer
         // between indices of the two factor sequences.
         index[k] = ind_L.k;
      }
      std::string to_str() {
         std::stringstream ss;
         ss << Indices_product::to_str() << " " << index[k];
         return ss.str();
      }
};


class Indices_sum : public Indices {
   // Concatenation of two sequences.
   // (Assumes that both represent multi-indices
   // of the same lenght (k).)
   // Meaning that, when we step over the sum,
   // we first exhaust the first component sequence,
   // and then step through the second one.
   protected:
      Indices& ind_L;
      Indices& ind_R;
      bool using_L;
   public:
      Indices_sum(Indices &ind_L_arg, Indices &ind_R_arg) :
         ind_L(ind_L_arg),
         ind_R(ind_R_arg),
         using_L(true) {
            assert(ind_L_arg.k == ind_R_arg.k);
            k = ind_L_arg.k;
      }
      void use_buff(int* const buff) {
         index = buff;
         ind_L.use_buff(buff);
         ind_R.set_exhausted(false);
      }
      bool get_exhausted() { return ind_R.get_exhausted(); }
      void set_exhausted(bool exhausted_var) { ind_R.set_exhausted(exhausted_var); }
      void up() {
         if (using_L) {
            ind_L.up();
            if(ind_L.get_exhausted()) {
               using_L = false;
               ind_R.use_buff(index);
            }
         }
         else {
            ind_R.up();
         }
      }
      void reset() {
         ind_L.use_buff(index);
         using_L = true;
         set_exhausted(false);
      }
      std::string to_str() {
         std::stringstream ss;
         ss << (using_L ? ind_L.to_str() : ind_R.to_str());
         return ss.str();
      }
      //friend std::ostream& operator<<(std::ostream&, const Indices_sum&);
};

//std::ostream& operator<<(std::ostream& out, const Indices_sum& indices) { 
//   out << indices.using_L ? indices.ind_L : indices.ind_R;
//   return out;
//}

//class Indices_sum_of_products_with_semicolon : public Indices_sum {
//   protected:
//      int semicolon_index_L;
//      int semicolon_index_R;
//   public:
//      Indices_sum_of_products_with_semicolon(Indices_product_with_semicolon& ind_L_arg, Indices_product_with_semicolon ind_R_arg) :
//         Indices_sum(ind_L_arg, ind_R_arg),
//         semicolon_index_L(ind_L_arg.ind_L.k),
//         semicolon_index_R(ind_R_arg.ind_L.k) {}
//      friend std::ostream& operator<<(std::ostream&, const Indices_sum_of_products_with_semicolon&);
//};
//
//std::ostream& operator<<(std::ostream& out, const Indices_sum_of_products_with_semicolon& indices) {
//   out << "(";
//   int semicolon_index = indices.using_R ? semicolon_index_R : semicolon_index_L;
//   for (int i = 0; i < semicolon_index - 1; i++)
//      out << indices.index[i] << ", ";
//   out << indices.index[semicolon_index - 1] << "; ";
//   for (int i = semicolon_index; i < indices.k - 1; i++)
//      out << indices.index[i] << ", ";
//   out << indices.index[indices.k - 1] << ") ";
//   out << indices.index[indices.k];
//   return out;
//}



Indices_product operator*(Indices &ind_L, Indices &ind_R) { return Indices_product(ind_L, ind_R); }
Indices_product_with_semicolon operator&(Indices &ind_L, Indices &ind_R) { return Indices_product_with_semicolon(ind_L, ind_R); }
Indices_sum operator+(Indices &ind_L, Indices &ind_R) { return Indices_sum(ind_L, ind_R); }
//Indices_sum operator+(Indices &ind_L, Indices &ind_R) { return Indices_sum(ind_L, ind_R); }


void index2columns(int index[], int columns[], int tile_width) {
   for (int i = 0; i < kDim; i+2) {
      columns[i] = index[i] * tile_width;
      columns[i+1] = (index[i] + 1) * tile_width;
   }
} 
//void test(const int k) {

int main() {

    
   Indices_triangle my_indices_1 = Indices_triangle(5, kDim - 3, true);
   Indices_triangle my_indices_2 = Indices_triangle(4, 3, true);
   Indices_triangle my_indices_3 = Indices_triangle(3, kDim - 1, true); 
   Indices_triangle my_indices_4 = Indices_triangle(2, 1, true); 
   
   Indices_product product_12 = my_indices_1 * my_indices_2;
   Indices_product_with_semicolon product_34 = my_indices_3 & my_indices_4;
   Indices_sum sum_1234 = product_12 + product_34;
   
   int indices_buff_1234[kDim + 1];
   sum_1234.use_buff(indices_buff_1234);
   //product_12.use_buff(indices_buff_1234);
   //std::cout << sum_1234 << std::endl;
   //std::cout << sum_1234.to_str() << std::endl;
   //product_12.print();
   sum_1234.print();
   


   
   //Indices_triangle gpu_scheduler = Indices_triangle(
   

   return 0;
}


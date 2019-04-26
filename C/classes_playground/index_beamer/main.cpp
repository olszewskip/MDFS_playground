// uncomment to disable assert()
// #define NDEBUG
#include <cassert>
#include <stdlib.h>
#include <iostream>

#define kDim 5

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
      int *indices_ptr = NULL;
      // connect pointer to indices
      virtual void use_buff(int* const) = 0;
      // number of indices (lenght of buffer)
      int k;
      // move the index one step forward in the sequence
      virtual void up() = 0; 
      // change the index back to the first element of the sequence
      virtual void reset() = 0;
      // various fields glued together to form
      // a contingous array inteded to be send
      // to another galaxy
      int *mpi_buff_ptr = NULL;
      int mpi_buff_size;
      // make the mpi_buff point at something useful
      virtual void prepare_mpi_buff() {
         mpi_buff_ptr = indices_ptr;
         mpi_buff_size = k;
      }
      // << operator for the current index
      friend std::ostream& operator<<(std::ostream &out, const Indices &ind) {
         out << "(";
         for (int i = 0; i < ind.k - 1; i++)
            out << ind.indices_ptr[i] << ", ";
         out << ind.indices_ptr[ind.k - 1] << ")";
         return out;
      }
      void print(unsigned int limit = -1) {
         unsigned int count = 0;
         while(!this->get_exhausted()) {
            if (count == limit)
               break;
            count++;
            std::cout << *this << std::endl;
            this->up();
         }
         this->reset();
      }
};

class Indices_triangle : public Indices {
   // Each index in the sequence is an increasing sequence,
   // e.g. (3,4,6,8) for k = 4 and n > 8
   private:
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
         indices_ptr = buff;
         reset();
      }
      bool get_exhausted() { return exhausted; }
      void set_exhausted(bool exhausted_var) { exhausted = exhausted_var; }
      void up() {
         assert(k > 0);
         bool increment = 1;
         for (int i = 0; i < k - 1; i++) {
            div_t division = div(indices_ptr[i] + increment, indices_ptr[i+1] + with_diag);
            increment = division.quot;
            indices_ptr[i] = division.rem + !with_diag * increment * i;
         }
         div_t division = div(indices_ptr[k - 1] + increment, n);
         indices_ptr[k - 1] = division.rem;
         exhausted = division.quot;
      }
      void reset() {
         if (!with_diag & n < k) {
            k = 0;
            exhausted = true;
            return;
         }
         for(int i = 0; i < k; i++) {
            indices_ptr[i] = 0 + !with_diag * i;
            exhausted = false;
         }
      }
};

class Indices_product : public Indices {   
   // Cartesian product of two sequences.
   // Again results in a sequence (the order
   // is "row-wise", i.e. second index changes
   // faster).
   private:
      Indices& ind_L;
      Indices& ind_R;
   public:
      Indices_product(Indices& ind_L_arg, Indices& ind_R_arg) :
         ind_L(ind_L_arg),
         ind_R(ind_R_arg) {
            k = ind_L_arg.k + ind_R_arg.k;
      }
      void use_buff(int* const buff) {
         indices_ptr = buff;
         ind_L.use_buff(buff);
         ind_R.use_buff(buff + ind_L.k);
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
};

class Indices_sum : public Indices {
   // Concatenation of two sequences.
   // (Assumes that both represent multi-indices
   // of the same lenght (k).)
   // Meaning that, when we step over the sum,
   // we first exhaust the first component sequence,
   // and then step through the second one.
   private:
      Indices& ind_L;
      Indices& ind_R;
      bool using_R;
   public:
      Indices_sum(Indices &ind_L_arg, Indices &ind_R_arg) :
         ind_L(ind_L_arg),
         ind_R(ind_R_arg),
         using_R(false) {
            assert(ind_L_arg.k == ind_R_arg.k);
            k = ind_L_arg.k;
      }
      void use_buff(int* const buff) {
         indices_ptr = buff;
         ind_L.use_buff(buff);
         ind_R.set_exhausted(false);
      }
      bool get_exhausted() { return ind_R.get_exhausted(); }
      void set_exhausted(bool exhausted_var) { ind_R.set_exhausted(exhausted_var); }
      void up() {
         if (!using_R) {
            ind_L.up();
            if(ind_L.get_exhausted()) {
               using_R = true;
               ind_R.use_buff(indices_ptr);
            }
         }
         else {
            ind_R.up();
         }
      }
      void reset() {
         ind_L.use_buff(indices_ptr);
         using_R = false;
         set_exhausted(false);
      }
};


Indices_product operator*(Indices &ind_L, Indices &ind_R) { return Indices_product(ind_L, ind_R); }
Indices_sum operator+(Indices &ind_L, Indices &ind_R) { return Indices_sum(ind_L, ind_R); }


int main() {
   Indices_triangle my_indices_1 = Indices_triangle(5, kDim - 3, true);
   Indices_triangle my_indices_2 = Indices_triangle(4, 3, true);
   Indices_triangle my_indices_3 = Indices_triangle(3, kDim - 1, true); 
   Indices_triangle my_indices_4 = Indices_triangle(2, 1, true); 
   
   Indices_product product_12 = my_indices_1 * my_indices_2;
   Indices_product product_34 = my_indices_3 * my_indices_4;
   Indices_sum sum_1234 = product_12 + product_34;
   
   int indices_buff_1234[kDim];
   sum_1234.use_buff(indices_buff_1234);
   sum_1234.print();

   return 0;
}


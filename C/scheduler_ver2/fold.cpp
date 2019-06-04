// uncomment to disable assert()
// #define NDEBUG
#include <cassert>
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <vector>
#include <chrono>
#include <thread>
#include <fstream>

// the dimensionality of k-tuples
const int kDim = 3;

class Indices {
   // Abstraction for a finite sequence of indices,
   // where each 'index' is a multi-index i.e. it 
   // is itself a fixed-sized sequence of integers,
   // say e.g. ( (0,0), (0,1), (0,2), (1,1), ...).
   // The integers of a current multi-index are
   // held in a buffer. This class serves to modify
   // the buffer to simulate sequence traversal.
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
   // e.g. (3,4,6,8) for k = 4 and n > 8.
   // With k = 0, the object is supposed to behave
   // like Indices_product_id
   protected:
      // upper bound (noninclusive) of values of all indices
      int n;
      // lower bound (inclusive) of values of all indices
      int offset;
      // is the index-sequence not strictly increasing
      bool with_diag;
      // is the whole sequence exhausted
      bool exhausted = false;
   public:
      // constructor
      Indices_triangle (int k_arg,
                        int n_arg,
                        int offset_arg,
                        bool with_diag_arg) :
         n(n_arg),
         offset(offset_arg),
         with_diag(with_diag_arg) {
         assert(offset >= 0);
         assert(n_arg >= offset);
         if (!with_diag_arg & (n_arg - offset) < k_arg) {
            k = 0;
         }
         else {
            k = k_arg;
         }
      }
      void use_buff(int* const buff) {
         index = buff;
         reset();
      }
      bool get_exhausted() { return exhausted; }
      void set_exhausted(bool exhausted_arg) { exhausted = exhausted_arg; }
      void up() {
         // F-style traversal
         if (k > 0) {
            bool increment = 1;
            for (int i = 0; i < k - 1; i++) {
               div_t division = div(index[i] + increment, index[i+1] + with_diag);
               increment = division.quot;
               index[i] = division.rem + increment * (offset + !with_diag * i);
            }
            index[k - 1] += increment;
            exhausted = index[k - 1] / n;
         }
         else {
            exhausted = true;
         }
      }
      void up_() {
         // C-style traversal
         if (k > 0) {
            int position = k - 1;
            bool shift_left = (index[position] + 1) / n;
            while (shift_left) {
               if (!position) {
                  exhausted = true;
                  return;
               }
               --position;
               shift_left = (index[position] + !with_diag) / index[position + 1];
            }
            ++index[position];
            for (int i = position + 1; i < k; i++) {
               index[i] = index[i - 1] + !with_diag;
            }  
         }
         else {
            exhausted = true;
         }
      }
      void reset() {
         for(int i = 0; i < k; i++)
            index[i] = offset + !with_diag * i;
         exhausted = false;
      }
      std::string to_str() {
         std::stringstream ss;
         if (k > 0) {
            ss << "(";
            for (int i = 0; i < k - 1; i++)
                ss << index[i] << ", ";
            ss << index[k - 1] << ")";
         }
         else {
            ss << "";
         }
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

class Indices_product_id : public Indices {
   // Multiplication's identity.
   protected:
      bool exhausted = false;
   public:
      Indices_product_id() { k = 0; }
      bool get_exhausted() { return exhausted; }
      void set_exhausted(bool exhausted_arg) { exhausted = exhausted_arg; }
      void use_buff(int* const buff) {}
      void up() { exhausted = true; }
      void reset() { exhausted = false; }
      std::string to_str() { return ""; }
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
         if (ind_L.get_exhausted()) { up(); }
      }
      bool get_exhausted() {
         return ind_L.get_exhausted() & ind_R.get_exhausted();
      }
      void set_exhausted(bool exhausted_var) {
         ind_L.set_exhausted(exhausted_var);
         ind_R.set_exhausted(exhausted_var);
      }
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
};

class Indices_sum_id : public Indices {
   // Summation's identity.
   // Parametrized by k.
   protected:
      bool exhausted = true;
   public:
      Indices_sum_id(int k_arg) { k = k_arg; }
      bool get_exhausted() { return true; }
      void set_exhausted(bool exhausted_arg) {}
      void use_buff(int* const buff) {}
      void up() {}
      void reset() {}
      std::string to_str() { return ""; }
};

Indices_sum& fold(std::vector<Indices_triangle>& triangle_vect_1,
                  std::vector<Indices_triangle>& triangle_vect_2,
                  std::vector<Indices_product>& triangle_1_x_2,
                  std::vector<Indices_sum>& triangle_1_x_2_sum,
                  Indices_sum_id& zero,
                  int index_1_L, int index_1_R,
                  int index_2_L, int index_2_R,
                  int lower_limit,
                  int upper_limit) {
   // \sum_{i = lower_limit}^{upper_limit}
   //    (index_1_L : index_1_R)^(upper_limit - i)
   //     \times
   //    (index_2_L : index_2_R)^i
   assert(upper_limit >= lower_limit);
   int component_count = upper_limit - lower_limit + 1;
   triangle_vect_1.clear();
   triangle_vect_2.clear();
   triangle_1_x_2.clear();
   triangle_1_x_2_sum.clear();
   triangle_vect_1.reserve(component_count);
   triangle_vect_2.reserve(component_count);
   triangle_1_x_2.reserve(component_count);
   triangle_1_x_2_sum.reserve(component_count);
   for (int i = lower_limit; i <= upper_limit; i++) {
      triangle_vect_1.emplace_back(upper_limit - i,
                                   index_1_R,
                                   index_1_L,
                                   true);
      triangle_vect_2.emplace_back(i,
                                   index_2_R,
                                   index_2_L,
                                   true);
      triangle_1_x_2.emplace_back(triangle_vect_1.back(),
                                  triangle_vect_2.back());
   }
   triangle_1_x_2_sum.emplace_back(zero, triangle_1_x_2[0]);
   for (int i = 1; i < component_count; i++) {
      triangle_1_x_2_sum.emplace_back(
            triangle_1_x_2_sum.back(),
            triangle_1_x_2[i]);
   }
   return triangle_1_x_2_sum.back();
}


int main() {

   int R[7] = {0, 3, 6, 12, 15, 18, 21};
   
   std::vector<Indices_triangle> signif_fat;
   std::vector<Indices_triangle> signif_lean;
   std::vector<Indices_product> signif_fat_x_lean;
   std::vector<Indices_sum> signif_fat_x_lean_sum;
   Indices_sum_id zero(kDim);

   Indices_sum& sum = fold(signif_fat,
                           signif_lean,
                           signif_fat_x_lean,
                           signif_fat_x_lean_sum,
                           zero,
                           R[0], R[1],
                           R[1], R[2],
                           1,
                           kDim);
   
   int buff[kDim];
   sum.use_buff(buff);
   sum.print();

   return 0;
}

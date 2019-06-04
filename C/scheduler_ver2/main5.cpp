// uncomment to disable assert()
// #define NDEBUG
#include <cassert>
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <vector>


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
         // assert(k_arg > 0); 
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




int main() {

   int buff[kDim];

   Indices_triangle triangle_1(1, 10, 0, true);
   std::vector<Indices_triangle> vect_1;
   vect_1.reserve(4);  // <======================== comment out
   vect_1.emplace_back(1, 10, 0, true);

   Indices_triangle triangle_2(2, 12, 10, true);
   vect_1.emplace_back(2, 12, 10, true);

   Indices_product prod_1(triangle_1, triangle_2);
   Indices_product prod_2(vect_1[0], vect_1[1]);

   Indices_triangle triangle_3(0, 10, 0, true);
   vect_1.emplace_back(0, 10, 0, true);

   Indices_triangle triangle_4(3, 12, 10, true);
   vect_1.emplace_back(3, 12, 10, true);

   Indices_product prod_3(triangle_3, triangle_4);
   Indices_product prod_4(vect_1[2], vect_1[3]);
   
   std::cout << "prod_1 = triangle_1 * triangle_2 =" << std::endl;
   prod_1.use_buff(buff);
   prod_1.print();

   std::cout << "prod_3 = triangle_3 * triangle_4 =" << std::endl;
   prod_3.use_buff(buff);
   prod_3.print();

   std::cout << "prod_2 = vect_1[0] * vect_1[1] =" << std::endl;
   prod_2.use_buff(buff);
   prod_2.print();

   std::cout << "prod_4 = vect_1[2] * vect_1[3] =" << std::endl;
   prod_4.use_buff(buff);
   prod_4.print();

  return 0;
} 

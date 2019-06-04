#include <mpi.h>
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

template <typename T>
class Matrix2D {
   // Selfexplanatory.
   private:
      T* data;
   public:
      int dim_0;
      int dim_1;
      Matrix2D(int dim_0_arg, int dim_1_arg) :
         dim_0(dim_0_arg),
         dim_1(dim_1_arg) {
         assert(dim_0 > 0 & dim_1 > 0);
         // zero-initialized
         data = new T[dim_0_arg * dim_1_arg]();
      }
      ~Matrix2D() {
         delete[] data;
      }
      Matrix2D(const Matrix2D&);
      Matrix2D& operator=(const Matrix2D&);
      T& operator()(int i, int j) {
         assert(i < dim_0 & j < dim_1);
         return data[dim_1 * i + j];
      }
      T* get_data() { return data; }
      T operator()(int i, int j) const {
         assert(i < dim_0 & j < dim_1);
         return data[dim_1 * i + j];
      }
      void print() {
         //std::cout << "|" << __PRETTY_FUNCTION__ << "\n";
         //std::cout << "| shape: " << dim_0 << ", " << dim_1 << "\n";
         for (int i = 0; i < dim_0; ++i) {
            for (int j = 0; j < dim_1; ++j)
               std::cout << (*this)(i, j) << " ";
            std::cout << "\n";
         }
         std::cout << "\n";
      }
};


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

class Indices_triangle_with_semicolon : public Indices_triangle {
   // It uses a buffer that is assumed to be one integer
   // longer than for the Indices_triangle class. Value
   // stored in this additional memory is passed in the 
   // constructor.
   public:
      int semicolon_index;
      Indices_triangle_with_semicolon(int k_arg, int n_arg,
                                      int offset_arg,
                                      bool with_diag_arg,
                                      int semicolon_index_arg) :
         semicolon_index(semicolon_index_arg),
         Indices_triangle (k_arg, n_arg, offset_arg, with_diag_arg) {}
      void use_buff(int* const buff) {
         Indices_triangle::use_buff(buff);
         // Append 0.
         index[k] = semicolon_index;
      }
      std::string to_str() {
         std::stringstream ss;
         ss << Indices_triangle::to_str() << semicolon_index;
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


class Indices_product_with_semicolon : public Indices_product {
   // It uses a buffer that is assumed to be one integer
   // longer than for the Indices_product class. Value
   // stored in this additional memory is the length of
   // first factor's buffer. 
   public:
      Indices_product_with_semicolon(Indices& ind_L_arg,
                                     Indices& ind_R_arg) :
         Indices_product(ind_L_arg, ind_R_arg) {}
      void use_buff(int* const buff) {
         Indices_product::use_buff(buff);
         // Append another value to the buffer, namely
         // the border-position within the buffer
         // between indices of the two factor sequences.
         index[k] = semicolon_index; // ind_L.k;
      }
      std::string to_str() {
         std::stringstream ss;
         ss << Indices_product::to_str() << semicolon_index; //index[k];
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

Indices_product operator*(Indices &ind_L, Indices &ind_R) {
   return Indices_product(ind_L, ind_R);
}
Indices_product_with_semicolon operator^(Indices &ind_L, Indices &ind_R) {
   return Indices_product_with_semicolon(ind_L, ind_R);
}
Indices_sum operator+(Indices &ind_L, Indices &ind_R) {
   return Indices_sum(ind_L, ind_R);
}


// 1
int main1() {
   // Some printing for demonstration.

   Matrix2D<int> m(2, 3);
   m(1,1) = 999;
   m.print();
   
   const int kDim = 5;

   Indices_triangle my_indices_1(kDim - 3, 7, 2, true);
   Indices_triangle my_indices_2(3, 4, 0, true);
   Indices_triangle my_indices_3(kDim - 1, 6, 3, true); 
   Indices_triangle my_indices_4(1, 10, 4, true);
   Indices_product_id one;

   std::cout << one.to_str() << std::endl;

   Indices_product product_01 = one * my_indices_1;
   Indices_product product_12 = product_01 * my_indices_2;
   Indices_product_with_semicolon product_34 = my_indices_3 ^ my_indices_4;
   Indices_sum sum_1234 = product_12 + product_34; 
   Indices_sum_id zero(5);
   Indices_sum sum_12340 = sum_1234 + zero;
   int indices_buff_12340[kDim + 1];

   sum_12340.use_buff(indices_buff_12340);
   sum_12340.print();

}


Matrix2D<int> gpu_index2pretty_columns(int gpu_index[], int tile_width) {
   // Bijection between tile-index and a 2d matrix.
   // that represent the tile by specifying column-indices.
   // Each two-element row in the k x 2 matrix contains beginning
   // and end (exclusive) of a column-indices-range.
   Matrix2D<int> columns(kDim, 2);
   for (int i = 0; i < kDim; i++) {
      columns(i, 0) = tile_width * gpu_index[i];
      columns(i, 1) = tile_width * (gpu_index[i] + 1);
   }
   return columns;
}

void gpu_tile_index2column_ranges(int gpu_column_ranges[],
                                  int gpu_tile_index[],
                                  int tile_width) {
   // Bijection between tile-index and a 2d matrix (row-major).
   // that represent the tile by specifying column-indices.
   // Each two-element row in the k * 2 matrix contains beginning
   // and end (exclusive) of a column-indices-range.
   for (int i = 0; i < kDim; i++) {
      gpu_column_ranges[2 * i] = tile_width * gpu_tile_index[i];
      gpu_column_ranges[2 * i + 1] = tile_width * (gpu_tile_index[i] + 1);
   }
}


struct Index_count {
   // Two integers.
   int index;
   int count;
   Index_count(int index_arg, int count_arg) :
      index(index_arg), count(count_arg) {}
   friend std::ostream& operator<<(std::ostream &out,
                                   const Index_count &ic) {
      out << ic.index << ": " << ic.count << std::endl;
      return out;
   }
};

std::vector<Index_count> index_counts(int buff[], int buff_size) {
   // Returns a vector that contains structs
   // representing unique indices in buff.
   std::vector<Index_count> unique_indices;
   if (buff_size > 0) {
      for (int i = 0; i < buff_size; i++) {
         int &index = buff[i];
         bool new_ = true;
         for (auto& unique : unique_indices) {
            if (index == unique.index) {
               ++unique.count;
               new_ = false;
               break;
            }
         }
         if (new_) { unique_indices.emplace_back(index, 1); }
      }
   }
   return unique_indices;
}


void cpu_tile_index2column_ranges(int cpu_column_ranges[],
                                  int cpu_tile_index[],
                                  int tile_width,
                                  std::vector<int> cum_n) {
   int semicolon = cpu_tile_index[kDim];
   for (int i = 0; i < semicolon; i++) {
      cpu_column_ranges[2 * i] = tile_width * cpu_tile_index[i];
      cpu_column_ranges[2 * i + 1] = tile_width * (cpu_tile_index[i] + 1);
   }
   for (int i = semicolon; i < kDim; i++) {
      cpu_column_ranges[2 * i] = cum_n[cpu_tile_index[i]];
      cpu_column_ranges[2 * i + 1] = cum_n[cpu_tile_index[i] + 1];
   }
}

void cpu_tuple_index2columns(int cpu_columns[],
                             int cpu_tile_index[],
                             int cpu_tuple_index[],
                             int tile_width,
                             std::vector<int> cum_n) {
   // Fill a vector with column indices
   // that represent a single k-tuple.
   int semicolon = cpu_tile_index[kDim];
   for (int i = 0; i < semicolon; i++)
      cpu_columns[i] = tile_width * cpu_tile_index[i] + cpu_tuple_index[i];
   for (int i = semicolon; i < kDim; i++)
      cpu_columns[i] = cum_n[cpu_tile_index[i]] + cpu_tuple_index[i];
}


// 2
int main2() {
   // // number of observations (rows in the input data table)
   // int N = 1000;
   // width of a square tile
   int tile_width = 100;
   // number of contigous column-subsets,
   // with the assumption that all columns
   // in one subset have equal basket-count;
   // by convention, the first subset is strictly
   // composed of the square tiles
   int m = 3;
   // number of square tiles (in one dimension)
   int M = 9;
   // column-counts in the column-subsets
   std::vector<int> n(m);
   n[0] = M * tile_width;
   n[1] = 97;
   n[2] = 3;
   // cumulative-column-counts
   std::vector<int> cum_n(m + 1);
   cum_n[0] = 0;
   for (int i = 0; i < n.size(); i++) {
      cum_n[i + 1] = cum_n[i] + n[i];
   }
   // the number of "divisions"
   int divisions = 4;
   // bucket-counts in column-subsets
   std::vector<int> s(m);
   s[0] = divisions + 1;
   s[1] = divisions + 1;
   s[2] = 2;

   // GPU scheduler
   
   Indices_triangle gpu_scheduler = Indices_triangle(kDim, M, 0, true);
   int gpu_buff[kDim];
   gpu_scheduler.use_buff(gpu_buff);
   std::cout << "gpu" << std::endl;
   std::cout << " ----- \n";
   
   while (!gpu_scheduler.get_exhausted()) {
      // The gpu_scheduler's index buffer is (presumably)
      // sent over with mpi.
      // std::cout << gpu_scheduler.to_str() << std::endl;
      
      // For each such buffer we do the following.
      int* gpu_tile_index = gpu_scheduler.index;
      
      int gpu_columns[2 * kDim];
      gpu_tile_index2column_ranges(gpu_columns,
                             gpu_tile_index,
                             tile_width);
      for (int i = 0; i < kDim; i++) {
         std::cout << gpu_columns[2 * i] << " : " <<
                      gpu_columns[2 * i + 1] << "\n";
      }

      // Matrix2D<int> gpu_columns = gpu_index2pretty_columns(tile_index, tile_width);
      //gpu_columns.print();

      std::cout << " ----- \n";     
      gpu_scheduler.up();
   }
   
   // CPU scheduler

   // Now, what I would want is something like:
   //
   // cpu_scheduler = triangle_with_semicolon(kDim, m, 1)
   // for i = 1 : kDim
   //    cpu_scheduler += triangle(i, M) ^ triangle(kDim - i, m, 1)
   //
   // But I don't know how to make that work in C++.
   // So instead You get the following.
   Indices_triangle_with_semicolon triangle_1 = Indices_triangle_with_semicolon(kDim, m, 1, true, 0);
   std::vector<Indices_triangle> triangle_2;
   std::vector<Indices_triangle> triangle_3;
   std::vector<Indices_product_with_semicolon> product_1;
   std::vector<Indices_sum> sum_1;
   triangle_2.reserve(kDim - 1);
   triangle_3.reserve(kDim - 1);
   product_1.reserve(kDim - 1);
   sum_1.reserve(kDim - 1);

   for (int i = 0; i < kDim - 1; i++) {
      triangle_2.emplace_back(i + 1, M, 0, true);
      triangle_3.emplace_back(kDim - i - 1, m, 1, true);
      product_1.emplace_back(triangle_2[i] ^ triangle_3[i]);
   }
   sum_1.emplace_back(triangle_1 + product_1[0]);
   for (int i = 1; i < kDim - 1; i++)
      sum_1.emplace_back(sum_1[i - 1] + product_1[i]);
   
   Indices_sum& cpu_scheduler = sum_1.back();
   // I guess it will do.

   int cpu_buff[kDim + 1];
   cpu_scheduler.use_buff(cpu_buff);
   
   std::cout << "cpu" << std::endl;
   std::cout << " ----- \n";

   while(!cpu_scheduler.get_exhausted()) {
      // The cpu_scheduler's index buffer is (presumably)
      // sent over with mpi.
      // std::cout << cpu_scheduler.to_str() << std::endl;
      
      // For each such buffer we do the following.
      int* cpu_tile_index = cpu_scheduler.index;
      
      // Aggregate the numbers in tile_index.
      // The boundary between square and rectangle
      // indices is denoted by semicolon.
      int semicolon = cpu_tile_index[kDim];
      std::vector<Index_count> counts_square
                = index_counts(cpu_tile_index, semicolon);
      std::vector<Index_count> counts_rectangle
                = index_counts(cpu_tile_index + semicolon, kDim - semicolon);
      // for (auto& index_count : counts_square) { std::cout << index_count; }
      // for (auto& index_count : counts_rectangle) { std::cout << index_count; }

      // Prepare the tuple-scheduler
      std::vector<Indices_triangle> triangle_4;
      triangle_4.reserve(counts_square.size() + counts_rectangle.size());
      for (auto& ic : counts_square) {
         triangle_4.emplace_back(ic.count, tile_width, 0, false);
      }
      for (auto& ic : counts_rectangle) {
         triangle_4.emplace_back(ic.count, n[ic.index], 0, false);
      }
      std::vector<Indices_product> product_2;
      product_2.reserve(triangle_4.size());
      Indices_product_id id; 
      product_2.emplace_back(id * triangle_4[0]);
      for (int i = 1; i < triangle_4.size(); i++) {
         product_2.emplace_back(product_2.back() * triangle_4[i]);
      }
      Indices_product &cpu_tuple_scheduler = product_2.back();
      int cpu_tuple_index[kDim];
      cpu_tuple_scheduler.use_buff(cpu_tuple_index);
     
      // Print a first few k-tuples of columns.
      int few = 4; 
      for (int i = 0; i < few; i++)
      //while (!cpu_tuple_scheduler.get_exhausted())
      {
         int cpu_columns[kDim];
         cpu_tuple_index2columns(cpu_columns,
                                 cpu_tile_index,
                                 cpu_tuple_index,
                                 tile_width,
                                 cum_n);
         for (int i = 0; i < kDim; i++)
            std::cout << cpu_columns[i] << " ";
         std::cout << std::endl;
         
         cpu_tuple_scheduler.up();
         if (cpu_tuple_scheduler.get_exhausted()) { break; }
      }

      std::cout << " ----- \n";
      cpu_scheduler.up();
   }
}


int dummy_work(int *column_ranges, int rank) {
   std::this_thread::sleep_for(std::chrono::milliseconds(random() % 100));
   std::cout << rank << " column-ranges: ";
   for (int i = 0; i < kDim - 1; i++) {
      std::cout << column_ranges[2 * i] << ":" << column_ranges[2 * i + 1] <<", ";
   }
   std::cout << column_ranges[2 * kDim - 2] << ":" << column_ranges[2 * kDim - 1] << std::endl;
   return 0;
}

void dummy_record(int result) {
   // Aggregate the incoming results on the master's side
}
void dummy_update_history(int rank, int* assignement) {
   // Keep track of what tile each worker is working on,
   // and also maybe about times per tile ?
}

// 3
int main3(int argc, char* argv[]) {
   // Printing to the screen


   // magic numbers

   // // number of observations (rows in the input data table)
   // int N = 1000;
   // width of a square tile
   int tile_width = 100;
   // number of contigous column-subsets,
   // with the assumption that all columns
   // in one subset have equal basket-count;
   // by convention, the first subset is strictly
   // composed of the square tiles
   int m = 3;
   // number of square tiles (in one dimension)
   int M = 9;
   // column-counts in the column-subsets
   std::vector<int> n(m);
   n[0] = M * tile_width;
   n[1] = 97;
   n[2] = 3;
   // cumulative-column-counts
   std::vector<int> cum_n(m + 1);
   cum_n[0] = 0;
   for (int i = 0; i < n.size(); i++) {
      cum_n[i + 1] = cum_n[i] + n[i];
   }
   // the number of "divisions"
   int divisions = 4;
   // bucket-counts in column-subsets
   std::vector<int> s(m);
   s[0] = divisions + 1;
   s[1] = divisions + 1;
   s[2] = 2;

   // end of magic numbers
   
   // time-start
   auto start = std::chrono::system_clock::now();

   MPI_Init(&argc, &argv);
   int rank, size;
   const int tag_ON = 1;
   const int tag_OFF = 0;
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   MPI_Comm_size(MPI_COMM_WORLD, &size);
   MPI_Status mpi_status;
   
   // output file
   std::string file_name = "out_" + std::to_string(rank) + ".txt";
   std::ofstream file(file_name);

   const int kRankTreshold = 3;
   assert(kRankTreshold < size);

   int gpu_tile_index[kDim + 1];
   int cpu_tile_index[kDim + 1];

   int result = 0; // identity from the point of view of recording the results
   
   // Scheduler
   if (rank == 0) {

      // make the 2 schedulers
      
      // gpu scheduler
      Indices_triangle_with_semicolon gpu_tile_scheduler = Indices_triangle_with_semicolon(kDim, M, 0, true, kDim);
      gpu_tile_scheduler.use_buff(gpu_tile_index);

      // cpu scheduler
      Indices_triangle_with_semicolon triangle_1 = Indices_triangle_with_semicolon(kDim, m, 1, true, 0);
      std::vector<Indices_triangle> triangle_2;
      std::vector<Indices_triangle> triangle_3;
      std::vector<Indices_product_with_semicolon> product_1;
      std::vector<Indices_sum> sum_1;
      triangle_2.reserve(kDim - 1);
      triangle_3.reserve(kDim - 1);
      product_1.reserve(kDim - 1);
      sum_1.reserve(kDim - 1);
      for (int i = 0; i < kDim - 1; i++) {
         triangle_2.emplace_back(i + 1, M, 0, true);
         triangle_3.emplace_back(kDim - i - 1, m, 1, true);
         product_1.emplace_back(triangle_2[i] ^ triangle_3[i]);
      }
      sum_1.emplace_back(triangle_1 + product_1[0]);
      for (int i = 1; i < kDim - 1; i++)
         sum_1.emplace_back(sum_1[i - 1] + product_1[i]);
      Indices_sum& cpu_tile_scheduler = sum_1.back();
      cpu_tile_scheduler.use_buff(cpu_tile_index);

      //
      int worker_count = size - 1;
      int rank_recv;

      // the main loop
      while(true) {
         MPI_Recv(&result, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &mpi_status);
         dummy_record(result); // is the process of a aggregating the results agnostic about the source?
         rank_recv = mpi_status.MPI_SOURCE;
         std::cout << rank << " received from " << rank_recv << std::endl;
         if (rank_recv < kRankTreshold) {
            // received from a gpu worker
            if (!gpu_tile_scheduler.get_exhausted()) {
               std::cout << rank << " sending " << gpu_tile_scheduler.to_str() << " to " << rank_recv << std::endl;
               // send gpu_index to the gpu worker with tag_ON
               MPI_Send(gpu_tile_index, kDim, MPI_INT, rank_recv, tag_ON, MPI_COMM_WORLD);
               gpu_tile_scheduler.up();
            }
            else {
               std::cout << rank << " sending OFF signal to " << rank_recv << std::endl;
               // terminate the gpu worker by sending tag_OFF
               MPI_Send(gpu_tile_index, kDim, MPI_INT, rank_recv, tag_OFF, MPI_COMM_WORLD);
               --worker_count;
            }
         }
         else {
            // received from a cpu worker
            if (!cpu_tile_scheduler.get_exhausted()) {
               std::cout << rank << " sending " << cpu_tile_scheduler.to_str() << " to " << rank_recv << std::endl;
               // send cpu_index to the cpu worker with tag_ON
               MPI_Send(cpu_tile_index, kDim + 1, MPI_INT, rank_recv, tag_ON, MPI_COMM_WORLD);
               cpu_tile_scheduler.up();
            }
            else if (!gpu_tile_scheduler.get_exhausted()) {
               std::cout << rank << " sending " << gpu_tile_scheduler.to_str() << " to " << rank_recv << std::endl;
               // send gpu_index tp cpu worker with tag_ON
               MPI_Send(gpu_tile_index, kDim + 1, MPI_INT, rank_recv, tag_ON, MPI_COMM_WORLD);
               gpu_tile_scheduler.up();
            }
            else{
               std::cout << rank << " sending OFF signal to " << rank_recv << std::endl;
               // terminate the cpu worker by sending tag_OFF
               MPI_Send(gpu_tile_index, kDim + 1, MPI_INT, rank_recv, tag_OFF, MPI_COMM_WORLD);
               --worker_count;
            }
         }
         if (worker_count == 0) { break; }
      }
      std::cout << rank << " says goodbye" << std::endl;
   }
   
   // Workers
   else if (rank < kRankTreshold) {
      // gpu workers
      MPI_Send(&result, 1, MPI_INT, 0, tag_ON, MPI_COMM_WORLD);
      while (true) {
         std::cout << rank << " waiting" << std::endl;
         MPI_Recv(gpu_tile_index, kDim, MPI_INT,
                  0, MPI_ANY_TAG, MPI_COMM_WORLD, &mpi_status);
         if(mpi_status.MPI_TAG == tag_OFF) {
            std::cout << rank << " says goodbye" << std::endl;
            break;
         }
         std::cout << rank << " received ";
         for (int i = 0; i < kDim; i++) {
            std::cout << gpu_tile_index[i];
         }
         auto end = std::chrono::system_clock::now();
         std::chrono::duration<double> diff = end - start;
         std::cout << " at " << diff.count() << std::endl;
         int gpu_column_ranges[2 * kDim];
         gpu_tile_index2column_ranges(gpu_column_ranges,
                                      gpu_tile_index,
                                      tile_width);
         result = dummy_work(gpu_column_ranges, rank); // print
         MPI_Send(&result, 1, MPI_INT,
                  0, tag_ON, MPI_COMM_WORLD);
      }
   }
   else {
      // cpu workers
      MPI_Send(&result, 1, MPI_INT, 0, tag_ON, MPI_COMM_WORLD);
      while (true) {
         std::cout << rank << " waiting" << std::endl;
         MPI_Recv(cpu_tile_index, kDim + 1, MPI_INT,
                  0, MPI_ANY_TAG, MPI_COMM_WORLD, &mpi_status);

         if(mpi_status.MPI_TAG == tag_OFF) {
            std::cout << rank << " says goodbye" << std::endl;
            break;
         }
         std::cout << rank << " received ";
         for (int i = 0; i < kDim + 1; i++) {
            std::cout << cpu_tile_index[i];
         }
         std::cout << std::endl;
         int cpu_column_ranges[2 * kDim];
         cpu_tile_index2column_ranges(cpu_column_ranges,
                                      cpu_tile_index,
                                      tile_width,
                                      cum_n);
         result = dummy_work(cpu_column_ranges, rank); // print
         
         // prepare k-tuple scheduler
         int semicolon = cpu_tile_index[kDim];
         std::vector<Index_count> counts_square
                   = index_counts(cpu_tile_index, semicolon);
         std::vector<Index_count> counts_rectangle
                   = index_counts(cpu_tile_index + semicolon, kDim - semicolon);
         std::vector<Indices_triangle> triangle_4;
         triangle_4.reserve(counts_square.size() + counts_rectangle.size());
         for (auto& ic : counts_square) {
            triangle_4.emplace_back(ic.count, tile_width, 0, false);
         }
         for (auto& ic : counts_rectangle) {
            triangle_4.emplace_back(ic.count, n[ic.index], 0, false);
         }
         std::vector<Indices_product> product_2;
         product_2.reserve(triangle_4.size());
         Indices_product_id id; 
         product_2.emplace_back(id * triangle_4[0]);
         for (int i = 1; i < triangle_4.size(); i++) {
            product_2.emplace_back(product_2.back() * triangle_4[i]);
         }
         Indices_product &cpu_tuple_scheduler = product_2.back();
         int cpu_tuple_index[kDim];
         cpu_tuple_scheduler.use_buff(cpu_tuple_index);
         // print out first few column k-tuples
         int few = 5;
         int print_count = 0;
         std::cout << rank << " column k-tuples: ";
         int cpu_tuple_columns[kDim];
         while (!cpu_tuple_scheduler.get_exhausted() & print_count < few) { 
            cpu_tuple_index2columns(cpu_tuple_columns,
                                    cpu_tile_index,
                                    cpu_tuple_index,
                                    tile_width,
                                    cum_n);
            for (int i = 0; i < kDim - 1; i++) {
               std::cout << cpu_tuple_columns[i] << ",";
            }
            std::cout << cpu_tuple_columns[kDim-1] <<";";
            cpu_tuple_scheduler.up();
            ++print_count;
         }
         std::cout << "..." << std::endl;
         // end of tuple business
         
         MPI_Send(&result, 1, MPI_INT,
                  0, tag_ON, MPI_COMM_WORLD);
      }
   }
   
   file.close();
   MPI_Finalize();
      
   return 0;
}


int dummy_work_2(int *column_ranges,
                 std::chrono::duration<double> time,
                 int rank,
                 std::ofstream& file) {
   std::this_thread::sleep_for(std::chrono::milliseconds(random() % 100));
   file << time.count() << ", ";
   for (int i = 0; i < 2 * kDim - 1; i++) {
      file << column_ranges[i] << ", ";
   }
   file << column_ranges[2 * kDim - 1] << std::endl;
   return 0;
}


// 4 
int main4(int argc, char* argv[]) {
   //Printing to files


   // magic numbers

   // // number of observations (rows in the input data table)
   // int N = 1000;
   // width of a square tile
   int tile_width = 100;
   // number of contigous column-subsets,
   // with the assumption that all columns
   // in one subset have equal basket-count;
   // by convention, the first subset is strictly
   // composed of the square tiles
   int m = 3;
   // number of square tiles (in one dimension)
   int M = 9;
   // column-counts in the column-subsets
   std::vector<int> n(m);
   n[0] = M * tile_width;
   n[1] = 97;
   n[2] = 3;
   // cumulative-column-counts
   std::vector<int> cum_n(m + 1);
   cum_n[0] = 0;
   for (int i = 0; i < n.size(); i++) {
      cum_n[i + 1] = cum_n[i] + n[i];
   }
   // the number of "divisions"
   int divisions = 4;
   // bucket-counts in column-subsets
   std::vector<int> s(m);
   s[0] = divisions + 1;
   s[1] = divisions + 1;
   s[2] = 2;

   // end of magic numbers
   
   // time-start
   auto start = std::chrono::system_clock::now();

   MPI_Init(&argc, &argv);
   int rank, size;
   const int tag_ON = 1;
   const int tag_OFF = 0;
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   MPI_Comm_size(MPI_COMM_WORLD, &size);
   MPI_Status mpi_status;
   
   // how many gpu workers?
   const int kRankTreshold = 3;
   assert(kRankTreshold < size);  

   // random seed
   srand(rank);
   
   // output file
   std::string file_name = "out_" + std::to_string(rank) + ".csv";
   std::ofstream file(file_name);

   int gpu_tile_index[kDim + 1];
   int cpu_tile_index[kDim + 1];

   int result = 0; // identity from the point of view of recording the results
   
   // Scheduler
   if (rank == 0) {

      // make the 2 schedulers
      
      // gpu scheduler
      Indices_triangle_with_semicolon gpu_tile_scheduler = Indices_triangle_with_semicolon(kDim, M, 0, true, kDim);
      gpu_tile_scheduler.use_buff(gpu_tile_index);

      // cpu scheduler
      Indices_triangle_with_semicolon triangle_1 = Indices_triangle_with_semicolon(kDim, m, 1, true, 0);
      std::vector<Indices_triangle> triangle_2;
      std::vector<Indices_triangle> triangle_3;
      std::vector<Indices_product_with_semicolon> product_1;
      std::vector<Indices_sum> sum_1;
      triangle_2.reserve(kDim - 1);
      triangle_3.reserve(kDim - 1);
      product_1.reserve(kDim - 1);
      sum_1.reserve(kDim - 1);
      for (int i = 0; i < kDim - 1; i++) {
         triangle_2.emplace_back(i + 1, M, 0, true);
         triangle_3.emplace_back(kDim - i - 1, m, 1, true);
         product_1.emplace_back(triangle_2[i] ^ triangle_3[i]);
      }
      sum_1.emplace_back(triangle_1 + product_1[0]);
      for (int i = 1; i < kDim - 1; i++)
         sum_1.emplace_back(sum_1[i - 1] + product_1[i]);
      Indices_sum& cpu_tile_scheduler = sum_1.back();
      cpu_tile_scheduler.use_buff(cpu_tile_index);

      //
      int worker_count = size - 1;
      int rank_recv;

      // the main loop
      while(true) {
         MPI_Recv(&result, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &mpi_status);
         rank_recv = mpi_status.MPI_SOURCE;
         auto end = std::chrono::system_clock::now();
         std::chrono::duration<double> diff = end - start;
         file << diff.count() << ", " << rank_recv << std::endl;
         if (rank_recv < kRankTreshold) {
            // received from a gpu worker
            if (!gpu_tile_scheduler.get_exhausted()) {
               // send gpu_index to the gpu worker with tag_ON
               MPI_Send(gpu_tile_index, kDim, MPI_INT, rank_recv, tag_ON, MPI_COMM_WORLD);
               gpu_tile_scheduler.up();
            }
            else {
               // terminate the gpu worker by sending tag_OFF
               MPI_Send(gpu_tile_index, kDim, MPI_INT, rank_recv, tag_OFF, MPI_COMM_WORLD);
               --worker_count;
            }
         }
         else {
            // received from a cpu worker
            if (!cpu_tile_scheduler.get_exhausted()) {
               // send cpu_index to the cpu worker with tag_ON
               MPI_Send(cpu_tile_index, kDim + 1, MPI_INT, rank_recv, tag_ON, MPI_COMM_WORLD);
               cpu_tile_scheduler.up();
            }
            else if (!gpu_tile_scheduler.get_exhausted()) {
               // send gpu_index tp cpu worker with tag_ON
               MPI_Send(gpu_tile_index, kDim + 1, MPI_INT, rank_recv, tag_ON, MPI_COMM_WORLD);
               gpu_tile_scheduler.up();
            }
            else{
               // terminate the cpu worker by sending tag_OFF
               MPI_Send(gpu_tile_index, kDim + 1, MPI_INT, rank_recv, tag_OFF, MPI_COMM_WORLD);
               --worker_count;
            }
         }
         dummy_record(result); // is the process of a aggregating the results agnostic about the source?
         if (worker_count == 0) { break; }
      }
   }
   
   // Workers
   else if (rank < kRankTreshold) {
      // gpu workers
      MPI_Send(&result, 1, MPI_INT, 0, tag_ON, MPI_COMM_WORLD);
      while (true) {
         MPI_Recv(gpu_tile_index, kDim, MPI_INT,
                  0, MPI_ANY_TAG, MPI_COMM_WORLD, &mpi_status);
         if(mpi_status.MPI_TAG == tag_OFF) {
            break;
         }
         auto end = std::chrono::system_clock::now();
         std::chrono::duration<double> time = end - start;
         int gpu_column_ranges[2 * kDim];
         gpu_tile_index2column_ranges(gpu_column_ranges,
                                      gpu_tile_index,
                                      tile_width);
         result = dummy_work_2(gpu_column_ranges, time, rank, file); // print to file
         MPI_Send(&result, 1, MPI_INT,
                  0, tag_ON, MPI_COMM_WORLD);
      }
   }
   else {
      // cpu workers
      MPI_Send(&result, 1, MPI_INT, 0, tag_ON, MPI_COMM_WORLD);
      while (true) {
         MPI_Recv(cpu_tile_index, kDim + 1, MPI_INT,
                  0, MPI_ANY_TAG, MPI_COMM_WORLD, &mpi_status);

         if(mpi_status.MPI_TAG == tag_OFF) {
            break;
         }
         auto end = std::chrono::system_clock::now();
         std::chrono::duration<double> time = end - start;
         int cpu_column_ranges[2 * kDim];
         cpu_tile_index2column_ranges(cpu_column_ranges,
                                      cpu_tile_index,
                                      tile_width,
                                      cum_n);
         result = dummy_work_2(cpu_column_ranges, time, rank, file); // print to file
         MPI_Send(&result, 1, MPI_INT,
                  0, tag_ON, MPI_COMM_WORLD);
      }
   }
   
   file.close();
   MPI_Finalize();
      
   return 0;
}

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


// 5
int main_5() {

   // Assumption about placement of the regions of columns along the dataset:

   // |          signif         |        nonsignif        |         contrast        |
   // |     fat    |    lean    |     fat    |    lean    |     fat    |    lean    |
   // |homog. tiles|heter. tiles|homog. tiles|heter. tiles|homog. tiles|heter. tiles|
   // R[0]         R[1]         R[2]         R[3]         R[4]         R[5]         R[6]
   // T[R[0]] ...  T[R[1]] ...  T[R[2]] ...  T[R[3]] ...  T[R[4]] ...  T[R[5]] ...  T[R[6]]

   // and the homogeanous tiles are (either mostly- or only-) square-shaped.
   // The array T(iles) contains column indeces that are the borders of tiles.
   // In other words, where one tile ends and other begins is given by
   // consequtive values in T.
   // The array R(egions) contains chosen indices of T with the meaning explained above.
   // The array D, not mentioned above, of size T.size - 1, will contain dof's of each tile.
   
  // int T[] = {0,   // 0 // fat // signif
  //            100, // 1 
  //            200, // 2
  //            222, // 3 // lean
  //            235, // 4
  //            247, // 5
  //            258, // 6 // fat // nonsignif
  //            358, // 7
  //            458, // 8
  //            558, // 9
  //            658, // 10
  //            758, // 11
  //            800, // 12 // lean
  //            842, // 13
  //            876, // 14
  //            876, // 15 // fat // contrast
  //            976, // 16
  //           1076, // 17
  //           1200, // 18 // lean
  //           1230, // 19
  //           1230, // 20
  //           1245};// 21 // end <= the total column-count

   int R[7] = {0, 3, 6, 12, 15, 18, 21};

   assert(R[2] - R[1] == R[4] - R[3] | R[4] - R[3] == 0);
   assert(R[2] - R[1] == R[6] - R[5] | R[6] - R[5] == 0);

  // // number of dicrete buckets per tile
  // int D[] = {5, // 0 // fat // signif
  //            5, // 1 
  //            5, // 2
  //            6, // 3 // lean
  //            7, // 4
  //            8, // 5
  //            5, // 6 // fat // nonsignif
  //            5, // 7
  //            5, // 8
  //            5, // 9
  //            5, // 10
  //            5, // 11
  //            6, // 12 // lean
  //            7, // 13
  //            8, // 14
  //            5, // 15 // fat // contrast
  //            5, // 14
  //            5, // 15
  //            6, // 16 // lean
  //            7, // 17
  //            8};// 18

   // Construct the series of homogenous k-dimensional tiles
   Indices_triangle signif_fat_1(kDim, R[1], R[0], true);
   Indices_triangle signif_fat_2(kDim - 1, R[1], R[0], true);
   Indices_triangle nonsignif_fat_2(1, R[3], R[2], true);
   Indices_triangle contrast_fat_2(1, R[5], R[4], true);
   Indices_sum nonsignif_plus_contrast_fat_2(nonsignif_fat_2, contrast_fat_2);
   Indices_product fat_2(signif_fat_2, nonsignif_plus_contrast_fat_2);
   Indices_sum homogen(signif_fat_1, fat_2);

   int buff[kDim];
   homogen.use_buff(buff);
   homogen.print();   
   
   std::cout << " ------ " << std::endl;

   // Construct the remaining series of k-dimensional tiles
   //
   std::vector<Indices_triangle> signif_fat_3;
   std::vector<Indices_triangle> signif_lean_3;
   std::vector<Indices_product> signif_fat_x_lean_3;
   std::vector<Indices_sum> signif_fat_x_lean_3_sum;
   Indices_sum_id zero_3(kDim);
   Indices_sum& heterogen_signif_3 = fold(signif_fat_3,
                                         signif_lean_3,
                                         signif_fat_x_lean_3,
                                         signif_fat_x_lean_3_sum,
                                         zero_3,
                                         R[0], R[1],
                                         R[1], R[2],
                                         1,
                                         kDim);
   //
   std::vector<Indices_triangle> signif_fat_4;
   std::vector<Indices_triangle> signif_lean_4;
   std::vector<Indices_product> signif_fat_x_lean_4;
   std::vector<Indices_sum> signif_fat_x_lean_4_sum;
   Indices_sum_id zero_4(kDim - 1);
   Indices_sum& heterogen_signif_4 = fold(signif_fat_4,
                                          signif_lean_4,
                                          signif_fat_x_lean_4,
                                          signif_fat_x_lean_4_sum,
                                          zero_4,
                                          R[0], R[1],
                                          R[1], R[2],
                                          1,
                                          kDim - 1);
   Indices_triangle nonsignif_and_contrast_4(1, R[6], R[2], true);
   Indices_product heterogen_mixed_4(heterogen_signif_4, nonsignif_and_contrast_4);
   //
   Indices_triangle signif_fat_5(kDim - 1, R[1], R[0], true);
   Indices_triangle nonsignif_lean_5(1, R[4], R[3], true);
   Indices_triangle contrast_lean_5(1, R[6], R[5], true);
   Indices_sum nonsignif_plus_contrast_lean_5(nonsignif_lean_5, contrast_lean_5);
   Indices_product heterogen_mixed_5(signif_fat_5, nonsignif_plus_contrast_lean_5);
   //
   Indices_sum heterogen_sum_1(heterogen_signif_3, heterogen_mixed_4);
   Indices_sum heterogen(heterogen_sum_1, heterogen_mixed_5);


   heterogen.use_buff(buff);
   heterogen.print();



   return 0;
}



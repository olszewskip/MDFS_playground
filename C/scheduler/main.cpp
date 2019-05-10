#include <mpi.h>
// uncomment to disable assert()
// #define NDEBUG
#include <cassert>
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <vector>

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
         //std::cout << "|shape: " << dim_0 << ", " << dim_1 << "\n";
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
      Indices_triangle (int k_arg, int n_arg, int offset_arg=0, bool with_diag_arg=true) :
         n(n_arg),
         offset(offset_arg),
         with_diag(with_diag_arg) {
         assert(k_arg > 0); 
         assert(offset >= 0);
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
      void set_exhausted(bool exhausted_var) { exhausted = exhausted_var; }
      void up() {
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
      void reset() {
         for(int i = 0; i < k; i++)
            index[i] = offset + !with_diag * i;
         exhausted = false;
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

class Indices_triangle_with_semicolon : public Indices_triangle {
   // It uses a buffer that is assumed to be one integer
   // longer than for the Indices_triangle class. Value
   // stored in this additional memory is 0 (which sounds
   // like it's not very useful, I know).
   public:
      Indices_triangle_with_semicolon(int k_arg, int n_arg,
                                      int offset_arg=0,
                                      bool with_diag_arg=true) :
         Indices_triangle (k_arg, n_arg, offset_arg, with_diag_arg) {}
      void use_buff(int* const buff) {
         Indices_triangle::use_buff(buff);
         // Append 0.
         index[k] = 0;
      }
      std::string to_str() {
         std::stringstream ss;
         ss << Indices_triangle::to_str() << " " << index[k];
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

class Indices_product_left_id : public Indices {
   // Multiplication's left identity.
   protected:
      bool exhausted = false;
   public:
      Indices_product_left_id() {k = 0;}
      bool get_exhausted() { return exhausted; }
      void set_exhausted(bool exhausted_arg) { exhausted = exhausted_arg; }
      void use_buff(int* const buff) {}
      void up() { exhausted = true; }
      void reset() { exhausted = false; }
      std::string to_str() { return "."; }
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
      bool get_exhausted() {
         return ind_R.get_exhausted();
      }
      void set_exhausted(bool exhausted_var) {
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
   Indices_triangle my_indices_1(kDim - 3, 7, 2);
   Indices_triangle my_indices_2(3, 4);
   Indices_triangle my_indices_3(kDim - 1, 6, 3); 
   Indices_triangle my_indices_4(1, 10, 4);
   Indices_product_left_id id;

   std::cout << id.to_str() << std::endl;

   Indices_product product_01 = id * my_indices_1;
   Indices_product product_12 = product_01 * my_indices_2;
   Indices_product_with_semicolon product_34 = my_indices_3 ^ my_indices_4;
   Indices_sum sum_1234 = product_12 + product_34;
   
   int indices_buff_1234[kDim + 1];
   sum_1234.use_buff(indices_buff_1234);
   sum_1234.print();
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

void gpu_index2columns(int gpu_columns[],
                       int gpu_index[], int tile_width) {
   // Bijection between tile-index and a 2d matrix (row-major).
   // that represent the tile by specifying column-indices.
   // Each two-element row in the k * 2 matrix contains beginning
   // and end (exclusive) of a column-indices-range.
   gpu_columns[2 * kDim];
   for (int i = 0; i < kDim; i++) {
      gpu_columns[2 * i] = tile_width * gpu_index[i];
      gpu_columns[2 * i + 1] = tile_width * (gpu_index[i] + 1);
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


void cpu_index2columns(int cpu_columns[],
                       int tile_index[], int semicolon, int tile_width,
                       int tuple_index[], std::vector<int> cum_n) {
   // Fill vector with column indices
   // that represents a single k-tuple.
   cpu_columns[kDim];
   for (int i = 0; i < semicolon; i++)
      cpu_columns[i] = tile_width * tile_index[i] + tuple_index[i];
   for (int i = semicolon; i < kDim; i++)
      cpu_columns[i] = cum_n[tile_index[i]] + tuple_index[i];
}


// 2
int main() {

   // the dimensionality of k-tuples
   const int kDim = 3;
   // // number of observations (rows in the input data table)
   // int N = 1000;
   // width of a square tile
   int tile_width = 100;
   // number of contigous column-subsets,
   // whith the assumption that all columns
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
   
   Indices_triangle gpu_scheduler = Indices_triangle(kDim, M);
   int gpu_buff[kDim];
   gpu_scheduler.use_buff(gpu_buff);
   std::cout << "gpu" << std::endl;
   std::cout << " ----- \n";
   
   while (!gpu_scheduler.get_exhausted()) {
      // The gpu_scheduler's index buffer is (presumably)
      // sent over with mpi.
      //std::cout << gpu_scheduler.to_str() << std::endl;
      
      // For each such buffer we do the following.
      int* tile_index = gpu_scheduler.index;
      
      int gpu_columns[2 * kDim];
      gpu_index2columns(gpu_columns,
                        tile_index, tile_width);
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
   Indices_triangle_with_semicolon triangle_1 = Indices_triangle_with_semicolon(kDim, m, 1);
   std::vector<Indices_triangle> triangle_2;
   std::vector<Indices_triangle> triangle_3;
   std::vector<Indices_product_with_semicolon> product_1;
   std::vector<Indices_sum> sum_1;
   triangle_2.reserve(kDim - 1);
   triangle_3.reserve(kDim - 1);
   product_1.reserve(kDim - 1);
   sum_1.reserve(kDim - 1);

   for (int i = 0; i < kDim - 1; i++) {
      triangle_2.emplace_back(i + 1, M);
      triangle_3.emplace_back(kDim - i - 1, m, 1);
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
      int* tile_index = cpu_scheduler.index;
      
      // Aggregate the numbers in tile_index.
      // The boundary between square and rectangle
      // indices is denoted by semicolon.
      int semicolon = tile_index[kDim];
      std::vector<Index_count> counts_square
                = index_counts(tile_index, semicolon);
      std::vector<Index_count> counts_rectangle
                = index_counts(tile_index + semicolon, kDim - semicolon);
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
      Indices_product_left_id id; 
      product_2.emplace_back(id * triangle_4[0]);
      for (int i = 1; i < triangle_4.size(); i++) {
         product_2.emplace_back(product_2.back() * triangle_4[i]);
      }
      Indices_product &tuple_scheduler = product_2.back();
      int tuple_index[kDim];
      tuple_scheduler.use_buff(tuple_index);
     
      // Print a first few k-tuples of columns.
      int few = 3; 
      for (int i = 0; i < few; i++)
      //while (!tuple_scheduler.get_exhausted())
      {
         int cpu_columns[kDim];
         cpu_index2columns(cpu_columns,
                           tile_index, semicolon, tile_width,
                           tuple_index, cum_n);
         for (int i = 0; i < kDim; i++)
            std::cout << cpu_columns[i] << " ";
         std::cout << std::endl;

         tuple_scheduler.up();
      }

      std::cout << " ----- \n";
      cpu_scheduler.up();
   }
}




int gpu_dummy_work(Matrix2D<int>& columns) {
   int sum = 0;
   for (int i = 0; i < columns.dim_0; i++)
      for (int j = 0; j < columns.dim_1; j++)
         sum += columns(i, j);
   return sum % 123;
}

void dummy_record(int result) {
   // Aggregate the incoming results on the master's side
}
void dummy_update_history(int rank, int* assignmenet) {
   // Keep track of what tile each worker is working on,
   // and also maybe about times per tile ?
}

// 3
int main3(int argc, char* argv[]) {
   
   // the dimensionality of k-tuples
   const int kDim = 2;
   // // number of observations (rows in the input data table)
   // int N = 1000;
   // width of a square tile
   int tile_width = 100;
   // number of contigous column-subsets,
   // whith the assumption that all columns
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
   // the number of "divisions"
   int divisions = 4;
   // bucket-counts in column-subsets
   std::vector<int> s(m);
   s[0] = divisions + 1;
   s[1] = divisions + 1;
   s[2] = 2;
   

   MPI_Init(&argc, &argv);
   int rank, size;
   const int tag_ON = 1;
   const int tag_OFF = 0;
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   MPI_Comm_size(MPI_COMM_WORLD, &size);

   int gpu_index[kDim];
   int result = 0; // identity from the point of view of recording the results
   MPI_Status mpi_status;
   int source;
   
   if (rank == 0) {
      // master
      Indices_triangle gpu_scheduler = Indices_triangle(kDim, M);
      gpu_scheduler.use_buff(gpu_index);
      while(!gpu_scheduler.get_exhausted()) {
         MPI_Recv(&result, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &mpi_status);
         dummy_record(result); // is the process of a aggregating the results agnostic about the source?
         source = mpi_status.MPI_SOURCE;
         if (source < 100)
            // send to gpu workers
            MPI_Send(gpu_index, kDim, MPI_INT, mpi_status.MPI_SOURCE, tag_ON, MPI_COMM_WORLD);
            gpu_scheduler.up();
         //else
         //   // send to cpu workers
         //   TODO
         //dummy_record_result(result);
         //dummy_update_history(mpi_status.MPI_SOURCE, gpu_index);
      }
      for (int i = 1; i < size; i++) {
         // Receive the last result per worker
         // and send terminating signal via the tag
         MPI_Recv(&result, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &mpi_status);
         dummy_record(result);
         source = mpi_status.MPI_SOURCE;
         MPI_Send(gpu_index, kDim, MPI_INT, source, tag_OFF, MPI_COMM_WORLD);
      }
   }
   
   else if (rank < 100) {
      // gpu workers
      MPI_Send(&result, 1, MPI_INT, 0, tag_ON, MPI_COMM_WORLD);
      while(true) {
         MPI_Recv(gpu_index, kDim, MPI_INT, 0, MPI_ANY_TAG, MPI_COMM_WORLD, &mpi_status);
         if(!mpi_status.MPI_TAG) { break; }
         Matrix2D<int> gpu_columns = gpu_index2pretty_columns(gpu_index, tile_width);
         gpu_columns.print();
         result = gpu_dummy_work(gpu_columns);
         MPI_Send(&result, 1, MPI_INT, 0, tag_ON, MPI_COMM_WORLD);
      }
   }

   // // cpu workers
   // else
   //    TODO

   MPI_Finalize();
      
   return 0;
}


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
      // The flag "exhausted" is not yet set when the
      // indices represent the last element of the sequence.
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
      void set_exhausted(bool exhausted_arg) {
         exhausted = exhausted_arg;
      }
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

int dummy_work(int *column_ranges, int rank) {
   // Print to screen
   std::this_thread::sleep_for(
         std::chrono::milliseconds(random() % 100));
   std::cout << rank << " column-ranges: ";
   for (int i = 0; i < kDim - 1; i++) {
      std::cout << column_ranges[2 * i] << ":"
                << column_ranges[2 * i + 1] <<", ";
   }
   std::cout << column_ranges[2 * kDim - 2] << ":"
             << column_ranges[2 * kDim - 1] << std::endl;
   return 0;
}

int dummy_work_2(int *column_ranges,
                 std::chrono::duration<double> time,
                 int rank,
                 std::ofstream& file) {
   // Print to file
   std::this_thread::sleep_for(
         std::chrono::milliseconds(random() % 100));
   file << time.count() << ", ";
   for (int i = 0; i < 2 * kDim - 1; i++) {
      file << column_ranges[i] << ", ";
   }
   file << column_ranges[2 * kDim - 1] << std::endl;
   return 0;
}

void dummy_record(int result) {
   // Aggregate the incoming results on the master's side
}
void dummy_update_history(int rank, int* assignement) {
   // Keep track of what tile each worker is working on,
   // and also maybe about times per tile ?
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

int main(int argc, char* argv[]) {

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

   int gpu_tile_index[kDim];
   int cpu_tile_index[kDim];
   
   int result = 0; // identity from the point of view of recording the results

   // Scheduler
   if (rank = 0) {

      // make the 2 schedulers

      // gpu scheduler
      //
      Indices_triangle signif_fat_1(kDim, R[1], R[0], true);
      //
      Indices_triangle signif_fat_2(kDim - 1, R[1], R[0], true);
      Indices_triangle nonsignif_fat_2(1, R[3], R[2], true);
      Indices_triangle contrast_fat_2(1, R[5], R[4], true);
      Indices_sum nonsignif_plus_contrast_fat_2(nonsignif_fat_2, contrast_fat_2);
      Indices_product fat_2(signif_fat_2, nonsignif_plus_contrast_fat_2);
      //
      Indices_sum gpu_tile_scheduler(signif_fat_1, fat_2);
      //
      gpu_tile_scheduler.use_buff(gpu_tile_index);

      // cpu scheduler
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
      Indices_sum cpu_tile_scheduler(heterogen_sum_1, heterogen_mixed_5);
      //
      cpu_tile_scheduler.use_buff(cpu_tile_index);

      // count your workers
      int worker_count = size - 1;
      int rank_recv;

      // the main loop
      while (true) {
         MPI_Recv(
    
      




   

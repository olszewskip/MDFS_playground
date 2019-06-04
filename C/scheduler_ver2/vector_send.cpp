//#include <mpi.h>
// uncomment to disable assert()
// #define NDEBUG
#include <cassert>
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <vector>
#include <chrono>
#include <thread>

const int kDim = 3;
const int kCuriosity = 2;

struct tuple_result {
   int columns[kDim];
   int IGs[kDim];
};

struct column_IG_context {
   int column;
   double IG;
   int context[kDim - 1];
   std::string to_str() {
      std::stringstream ss;
      ss << "column: " << column << ", "
         << "IG: " << IG << ", "
         << "context: ";
      for (int i = 0; i < kDim - 2; i++)
         ss << context[i] << ", ";
      ss << context[kDim - 2];
      return ss.str();
    }
};

void copy_with_omission(int* longer,
                        int* shorter,
                        int index_to_omit) {
   // longer = {1, 2, 3}; index_to_omit = 1
   // => shorter = {1, 3}
   assert(index_to_omit < kDim);
   bool thrown_away = false;
   for (int i = 0; i < kDim; i++) {
      if (i == index_to_omit) {
         thrown_away = true;
         continue;
      }
      else {
         shorter[i - thrown_away] = longer[i];
      }
   }
}

void fill_column_IG_context_array(column_IG_context array[],
                                  tuple_result result) {
   // result = {{1, 2, 3}, {123, 234, 345}} =>
   // array = {{1, 123, {2,3}}, {2, 234, {1,3}}, {3, 345, {1,2}}}
   for (int i = 0; i < kDim; i++) {
      array[i].column = result.columns[i];
      array[i].IG = result.IGs[i];
      copy_with_omission(result.columns, array[i].context, i);
   }
} 

int argmin_of_IGs(double IGs[]) {
   int i = 0;
   double min_IG = IGs[0];
   for (int j = 1; j < kCuriosity; j++) {
       if (IGs[j] < min_IG) { 
	   min_IG = IGs[j];
	   i = j;
       }
   }
   return i;
}

struct tile_result_entry {
   int column;
   double IGs[kCuriosity];
   int contexts[kCuriosity][kDim - 1];
   tile_result_entry(column_IG_context arg):
      // initialize with a first pair
      // of the IG and context[] arguments
      column(arg.column) {
      IGs[0] = arg.IG;
      for (int j = 0; j < kDim - 1; j++)
         contexts[0][j] = arg.context[j];
      for (int i = 1; i < kCuriosity; i++) {
         IGs[i] = 0;
         for (int j = 0; j < kDim - 1; j++)
             contexts[i][j] = -1;
      }
   }
   void update (double new_IG, int new_context[]) {
       int argmin = argmin_of_IGs(IGs);
       if (IGs[argmin] < new_IG) {
           IGs[argmin] = new_IG;
           for (int i = 0; i < kDim - 1; i++)
               contexts[argmin][i] = new_context[i];
       }
   }
   std::string to_str() {
      std::stringstream ss;
      ss << "column: " << column << " | ";
      for (int i = 0; i < kCuriosity; i++) {
	 if (IGs[i] == 0) { break; }
         ss << "IG: " << IGs[i] << ", context: ";
         for (int j = 0; j < kDim - 2; j++)
             ss << contexts[i][j] << ", ";
         ss << contexts[i][kDim - 2] << " | ";
      }
      return ss.str();
   }

};

struct tile_result {
   std::vector<tile_result_entry> vect;
   void record(column_IG_context arg) {
      bool not_found = true;
      for (auto& entry : vect) {
          if (entry.column == arg.column) {
              not_found = false;
              entry.update(arg.IG, arg.context);
              break;
          }
      }
      if (not_found) {
          vect.emplace_back(arg);
      }
   }
   void record(tuple_result result) {
      column_IG_context array[kDim];
      fill_column_IG_context_array(array, result);
      for (int i = 0; i < kDim; i++)
          record(array[i]);
   }
   void print() {
      for (auto& entry : vect)
          std::cout << entry.to_str() << std::endl;
   }
};



int main() {
  
   //tuple_result some_tuple_result {{1, 2, 3}, {101, 102, 103}};
   //tuple_result some_other_tuple_result {{1, 3, 5}, {104, 105, 106}}; 
  // column_IG_context array[kDim];
  // fill_column_IG_context_array(array, some_tuple_result); 
  // for (int i = 0; i < kDim; i++)
  //    std::cout << array[i].to_str() << std::endl;
  // 
  // tile_result_entry entry(array[0]);
  // std::cout << entry.to_str() << std::endl;
   tile_result some_tile_result;
  // some_tile_result.record(array[0]);
   //some_tile_result.record(some_tuple_result);
   //some_tile_result.record(some_other_tuple_result);
   for (int i = 0; i < 1000; i++) {
      some_tile_result.record({{rand()%10, rand()%10, rand()%10},
                               {rand()%100, rand()%100, rand()%100}});
   }
   some_tile_result.print();
   
   

   //column_IG_context array[kDim];
   //fill_column_IG_context_array(array, result);
   //for (int i = 0; i < kDim; i++) {
   //	std::cout << array[i].to_str() << std::endl;
   //}

   //std::cout << sizeof(tile_result_entry) << std::endl;



  // MPI_Init();
  // int rank, size;
  // MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  // MPI_Comm_size(MPI_COMM_WORLD, &size);
  // MPI_Status mpi_status;
  // const int TAG = 0;

  // if (rank == 1) {
  //    // make a random vector that represents
  //    // results from a single tile
  //    std::vector<tile_result_entry> tile_result;
  //    // let's say there are 10 entries in the vector
  //    for (int i = 0; i < 10; i++) {
  //       // each consists of a key-value pair,
  //       // the key being a single column-index,
  //       // and the value being a list of length "curiosity".
  //       // The list contains pairs: column-context + information-gain.
  //       int column = (int)random()%10;
  //       int context_IG[curiosity];
  //       for (int = 0; i < curiosity; i++) {
  //          int context[kDim - 1];
  //          for (int i = 0; i < kDim - 1; i++) {
  //             context[i] = (int)random%10;
  //          }
  //          int IG = (int)random%10;
  //          context_IG[i] = {context, IG};
  //       }
  //       tile_result.push_back(column, context_IG);
  //    }

  //    std::cout << rank << std::endl;
  //    for (auto& some_result : vect) {
  //       std::cout << some_result.column << ": ";
  //       for (int i = 0; i < curiosity; i++)
  //          std::cout << some_result.info_gains[i];
  //       std::cout << std::endl;
  //    }
  //    std::cout << " ---- " << std::endl;
  //    MPI_Send (vect.data(), sizeof(*vect.data()) * vect.size(), MPI_INT,
  //              1, TAG, MPI_COMM_WORLD);
  // }

  // else if (rank == 1) {
  //    std::this_thread::sleep_for(std::chrono::milliseconds(10));

  //    const int static_buffer_size = 500;
  //    int static_buffer[static_buffer_size];
  //    int* dynamic_buffer;
  //    bool large = false; // does the dynamic buffer need cleaning?
  //    dict_entry* recieved_dict_entries;
  //    int incoming_count;

  //    MPI_Get_count (mpi_status, MPI_INT, &incoming_count);
  //    if (incoming_count > static_buffer) {
  //       large = true;
  //       dynamic_buffer = new int[incoming_count];
  //       MPI_Recv (dynamic_buffer, incoming_count, MPI_INT,
  //                 0, TAG, MPI_COMM_WORLD, mpi_status);
  //       recieved_dict_entries;
  //    }
  //    else {
  //       MPI_Recv (static_buffer, incoming_count, MPI_INT,
  //                 0, TAG, MPI_COMM_WORLD, mpi_status);
  //    }

  //    std::vector<dict_entry> vect 





   return 0;
}



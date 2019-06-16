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
#include <random>
#include <stddef.h>

// the dimensionality of k-tuples
const int kDim = 3;
// number of largest IG per column to remember
const int kCuriosity = 2;

#include "Indices.h"
#include "data_structures.h"

int main(int argc, char* argv[]) {

   // MPI initialization
   MPI_Init(&argc, &argv);
   int rank, size;
   const int TAG_ON = 1;
   const int TAG_OFF = 0;
   // MPI_Comm COMM = MPI_COMM_WORLD;
   MPI_Comm_rank(MPI_COMM_WORLD, &rank);
   MPI_Comm_size(MPI_COMM_WORLD, &size);
   MPI_Status mpi_status;

   // create an mpi-type for sending
   // arrays of tile_result_entries
   int type_count = 2;
   int blocklengths[] = {1 + kCuriosity * (kDim - 1), kCuriosity};
   MPI_Aint displacements[] = {offsetof(tile_result_entry, column),
                               offsetof(tile_result_entry, IGs)};
   MPI_Datatype types[] = {MPI_INT, MPI_DOUBLE};
   MPI_Datatype tmp_type, mpi_tile_result_entry;
   MPI_Aint lower_bound, extent;
   MPI_Type_create_struct(type_count, blocklengths, displacements,
                          types, &tmp_type);
   MPI_Type_get_extent(tmp_type, &lower_bound, &extent);
   MPI_Type_create_resized(tmp_type, lower_bound, extent,
                           &mpi_tile_result_entry);
   MPI_Type_commit(&mpi_tile_result_entry);
   //
   
   tile_result the_tile_result;
   std::cout << rank << " " << the_tile_result.vect.size() << std::endl;

   if (rank == 0) {
      tile_result_entry* recv_buffer = new tile_result_entry[100];
      int incoming_count;
      MPI_Probe(MPI_ANY_SOURCE, TAG_ON, MPI_COMM_WORLD, &mpi_status);
      MPI_Get_count(&mpi_status, mpi_tile_result_entry, &incoming_count);
      std::cout << rank << " " << incoming_count << std::endl;
      MPI_Recv(recv_buffer, incoming_count, mpi_tile_result_entry,
               MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &mpi_status);
      std::cout << rank << " received" << std::endl;
      for(int i = 0; i < incoming_count; i++) {
         std::cout << recv_buffer[i].to_str() << std::endl;
      }
   }
   else if (rank == 1) {
      //the_tile_result.record({{6, 5, 4},
      //                        {123, 234, 345}});
      std::cout << rank << " size of the vect: " << the_tile_result.vect.size() << std::endl;
      MPI_Send(the_tile_result.vect.data(), the_tile_result.vect.size(), mpi_tile_result_entry,
               0, TAG_ON, MPI_COMM_WORLD);
   }

   MPI_Finalize();
   return 0;
}
               









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

   std::uniform_real_distribution<double> unif(0, 200);
   std::default_random_engine re;

   // Assumption about placement of the regions of columns along the dataset:

   // |          signif         |        nonsignif        |         contrast        |
   // |     fat    |    lean    |     fat    |    lean    |     fat    |    lean    |
   // |homog. tiles|heter. tiles|homog. tiles|heter. tiles|homog. tiles|heter. tiles|
   // R[0]         R[1]         R[2]         R[3]         R[4]         R[5]         R[6]
   // T[R[0]] ...  T[R[1]] ...  T[R[2]] ...  T[R[3]] ...  T[R[4]] ...  T[R[5]] ...  T[R[6]]

   // and the homogeanous tiles are only square-shaped (the margin that You'd expect
   // to be there, is assigned to the lean region).
   // The array T(iles) contains column indeces that are the borders of tiles.
   // In other words, where one tile ends and other begins is given by
   // consequtive values in T.
   // The array R(egions) contains chosen indices of T with the meaning explained above.
   // The array D, not mentioned above, of size T.size - 1, will contain dof's of each tile.
   
   const int tile_width = 100;  
   
   int T[] = {0,   // 0 // fat // signif
              100, // 1 
              200, // 2 // lean
              222, // 3
              235, // 4
              247, // 5
              258, // 6 // fat // nonsignif
              358, // 7
              458, // 8
              558, // 9
              658, // 10
              758, // 11 // lean
              800, // 12
              842, // 13
              876, // 14
              876, // 15 // fat // contrast
              976, // 16
             1076, // 17 // lean
             1200, // 18
             1230, // 19
             1230, // 20
             1245};// 21 // end <= the total column-count

   int R[7] = {0, 2, 6, 11, 15, 17, 21};

   assert(R[2] - R[1] == R[4] - R[3] | R[4] - R[3] == 0);
   assert(R[2] - R[1] == R[6] - R[5] | R[6] - R[5] == 0);

   // number of dicrete buckets per tile
   int D[] = {5, // 0 // fat // signif
              5, // 1 
              5, // 2
              6, // 3 // lean
              7, // 4
              8, // 5
              5, // 6 // fat // nonsignif
              5, // 7
              5, // 8
              5, // 9
              5, // 10
              5, // 11
              6, // 12 // lean
              7, // 13
              8, // 14
              5, // 15 // fat // contrast
              5, // 14
              5, // 15
              6, // 16 // lean
              7, // 17
              8};// 18


   // time-start
   auto start = std::chrono::system_clock::now();

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

   // how many gpu workers?
   // TODO assume there may not be any gpu_workers
   const int kRankTreshold = 2;

   // random seed
   srand(rank);
   
   // output file
   std::string file_name = "out_" + std::to_string(rank) + ".csv";
   std::ofstream file(file_name);  
   
   // buffers for the MPI
   int gpu_tile_index[kDim];
   int cpu_tile_index[kDim];
   
   // an initialy empty tile_result,
   // identity from the point of
   // view of recording the results
   tile_result the_tile_result;

   // Scheduler
   if (rank == 0) {
      // make the 2 schedulers
      Indices_sum& cpu_tile_scheduler = build_cpu_tile_scheduler(R);
      cpu_tile_scheduler.use_buff(cpu_tile_index);
      //
      Indices_sum& gpu_tile_scheduler = build_gpu_tile_scheduler(R);
      gpu_tile_scheduler.use_buff(gpu_tile_index);

      //
      int worker_count = size - 1;
      int incoming_count;
      int rank_recv;

      // Scheduler's receive-buffer
      tile_result_entry* recv_buffer = new tile_result_entry[10 * tile_width * kDim];

      // the main loop
      while (true) {
         std::cout << rank << " about to probe" << std::endl;
         MPI_Probe(MPI_ANY_SOURCE, TAG_ON, MPI_COMM_WORLD, &mpi_status);
         rank_recv = mpi_status.MPI_SOURCE;
         std::cout << rank << " bumped into " << rank_recv << " when probing" << std::endl;
         MPI_Get_count(&mpi_status, mpi_tile_result_entry, &incoming_count);
         std::cout << rank << " preparing to receive " << incoming_count <<"-long array from " << rank_recv << std::endl;
         MPI_Recv(recv_buffer, incoming_count, mpi_tile_result_entry,
                  rank_recv, MPI_ANY_TAG, MPI_COMM_WORLD, &mpi_status);
         //std::cout << rank << " received " << incoming_count << " from " << rank_recv << std::endl;
         if (rank_recv < kRankTreshold) {
            // received from a cpu worker
            if (!cpu_tile_scheduler.get_exhausted()) {
               // send cpu_index to the cpu worker with tag_ON
               MPI_Send(cpu_tile_index, kDim, MPI_INT, rank_recv, TAG_ON, MPI_COMM_WORLD);
               std::cout << rank << " sent " << cpu_tile_scheduler.to_str() << " to cpu_worker " << rank_recv << std::endl;
               cpu_tile_scheduler.up();
            }
            else if (!gpu_tile_scheduler.get_exhausted()) {
               // send gpu_index to cpu worker with tag_ON
               MPI_Send(gpu_tile_index, kDim, MPI_INT, rank_recv, TAG_ON, MPI_COMM_WORLD);
               std::cout << rank << " sent " << gpu_tile_scheduler.to_str() << " to cpu_worker " << rank_recv << std::endl;
               gpu_tile_scheduler.up();
            }
            else{
               // terminate the cpu worker by sending tag_OFF
               std::cout << rank << " about to terminate cpu_worker " << rank_recv << std::endl;
               MPI_Send(gpu_tile_index, kDim, MPI_INT, rank_recv, TAG_OFF, MPI_COMM_WORLD);
               --worker_count;
               std::cout << rank << " number of workers left: " << worker_count << std::endl;
            }
         }
         else {
            // received from a gpu worker
            if (!gpu_tile_scheduler.get_exhausted()) {
               // send gpu_index to the gpu worker with tag_ON
               MPI_Send(gpu_tile_index, kDim, MPI_INT, rank_recv, TAG_ON, MPI_COMM_WORLD);
               std::cout << rank << " sent " << gpu_tile_scheduler.to_str() << " to gpu_worker " << rank_recv << std::endl;
               gpu_tile_scheduler.up();
            }
            else {
               // terminate the gpu worker by sending tag_OFF
               std::cout << rank << " about to terminate gpu_worker " << rank_recv << std::endl;
               MPI_Send(gpu_tile_index, kDim, MPI_INT, rank_recv, TAG_OFF, MPI_COMM_WORLD);
               --worker_count;
            }
         }
         // TODO
         // somehow record the recv_buffer knowing the incoming_count

         if (worker_count == 0) {
            std::cout << rank << " says goodbye" << std::endl;
            break;
         }
      }
   }

   // Workers
   else if (rank < kRankTreshold) {
      // cpu workers
      std::cout << rank << " about to send " << the_tile_result.vect.size() << "-long array" << std::endl;
      MPI_Send(the_tile_result.vect.data(), the_tile_result.vect.size(), mpi_tile_result_entry,
               0, TAG_ON, MPI_COMM_WORLD);
      std::cout << rank << " had sent " << the_tile_result.vect.size() << " elements" << std::endl;
      // the main loop
      while (true) {
         
         std::cout << rank << " about to receive" << std::endl;
         MPI_Recv(cpu_tile_index, kDim, MPI_INT,
                  0, MPI_ANY_TAG, MPI_COMM_WORLD, &mpi_status);
         std::cout << rank << " received ";
         for (int i = 0; i < kDim - 1; i++)
            std::cout << cpu_tile_index[i] <<", ";
         std::cout << cpu_tile_index[kDim - 1] << std::endl;
         
         if (mpi_status.MPI_TAG == TAG_OFF) {
            std::cout << rank << " says goodbye" << std::endl;
            break;
         }


         
         // assume kDim == 3 ?
         // TODO
         // come up with something so that a variable kDim would work
         int column_0_L = T[cpu_tile_index[0]];
         int column_0_R = T[cpu_tile_index[0] + 1];
         int column_1_L = T[cpu_tile_index[1]];
         int column_1_R = T[cpu_tile_index[1] + 1];
         int column_2_L = T[cpu_tile_index[2]];
         int column_2_R = T[cpu_tile_index[2] + 1];
         
         auto start = std::chrono::system_clock::now();
         
         // update the_tile_result object
         for (int i = column_0_L; i < column_0_R; i++) {
            for ( int j = std::max(column_1_L, i); j < column_1_R; j++) {
               for ( int k = std::max(column_2_L, j); k < column_2_R; k++) {
                  // work
                  // kDim random IGs
                  the_tile_result.record({{i, j, k},
                                          {unif(re), unif(re), unif(re)}});
               }
            }
         }
         //std::this_thread::sleep_for(std::chrono::milliseconds(random() % 10));
         auto end = std::chrono::system_clock::now();
         auto time = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
         std::cout << rank << " has worked for " << time.count() << " millisec. on a single tile" << std::endl;

         std::cout << rank << " about to send " << the_tile_result.vect.size() << "-long array" << std::endl;
         MPI_Send(the_tile_result.vect.data(), the_tile_result.vect.size(), mpi_tile_result_entry,
                  0, TAG_ON, MPI_COMM_WORLD);
         std::cout << rank << " had sent " << the_tile_result.vect.size() << " elements" << std::endl;

         // reset the_tile_result object
         the_tile_result.vect.clear();
      }
   }
   else {
      // gpu workers
      std::cout << rank << " about to send " << the_tile_result.vect.size() << "-long array" << std::endl;
      MPI_Send(the_tile_result.vect.data(), the_tile_result.vect.size(), mpi_tile_result_entry,
               0, TAG_ON, MPI_COMM_WORLD);
      std::cout << rank << " had sent " << the_tile_result.vect.size() << " elements" << std::endl;
      // the main loop
      while (true) {
         std::cout << rank << " about to receive" << std::endl;
         MPI_Recv(gpu_tile_index, kDim, MPI_INT,
                  0, MPI_ANY_TAG, MPI_COMM_WORLD, &mpi_status);
         std::cout << rank << " received ";
         for (int i = 0; i < kDim - 1; i++)
            std::cout << gpu_tile_index[i] <<", ";
         std::cout << gpu_tile_index[kDim - 1] << std::endl;
         
         if (mpi_status.MPI_TAG == TAG_OFF) {
            std::cout << rank << " says goodbye" << std::endl;
            break;
         }
         
         // kDim == 3 ?
         int column_0_L = T[gpu_tile_index[0]];
         int column_0_R = T[gpu_tile_index[0] + 1];
         int column_1_L = T[gpu_tile_index[1]];
         int column_1_R = T[gpu_tile_index[1] + 1];
         int column_2_L = T[gpu_tile_index[2]];
         int column_2_R = T[gpu_tile_index[2] + 1];
         
         auto start = std::chrono::system_clock::now();
         
         // update the_tile_result object
         for (int i = column_0_L; i < column_0_R; i++) {
            for ( int j = std::max(column_1_L, i); j < column_1_R; j++) {
               for ( int k = std::max(column_2_L, j); k < column_2_R; k++) {
                  // kDim random IGs
                  the_tile_result.record({{i, j, k},
                                          {unif(re), unif(re), unif(re)}});
               }
            }
         }
         //std::this_thread::sleep_for(std::chrono::milliseconds(random() % 10));
         auto end = std::chrono::system_clock::now();
         auto time = std::chrono::duration_cast<std::chrono::milliseconds>(end - start);
         std::cout << rank << " has worked for " << time.count() << " millisec. on a single tile" << std::endl;

         std::cout << rank << " about to send " << the_tile_result.vect.size() << "-long array" << std::endl;
         MPI_Send(the_tile_result.vect.data(), the_tile_result.vect.size(), mpi_tile_result_entry,
                  0, TAG_ON, MPI_COMM_WORLD);
         std::cout << rank << " had sent " << the_tile_result.vect.size() << " elements" << std::endl;

         // reset the_tile_result object
         the_tile_result.vect.clear();
      }
   }  

   MPI_Finalize();
   return 0;
}





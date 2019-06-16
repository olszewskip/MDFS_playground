// uncomment to disable assert()
// #define NDEBUG

// uncomment to change mode of IG aggregation
// #define AGGREGATE_ABOVE_TRESHOLD_IG
#define AGGREGATE_MAX_IG
#define REMEMBER_CONTEXT

// the dimensionality of k-tuples
#define TUPLE_DIMENSIONALITY 3
const int kDim = TUPLE_DIMENSIONALITY;

#include <mpi.h>
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
#include <algorithm>

#include "Indices.h"
#include "data_structures.h"
//#include "dummy_IG_computation.h"

int main(int argc, char* argv[]) {

   std::uniform_real_distribution<double> unif(0, 200);
   std::default_random_engine re;

   // Assumption about placement of the regions of columns along the dataset (assumed to be a 2D table):

   // |          signif         |        nonsignif        |         contrast        |
   // |     fat    |    lean    |     fat    |    lean    |     fat    |    lean    |
   // |homog. tiles|heter. tiles|homog. tiles|heter. tiles|homog. tiles|heter. tiles|
   // R[0]         R[1]         R[2]         R[3]         R[4]         R[5]         R[6]
   // T[R[0]] ...  T[R[1]] ...  T[R[2]] ...  T[R[3]] ...  T[R[4]] ...  T[R[5]] ...  T[R[6]]

   // and the homogeanous tiles are only square-shaped (the margin that You'd expect
   // to be there, is assigned to the lean region).
   // The array T(ile_bouds) contains column indeces that are the borders of tiles.
   // In other words, where one tile ends and other begins is given by
   // consequtive values in T(ile_bounds).
   // The array R(egion_bounds) contains chosen indices of T(ile_bouds) with the meaning explained above.
   // The array Dofs, not mentioned above, of size "Tile_bounds.size" - 1, will contain the "numbers of degrees
   //  of freedom" (i.e. the bucket count aggregated by multiplication) of each tile.
   
   // total column count of the dataset
   int total_col_count = 1245;

   int Tile_bounds[] = {
                0, // 0 // fat // signif
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
             1245};// 21 // end

   // defined above
   int Region_bounds[7] = {0, 2, 6, 11, 15, 17, 21};
   
   //the last index in Tile_bounds should be the total column count
   assert(total_col_count == Tile_bounds[Region_bounds[6]]);

   // pick up the maximum tile-width along the way
   int max_tile_width = 0;
   for (int i = 0; i < Region_bounds[6]; i++) {
      int tile_width = Tile_bounds[i + 1] - Tile_bounds[i];
      max_tile_width = std::max(max_tile_width, tile_width);
   }
      
   // check that the number of lean tiles matches across lean-regions,
   // assuming that there are lean-nonsignificant- or lean-contrast-regions 
   assert(Region_bounds[2] - Region_bounds[1] == Region_bounds[4] - Region_bounds[3] |
          Region_bounds[4] - Region_bounds[3] == 0);
   assert(Region_bounds[2] - Region_bounds[1] == Region_bounds[6] - Region_bounds[5] |
          Region_bounds[6] - Region_bounds[5] == 0);

   // number of dicrete buckets per tile
   int Dofs[] = {5, // 0 // fat // signif
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
   
   // TODO
   // assert that dofs of fat tiles are all equal
   
   #ifdef AGGREGATE_ABOVE_TRESHOLD_IG
   // we will record all IGs greater than that
   double IG_treshold = 123.;
   // see tile_results_buffer variable below
   int buffer_size_factor = 10;
   #endif

   // // time-measuring
   // std::chrono::high_resolution_clock::time_point time_point_0 = std::chrono::high_resolution_clock::now();
   // // get the elapsed time by doing
   // std::chrono::high_resolution_clock::time_point time_point_1 = std::chrono::high_resolution_clock::now();
   // auto time_diff = time_point_1 - time_point_0;
   // time_point_0 = time_point_1;
   // auto time_diff_millisec = std::chrono::duration_cast<std::chrono::milliseconds>(time_diff);
   // unsigned int milliseconds = time_diff_millisec.count();

   // MPI initialization
   MPI_Init(&argc, &argv);
   int rank, size;
   const int TAG_ON = 1;
   const int TAG_OFF = 0;
   MPI_Comm COMM_0 = MPI_COMM_WORLD;
   MPI_Comm_rank(COMM_0, &rank);
   MPI_Comm_size(COMM_0, &size);
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
   // TODO make sure everything works when there are no gpu_workers
   const int kRankTreshold = 2;

   // random seed
   srand(rank);
   re.seed(rank);
   
   // output file (in case You need it),
   // hopefully closed somewhere at the end
   std::string file_name = "out_" + std::to_string(rank) + ".csv";
   std::ofstream file(file_name);
   //file << "Hello from rank " << rank << "!" << std::endl;
   
   // buffers for the scheduler -> cpu-/gpu- worker communication
   int gpu_tile_index[kDim];
   int cpu_tile_index[kDim];

   // buffer for the worker -> scheduler communication
   #ifdef AGGREGATE_MAX_IG 
   column_result_t* column_results_buffer = new column_result_entry_t[kDim * max_tile_width];
   #endif
   #ifdef AGGREGATE_ABOVE_TRESHOLD_IG
   column_result_entry_t* column_results_buffer = new tile_result_entry_t[kDim * max_tile_width * buffer_size_factor];
   #endif

   
   // an initialy empty tile_result,
   // identity from the point of
   // view of recording the results
   // tile_result the_tile_result;

   // Scheduler
   if (rank == 0) {

      // make the 2 schedulers
      Indices_sum& cpu_tile_scheduler = build_cpu_tile_scheduler(Region_bounds);
      cpu_tile_scheduler.use_buff(cpu_tile_index);
      //
      Indices_sum& gpu_tile_scheduler = build_gpu_tile_scheduler(Region_bounds);
      gpu_tile_scheduler.use_buff(gpu_tile_index);

      // keep track of how many workers are there
      // and what tile each one is working on
      int worker_count = size - 1;
      worker_files files(worker_count);

      // buffer for the final aggregation
      column_result* column_results = new columns_result[total_col_count];
      
      // variables for the MPI calls
      int incoming_count;
      int rank_recv;
      int tag_recv;

      // the main loop
      while (true) {
         std::cout << rank << " about to probe" << std::endl;
         MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, COMM_0, &mpi_status);
         rank_recv = mpi_status.MPI_SOURCE;
         tag_recv = mpi_status.MPI_TAG;
         std::cout << rank << " bumped into " << rank_recv << "with tag " << tag_rec << " when probing" << std::endl;
         MPI_Get_count(&mpi_status, mpi_tile_result_entry, &incoming_count);
         std::cout << rank << " preparing to receive " << incoming_count <<"-long array from " << rank_recv << std::endl;
         MPI_Recv(tile_result_buffer, incoming_count, mpi_tile_result_entry,
                  rank_recv, tag_recv, COMM_0, &mpi_status);
         std::cout << rank << " received " << incoming_count << " from " << rank_recv << std::endl;
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

   file.close();
   MPI_Finalize();
   return 0;
}





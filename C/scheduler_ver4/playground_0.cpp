#define AGGREGATE_MAX_IG
#define REMEMBER_CONTEXT

#define TUPLE_DIMENSIONALITY 3
const int kDim = TUPLE_DIMENSIONALITY;

//#include <mpi.h>
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
std::uniform_real_distribution<double> unif(0, 200);
std::default_random_engine re;
#include "data_structures.h"
//#include "do_dummy_tuple_computation.h"

int main(int argc, char* argv[]) {


   
   int max_tile_width = 3;
   int Tile_bounds[] = {0, 3, 5};
   int tile_index[kDim] = {0, 0, 1};
   
   #ifdef AGGREGATE_MAX_IG 
   column_result_t* column_results_buffer = new column_result_t[kDim * max_tile_width];
   #endif

   tile_handler_t tile_handler(column_results_buffer, tile_index);
   tile_handler.update_tile_info(Tile_bounds);
   tile_handler.print_tile_info();
   tile_handler.clean_buffer();
   tile_handler.print();
   std::cout << " --- " << std::endl;
   tile_handler.compute_tile();
   tile_handler.print(); 

   delete [] column_results_buffer;
   return 0;
}
   




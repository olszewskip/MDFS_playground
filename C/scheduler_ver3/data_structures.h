struct min_work_time_t {
   // remember the minimum among
   // incoming positive integers
   unsigned int min = 0;
   void update(unsigned int& new_min) {
      min = std::min(min, new_min);
   }
};

struct max_work_time_t {
   // remember the maximum among
   // incoming positive integers
   unsigned int max = 0;
   void update(unsigned int& new_max) {
      max = std::max(max, new_max);
   }
};

struct mean_work_time_t {
   // remeber the mean (a double)
   // among incoming positive values
   int count = 0;
   double mean = 0.;
   void update(unsigned int& new_time) {
      ++count;
      mean = (double) count / (count + 1) * mean + new_time / (double)(count + 1);
   }
};

struct worker_file {
   // Presumed to be packed into an array-like
   // structure, and referenced by an index that
   // corresponds to a worker's rank.
   // Serves to hold current information and statistics
   // about what a given worker has worked or is working on.
   // Meant to be used by the scheduler (rank=0 process).
   int current_tile_index[kDim];
   min_work_time_t min_work_time;
   max_work_time_t max_work_time;
   mean_work_time_t mean_work_time;
   std::chrono::high_resolution_clock::time_point last_timestamp;
   
   void record_assignement(int tile_index[]) {
      auto last_timestamp = std::chrono::high_resolution_clock::now();
      for (int i = 0; i < kDim; i++)
        current_tile_index[i] = tile_index[i]; 
   }

   void record_time_diff() {
      auto now = std::chrono::high_resolution_clock::now();
      unsigned int new_work_time =
         std::chrono::duration_cast<std::chrono::milliseconds>(now - last_timestamp).count();
      last_timestamp = now;
      min_work_time.update(new_work_time);
      max_work_time.update(new_work_time);
      mean_work_time.update(new_work_time);
   }
};

struct worker_files {
   worker_file* files = NULL;
   worker_files(int worker_count) {
    files = new worker_file[worker_count];
   }
   worker_file* worker(int rank) { return files + rank - 1; }
};  

#ifdef AGGREGATE_MAX_IG
struct column_result_t {
   // Either just a single double, the "IG",
   // or the IG with (kDim - 1)-long tuple
   // of column indices that were formally
   // conditioned on (the context).
   // The column index, that an instance of
   // this structure implicitly corrensponds
   // to, is not encoded by any of the fields
   // but rather by placement of that instance
   // in an external array-like structure.
   double IG = -1.;
   //
   // Presence of the field "context", the
   // constructor and the "update" method
   // depend on the mode of compilation.
   #ifdef REMEMBER_CONTEXT
   int context[kDim - 1];
   //
   column_result_t() {
      for (int i = 0; i < kDim - 1; i++)
         context[i] = -1;
   }
   // The relevant IG is assumed to sit in a
   // kDim-long array "IGs". "tuple_index[k_index]" is
   // the column that is being described by "IGs[k_index]"
   // and the other column-indices in "tuple_index"
   // form the context.
   void update(int k_index, double IGs[], int tuple_index[]) {
      if (IGs[k_index] > IG) {
         IG = IGs[k_index];
         bool omitted = false;
         for (int i = 0; i < kDim; i++) {
            if (i == k_index) { omitted = true; continue; }
            context[i - omitted] = tuple_index[i];
         }
      }
   }
   // Alternatively the IG may be updated from an existing
   // column_result.
   void update (column_result_t& arg) {
      if (arg.IG > IG) {
         IG = arg.IG;
         for (int i = 0; i < kDim - 1; i++)
            context[i] = arg.context[i];
      }
   }
   #else
   void update (int k_index, double IGs[]) {
      if (IGs[k_index] > IG) { IG = IGs[k_index]; }
   }
   void update (column_result_t& arg) {
      if (arg.IG > IG) { IG = arg.IG; }
   }
   #endif //REMEMBER_CONTEXT
   //
   std::string to_str() {
      std::stringstream ss;
      ss << "IG: " << IG;
      #ifdef REMEMBER_CONTEXT
         ss << "; context: ";
         for (int i = 0; i < kDim - 2; i++)
            ss << context[i] << ", ";
         ss << context[kDim - 2];
      #endif // REMEMBER_CONTEXT
      return ss.str();
   }
};

typedef column_result_t archive_entry_t;

struct tile_handler_t {
   // Class whose main purpose is to correctly
   // update an external buffer with objects
   // of type "column_result_t" based on the IGs
   // computed for each k-tuple of columns in a tile.
   // The slight complication arises from the fact that
   // the (dim_0 x dim_1 x ...) results from the inherently
   // kDim-dimensional tile are being flattened into a
   // (dim_0 + dim_1 + ... )-long array that is the buffer.
   // 
   // The class starts by having two external buffers
   // connected to its internal pointers.
   // Then it functions by
   // a) updating another two of its fields with
   //    information about the current tile (it
   //    uses the "Tile_bounds" array for that),
   // b) cleaning the buffer,
   // c) computing and  recording the kDim IGs per tuple
   //    at right  positions in the external buffer mentioned
   //    at the beginnig. It uses the "compute_tile" funciton
   //    defined outside this class definition.
   // 
   column_result_t* column_results_buffer = NULL;
   int* tile_index = NULL;
   //
   int NW_corner_columns[kDim];
   int SE_corner_columns[kDim];
   int tile_dims[kDim]; // (dim_0, dim_1, ...)
   int buffer_offsets[kDim]; // {0, dim_0, dim_0 + dim_1, ...}
   int buffer_len; // dim_0 + dim_1 + ...
   //
   tile_handler_t (column_result_t* column_results_buffer_arg,
                 int tile_index_arg[]) :
      column_results_buffer(column_results_buffer_arg),
      tile_index(tile_index_arg) {}
   //
   void update_tile_info(int Tile_bounds[]) {
      for (int i = 0; i < kDim; i++) {
         NW_corner_columns[i] = Tile_bounds[tile_index[i]];
         SE_corner_columns[i] = Tile_bounds[tile_index[i] + 1];
         tile_dims[i] = SE_corner_columns[i] - NW_corner_columns[i];
      }
      buffer_offsets[0] = 0;
      for (int i = 1; i < kDim; i++)
         buffer_offsets[i] = buffer_offsets[i - 1] + tile_dims[i - 1];
      buffer_len  = buffer_offsets[kDim - 1] + tile_dims[kDim - 1];
   }
   void print_tile_info() {
      std::cout << "NW corner columns: ";
      for (int i = 0; i < kDim; i++)
         std::cout << NW_corner_columns[i] << ", ";
      std::cout << std::endl;
      std::cout << "SE  corner columns: ";
      for (int i = 0; i < kDim; i++)
         std::cout << SE_corner_columns[i] << ", ";
      std::cout << std::endl;
      std::cout << "tile dims: ";
      for (int i = 0; i < kDim; i++)
         std::cout << tile_dims[i] << ", ";
      std::cout << std::endl;
   } 
   void clean_buffer() {
      // Wipes a starting portion of the "column_results_buffer".
      // Uses the field "buffer_len" set by the "update_tile_info" method
      // and relies on the constructor of "column_result_t".
      for (int buffer_index = 0; buffer_index < buffer_len; buffer_index++)
         column_results_buffer[buffer_index] = column_result_t();
   }
   void compute_tuple(double IGs[], int tuple_index[]);
   //
   void record_tuple(double IGs[], int tuple_index[]) {
      // Put results of a single kDim-dimensional tuple into the buffer.
      int buffer_index;
      for (int dim_index = 0; dim_index < kDim; dim_index++) {
         // This is what the whole fuss with this class
         // is about: computing the buffer_index.
         buffer_index = buffer_offsets[dim_index];
         buffer_index += tuple_index[dim_index] - NW_corner_columns[dim_index];
         column_results_buffer[buffer_index].update(dim_index, IGs
                                                    #ifdef REMEMBER_CONTEXT
                                                    , tuple_index
                                                    #endif
                                                    );
      }
   }
   void compute_tile() {
      //
      double IGs[kDim];
      int tuple_index[kDim];
      #if TUPLE_DIMENSIONALITY == 3
      for (int col_0 = NW_corner_columns[0];
           col_0 < SE_corner_columns[0];
           col_0++) {
         tuple_index[0] = col_0;
         for (int col_1 = std::max(col_0 + 1, NW_corner_columns[1]);
              col_1 < SE_corner_columns[1];
              col_1++) {
            tuple_index[1] = col_1;
            for (int col_2 = std::max(col_1 + 1, NW_corner_columns[2]);
                 col_2 < SE_corner_columns[2];
                 col_2++) {
               tuple_index[2] = col_2;
               compute_tuple(IGs, tuple_index);
               record_tuple(IGs, tuple_index);
            }
         }
      }
      #endif
   }
   void print() {
      // There are two equivalent indices being iterated over here:
      // On one hand the "column_index" is taken from one of the
      // ranges spanned by tile_dims[dim_index] consequtive column
      // indices. On the other hand the "buffer_index" is indepentdently
      // incremented after each "column_index" has been computed.
      int buffer_index = 0;
      for (int dim_index = 0; dim_index < kDim; dim_index++) {
         for (int column_index = NW_corner_columns[dim_index];
              column_index < SE_corner_columns[dim_index];
              column_index++) {
            std::cout << "buffer index: " << buffer_index << std::endl;
            column_result_t& column_result = column_results_buffer[buffer_index];
            std::cout << "column: " << column_index << " "
                      << column_result.to_str() << std::endl;
            ++buffer_index;
         }
      }
   }
   void archive(archive_entry_t archives[]) {
      // A method only for the scheduler (rank=0 process)
      // to iterate over the mpi-received buffer and update
      // the "archives" indexed by all column indeces.
      // See the "print" method which works in similar fashion.
      int buffer_index = 0;
      for (int dim_index = 0; dim_index < kDim; dim_index++) {
         for (int column_index = NW_corner_columns[dim_index];
              column_index < SE_corner_columns[dim_index];
              column_index++) {
            column_result_t& column_result = column_results_buffer[buffer_index];
            archives[column_index].update(column_result);
            ++buffer_index;
         }
      }
   }
};
                      
#endif // AGGREGATE_MAX_IG

//#ifdef AGGREGATE_ABOVE_TRESHOLD_IG
//struct IG_context {
//   double IG;
//   int context[kDim - 1];
//   IG_context(double IG_arg, int context_arg[]) :
//      IG(IG_arg) {
//      for (int i = 0; i < kDim - 1; i++) {
//         context[i] = context_arg[i];
//      }
//   }
//}
//struct tile_result_entry {
//   std::vector<IG_context> IGs_with_contexts;
//   void update(double new_IG, new_context[]) {
//      IGs_with_contexts.emplace_back(new_IG, new_context);
//   }
//}
//#endif


void tile_handler_t::compute_tuple(double IGs[], int tuple_index[]) {
   for (int i = 0; i < kDim; i++)
      IGs[i] = unif(re);
}  


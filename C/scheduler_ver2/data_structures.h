struct tuple_result {
   int columns[kDim];
   double IGs[kDim];
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
   int contexts[kCuriosity][kDim - 1];
   double IGs[kCuriosity];
   tile_result_entry() {}
   tile_result_entry(column_IG_context arg):
      // initialize filling in the first pair
      // of IG and context[] fields
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
   bool absorb (tile_result_entry& arg) {
      if (arg.column != column) { return false; }
      else {
         for (int i = 0; i < kCuriosity; i++)
            update(arg.IGs[i], arg.contexts[i]);
      }
      return true;
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
   void record(tile_result& new_result) {
      for (auto& new_entry : new_result.vect)
         record(new_entry);
   }
   void record(tile_result_entry& new_entry) {
      for (auto& entry : vect) {
         if (entry.absorb(new_entry)) { return; }
      }
      vect.push_back(new_entry);
   }
   void record(column_IG_context arg) {
      for (auto& entry : vect) {
          if (entry.column == arg.column) {
              entry.update(arg.IG, arg.context);
              return;
          }
      }
      vect.emplace_back(arg);
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



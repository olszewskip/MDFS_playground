/*

*/

//#include <iostream>
//#include <chrono>
//#include <thread>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

namespace py = pybind11;

void example0(py::array_t<int>& input) {
    
    py::buffer_info buf = input.request();
    py::print("ptr:", buf.ptr);
    py::print("itemsize:", buf.itemsize);
    py::print("format:", buf.format);
    py::print("ndim:", buf.ndim);
    for (int i = 0; i < buf.shape.size(); i++) {
      py::print(i, "shape:", buf.shape[i]);
    }
    for (int i = 0; i < buf.strides.size(); i++) {
      py::print(i, "stride:", buf.strides[i]);
    }
    
    int *ptr = (int *) buf.ptr;
    int element_count = 1;
    for (auto r: buf.shape) {
      element_count *= r;
    }
    for (int i = 0; i < element_count; i++) {
      py::print(i, "element:", *ptr++);
    }
    
    int NBUCKETS = 100;
    
    lock_t my_lock;
    
    //#pragma omp parallel for
    //for(i=0; i<NBUCKETS; i++){
    //    omp_init_lock(&hist_locks[i]);
    //}


}
        
PYBIND11_MODULE(locks, module) {
    module.def("example0", &example0);
}
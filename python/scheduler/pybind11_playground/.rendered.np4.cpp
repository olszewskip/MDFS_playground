/*

*/
#include <omp.h>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#define NUM_THREADS 2
namespace py = pybind11;

void double_elems(py::array_t<int> input){
     
    py::buffer_info buf = input.request();
    int *ptr = (int *) buf.ptr;
    
    int element_count = 1;
    for (auto r: buf.shape) {
      element_count *= r;
    }
    
    omp_set_num_threads(NUM_THREADS);
    #pragma omp parallel for schedule(static, 1)
    for (int i = 0; i < element_count; i++) {
        ptr[i] *= 2;
    }
     
}
        
PYBIND11_MODULE(np4, module) {
    module.def("double_elems", &double_elems);
}
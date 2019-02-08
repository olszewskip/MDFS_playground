/*

*/
#include <omp.h>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#define NUM_THREADS 2
namespace py = pybind11;

int scalar_prod(py::array_t<int> input1, py::array_t<int> input2){
    
    py::buffer_info buf1 = input1.request();
    py::buffer_info buf2 = input2.request();
    int *ptr1 = (int *) buf1.ptr;
    int *ptr2 = (int *) buf2.ptr;
    
    int element_count = 1;
    for (auto r: buf1.shape) {
      element_count *= r;
    }
    
    omp_set_num_threads(NUM_THREADS);
    
    int prod = 0;
    
    #pragma omp parallel
    {
        #pragma omp for reduction(+:prod)
        for (int i = 0; i < element_count; i++) {
            prod += ptr1[i] * ptr2[i];
        }
    }
    
    return prod;
}
        
PYBIND11_MODULE(np6, module) {
    module.def("scalar_prod", &scalar_prod);
}
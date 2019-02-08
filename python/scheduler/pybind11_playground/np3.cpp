/*
<%
setup_pybind11(cfg)
cfg['compiler_args'] = ['-fopenmp']
cfg['linker_args'] = ['-fopenmp']
%>
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
    #pragma omp parallel
    {
        int id = omp_get_thread_num();
        int nthrds = omp_get_num_threads();
        for (int i = id; i < element_count; i += nthrds) {
            ptr[i] *= 2;
        }
    }
     
}
        
PYBIND11_MODULE(np3, module) {
    module.def("double_elems", &double_elems);
}
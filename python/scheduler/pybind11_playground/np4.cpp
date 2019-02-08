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

void scalar_prod(int N, py::array_t<int> input1, py::array_t<int> input2){
     
    py::buffer_info buf1 = input1.request();
    py::buffer_info buf2 = input2.request();
    int *ptr1 = (int *) buf1.ptr;
    int *ptr2 = (int *) buf2.ptr;
    
    
    omp_set_num_threads(NUM_THREADS);
    
    int prod = 0;
    
    #pragma omp parallel
    {
        local_prod = 0;
        #pragma omp for
        for (int i = 0; i < N; i++) {
            local_prod += ptr1[i] + ptr2[i];
        }
        
        #pragma omp atomic
        prod += local_prod
    }
     
}
        
PYBIND11_MODULE(np4, module) {
    module.def("scalar_prod", &scalar_prod);
}
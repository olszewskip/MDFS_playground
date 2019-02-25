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

namespace py = pybind11;

#define NUM_THREADS 16

double scalar_prod(py::array_t<double> input1, py::array_t<double> input2){
    
    py::buffer_info buf1 = input1.request();
    py::buffer_info buf2 = input2.request();
    auto *ptr1 = (double *) buf1.ptr;
    auto *ptr2 = (double *) buf2.ptr;
    
    int element_count = 1;
    for (auto r: buf1.shape) {
      element_count *= r;
    }
    
    omp_set_num_threads(NUM_THREADS);
    
    double prod = 0.;
    
    int nthreads;
    #pragma omp parallel
    {
        int nthrds = omp_get_num_threads();
        int id = omp_get_thread_num();
        if (id==0) nthreads = nthrds;
        #pragma omp for reduction(+:prod)
        for (int i = 0; i < element_count; i++) {
            prod += ptr1[i] * ptr2[i];
        }
    }
    
    py::print("I got", nthreads, "threads!");
    
    return prod;
}
        
PYBIND11_MODULE(dummy_work1, module) {
    module.def("scalar_prod", &scalar_prod);
}
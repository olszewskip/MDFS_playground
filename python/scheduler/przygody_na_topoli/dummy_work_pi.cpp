/*
<%
setup_pybind11(cfg)
cfg['compiler_args'] = ['-std=c++11', '-fopenmp']
cfg['linker_args'] = ['-fopenmp']
%>
*/
#include <omp.h>
#include <pybind11/pybind11.h>

namespace py = pybind11;
using namespace pybind11::literals;


double omp_pi(int n, int rank=0, int size=1) {
    
    double step = 1.0 / n;
    double pi = 0;
    
    int nthreads = 0;
    
    #pragma omp parallel
    {
        int id = omp_get_thread_num();
        if (id==0) nthreads = omp_get_num_threads();
        
        double x;
        #pragma omp for reduction(+:pi)
        for(int i = rank ; i < n ; i+=size){
            x = (i + 0.5) * step;
            pi += 4.0 / (1 + x*x);
        }
    }    
    
    py::print("I got", nthreads, "threads!");
    
    pi *= step;
    return pi;
}

        
PYBIND11_MODULE(dummy_work_pi, module) {
    module.def("omp_pi", &omp_pi, "n"_a, "rank"_a=0, "size"_a=1);
}

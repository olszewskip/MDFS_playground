/*
<%
setup_pybind11(cfg)
cfg['compiler_args'] = ['-fopenmp']
cfg['linker_args'] = ['-fopenmp']
%>
*/
#include <pybind11/pybind11.h>
#include <omp.h>

namespace py = pybind11;

double compute_pi(int n) {
    
    double step = 1.0/n;
    double pi = 0;
    
    omp_set_num_threads(2);
    #pragma omp parallel
    {
        double x;
        #pragma omp for reduction(+:pi)
        for(int i = 0 ; i < n ; i++){
            x = (i + 0.5) * step;
            pi += 4.0 / (1 + x*x);
        }
    }
    
    pi *= step;
    return pi;
}


PYBIND11_MODULE(openmp_pi, m) {
    m.def("compute_pi", &compute_pi);
}
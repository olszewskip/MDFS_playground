/*
<%
cfg['compiler_args'] = ['-fopenmp']
cfg['linker_args'] = ['-fopenmp']
setup_pybind11(cfg)
%>
*/

#include <omp.h>
#include <pybind11/pybind11.h>

namespace py = pybind11;

double sequential(int n) {
    
    double step = 1.0 / n;
    double pi = 0;

    double x;
    for(int i = 0 ; i < n ; i++){
        x = (i + 0.5) * step;
        pi += 4.0 / (1 + x*x);
    }
    
    pi *= step;
    return pi;
}

// #define NUM_THREADS 1

double parallel(int n) {
    
    double step = 1.0 / n;
    double pi = 0;
    
    int nthreads = 0;

    //omp_set_num_threads(NUM_THREADS);
    
    #pragma omp parallel
    {
        int id = omp_get_thread_num();
        if (id==0) nthreads = omp_get_num_threads();
        
        double x;
        #pragma omp for reduction(+:pi)
        for(int i = 0 ; i < n ; i++){
            x = (i + 0.5) * step;
            pi += 4.0 / (1 + x*x);
        }
    }    
    
    py::print("I got", nthreads, "threads!");
    
    pi *= step;
    return pi;
}


PYBIND11_MODULE(my_pi, m) {
    m.def("sequential", &sequential);
    m.def("parallel", &parallel);
}
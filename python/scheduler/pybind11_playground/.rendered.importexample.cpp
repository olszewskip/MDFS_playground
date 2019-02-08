/*

*/
#include <pybind11/pybind11.h>
#include <omp.h>

namespace py = pybind11;

double do_something(int n) {
    
    double step = 1.0/n;
    double pi = 0;
    
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


PYBIND11_MODULE(importexample, m) {
    m.def("do_something", &do_something);
}
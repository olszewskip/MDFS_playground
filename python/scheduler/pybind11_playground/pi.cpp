/*
<%
setup_pybind11(cfg)
%>
*/
#include <pybind11/pybind11.h>

namespace py = pybind11;

double compute_pi(int n) {
    
    double step = 1.0/n;
    double pi = 0;

    double x;
    for(int i = 0 ; i < n ; i++){
        x = (i + 0.5) * step;
        pi += 4.0 / (1 + x*x);
    }
    
    pi *= step;
    return pi;
}


PYBIND11_MODULE(pi, m) {
    m.def("compute_pi", &compute_pi);
}
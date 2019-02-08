/*

*/
#include <pybind11/pybind11.h>
#include <omp.h>

namespace py = pybind11;

int do_something(int n) {
    int m = 13;
    #pragma omp parallel
    {
    for(int i=0; i<n; i++) {
        m += 1;
    }
    }
    return m;
}


PYBIND11_MODULE(importexample2, m) {
    m.def("do_something", &do_something, py::call_guard<py::gil_scoped_release>());
}
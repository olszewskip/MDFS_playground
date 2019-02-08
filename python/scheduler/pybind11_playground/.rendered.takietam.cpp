/*

*/
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
namespace py = pybind11;

void sayhi() {
    py::print("hi!");
}
        
PYBIND11_MODULE(takietam, module) {
    module.def("sayhi", &sayhi);
}
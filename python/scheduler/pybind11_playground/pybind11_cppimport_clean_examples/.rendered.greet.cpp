/*

*/

#include <pybind11/pybind11.h>

namespace py = pybind11;

void sayhi() {
    py::print("hi!");
}
        
PYBIND11_MODULE(greet, module) {
    module.def("sayhi", &sayhi);
}
/*
<%
setup_pybind11(cfg)
cfg['linker_args'] = ['my_math.cpp']
%>
*/

#include <pybind11/pybind11.h>
#include "my_math.hpp"

namespace py = pybind11;
using namespace pybind11::literals;

PYBIND11_MODULE(wrap_my_math, m) {
    m.doc() = "The best mathematical library in the universe";
    m.def("add", &add, "A function which adds two numbers", "i"_a, "j"_a=0);
    m.def("multiply", &multiply, "A function which multiplies two numbers", "i"_a, "j"_a=1);
}
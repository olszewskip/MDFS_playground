/*

*/

#include <pybind11/pybind11.h>
#include "my_add.hpp"

namespace py = pybind11;
using namespace pybind11::literals;

PYBIND11_MODULE(my_add_wrap, m) {
    m.doc() = "pybind11 example plugin";
    m.def("add", &add, "A function which adds two numbers", "i"_a, "j"_a=2);
}
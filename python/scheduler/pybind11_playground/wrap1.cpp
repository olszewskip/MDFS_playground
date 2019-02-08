#include <pybind11/pybind11.h>
#include "funcs.hpp"

namespace py = pybind11;
using namespace pybind11::literals;

PYBIND11_MODULE(example, m) {
    m.doc() = "pybind11 example plugin"; // optional module docstring

    m.def("add", &add, "A function which adds two numbers",
          "i"_a, "j"_a=2);
}
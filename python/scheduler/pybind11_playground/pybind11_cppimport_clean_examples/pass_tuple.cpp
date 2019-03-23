/*
<%
setup_pybind11(cfg)
%>
*/

#include <tuple>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

namespace py = pybind11;

void print_tuple(std::tuple<int, int> tuple_) {
    
    py::print(tuple_);
    py::print(std::get<0>(tuple_), std::get<1>(tuple_));

}
        
PYBIND11_MODULE(pass_tuple, module) {
    module.def("print_tuple", &print_tuple);
}
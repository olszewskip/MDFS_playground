/*

*/
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
namespace py = pybind11;

void print_ndim(py::array_t<double> input1) {
    py::buffer_info buf1 = input1.request();
    int a = 1;
    py::print(buf1.ndim);
}
        
PYBIND11_MODULE(ndim, module) {
    module.def("print_ndim", &print_ndim);
}
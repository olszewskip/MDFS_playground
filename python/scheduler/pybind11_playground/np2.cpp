/*
<%
setup_pybind11(cfg)
%>
*/
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
namespace py = pybind11;

void double_elems(py::array_t<double> input) {
    
    py::buffer_info buf = input.request();
    double *ptr = (double *) buf.ptr;
    
    int element_count = 1;
    for (auto r: buf.shape) {
      element_count *= r;
    }
    for (int i = 0; i < element_count; i++) {
        *ptr++ *= 2.0;
    }
}
        
PYBIND11_MODULE(np2, module) {
    module.def("double_elems", &double_elems);
}
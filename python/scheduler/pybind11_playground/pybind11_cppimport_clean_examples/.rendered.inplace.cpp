/*

*/
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
namespace py = pybind11;

void twice(py::array_t<double> input) {
    
    py::buffer_info buf = input.request();
    auto *ptr = (double *) buf.ptr;
    
    int element_count = 1;
    for (auto r: buf.shape) {
      element_count *= r;
    }
    
    for (int i = 0; i < element_count; i++) {
        *ptr++ *= 2;
    }
}
        
PYBIND11_MODULE(inplace, module) {
    module.def("twice", &twice);
}
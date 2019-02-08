/*

*/
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
namespace py = pybind11;

void double_elems(py::array_t<int> input) {
    
    py::buffer_info buf = input.request();
    int *ptr = (int *) buf.ptr;
    
    int element_count = 1;
    for (auto r: buf.shape) {
      element_count *= r;
    }
    for (int i = 0; i < element_count; i++) {
        *ptr++ *= 2;
    }
}
        
PYBIND11_MODULE(np2, module) {
    module.def("double_elems", &double_elems);
}
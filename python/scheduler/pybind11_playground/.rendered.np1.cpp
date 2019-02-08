/*

*/
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
namespace py = pybind11;

void print_info(py::array_t<int> input1) {
    py::buffer_info buf1 = input1.request();
    py::print("ptr:", buf1.ptr);
    py::print("itemsize:", buf1.itemsize);
    py::print("format:", buf1.format);
    py::print("ndim:", buf1.ndim);
    for (int i = 0; i < buf1.shape.size(); i++) {
      py::print(i, "shape:", buf1.shape[i]);
    }
    for (int i = 0; i < buf1.strides.size(); i++) {
      py::print(i, "stride:", buf1.strides[i]);
    }
    
    int *ptr = (int *) buf1.ptr;
    int element_count = 1;
    for (auto r: buf1.shape) {
      element_count *= r;
    }
    for (int i = 0; i < element_count; i++) {
        py::print(i, "element:", *ptr++);
    }

}
        
PYBIND11_MODULE(np1, module) {
    module.def("print_info", &print_info);
}
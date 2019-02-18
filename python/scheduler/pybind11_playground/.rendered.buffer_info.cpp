/*

*/

//#include <iostream>
//#include <chrono>
//#include <thread>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

namespace py = pybind11;

void print_info(py::array_t<int>& input) {
    
    py::buffer_info buf = input.request();
    py::print("ptr:", buf.ptr);
    py::print("itemsize:", buf.itemsize);
    py::print("format:", buf.format);
    py::print("ndim:", buf.ndim);
    for (int i = 0; i < buf.shape.size(); i++) {
      py::print(i, "shape:", buf.shape[i]);
    }
    for (int i = 0; i < buf.strides.size(); i++) {
      py::print(i, "stride:", buf.strides[i]);
    }
    
    int *ptr = (int *) buf.ptr;
    int element_count = 1;
    for (auto r: buf.shape) {
      element_count *= r;
    }
    for (int i = 0; i < element_count; i++) {
      py::print(i, "element:", *ptr++);
    }
    
//    std::this_thread::sleep_for(3s);

}
        
PYBIND11_MODULE(buffer_info, module) {
    module.def("print_info", &print_info);
}
/*

*/

#include <iostream>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
namespace py = pybind11;

void print_dict(py::dict dict) {
    
    // print to console
    for (auto item : dict)
        std::cout << "key=" << std::string(py::str(item.first)) << ", "
                  << "value=" << std::string(py::str(item.second)) << std::endl;
}

void py_print_dict(py::dict dict) {
    
    for (auto item : dict)
        py::print(item.first, "=>", item.second);
}

void py_print_map_as_dict(std::map<std::string, int> &dict) {
    
    for (auto item : dict)
        py::print(item.first, "=>", item.second);
}
        
PYBIND11_MODULE(pass_dict, module) {
    module.def("print_dict", &print_dict);
    module.def("py_print_dict", &py_print_dict);
    module.def("py_print_map_as_dict", &py_print_map_as_dict);
}
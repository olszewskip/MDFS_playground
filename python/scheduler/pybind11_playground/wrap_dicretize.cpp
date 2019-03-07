/*
<%
setup_pybind11(cfg)
cfg['linker_args'] = ['discretize.cpp']
%>
*/

#include <pybind11/pybind11.h>
#include "my_math.hpp"

namespace py = pybind11;
using namespace pybind11::literals;

PYBIND11_MODULE(wrap_dicretize, module) {
    module.doc() = "Wrapper of https://bitbucket.org/mdfs/mdfs/src/master/src/cpu/discretize.cpp";
    module.def("discretize", &discretize,
               "In place discretization. ",
               "seed"_a,
               "discretization_index"_a,
               "feature_id"_a,
               "divisions"_a,
               "object_count"_a,
               "in_data"_a,
               "sorted_in_data"_a,
               "out_data"_a,
               "range"_a
              );
}
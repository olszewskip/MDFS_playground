
#ifndef DISCRETIZE_H
#define DISCRETIZE_H

#include <cstddef>
#include <cstdint>
#include <vector>

#include <pybind11/pybind11.h>
//#include <pybind11/stl.h>
#include <pybind11/numpy.h>

namespace py = pybind11;

void discretize(
    uint32_t seed,
    uint32_t discretization_index,
    uint32_t feature_id,
    std::size_t divisions,
    std::size_t object_count,
    py::array_t<double> py_in_data,
    const std::vector<double>& sorted_in_data, // Note: why different type?
    uint8_t* out_data,
    double range
);

#endif
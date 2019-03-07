/*

*/

#include <algorithm>
#include <cstdint>
#include <cstring>
#include <random>
#include <vector>
#include <pybind11/pybind11.h>
//#include <pybind11/stl.h>
#include <pybind11/numpy.h>

namespace py = pybind11;
using namespace pybind11::literals;

void discretize(
    uint32_t seed,
    uint32_t discretization_index,
    uint32_t feature_id,
    std::size_t divisions,
    std::size_t object_count,
    py::array_t<double> py_in_data,
    py::array_t<double> py_sorted_in_data,
    py::array_t<uint8_t> py_out_data,
    double range_
) {
    
    // the python part
    py::buffer_info py_in_data_buf = py_in_data.request();
    auto *in_data = static_cast<const double *>(py_in_data_buf.ptr);    
    
    py::buffer_info py_sorted_in_data_buf = py_sorted_in_data.request();
    auto *sorted_in_data = static_cast<const double *>(py_sorted_in_data_buf.ptr);
    
    py::buffer_info py_out_data_buf = py_out_data.request();
    auto *out_data = static_cast<uint8_t *>(py_out_data_buf.ptr);
    
    // end of the python part
    
    double* thresholds = new double[divisions];

    // brackets to limit scope
    {
        double sum = 0.0f;
        // brackets to limit scope
        {
            std::mt19937 seed_random_generator0(seed);
            std::mt19937 seed_random_generator1(seed_random_generator0() ^ discretization_index);
            std::mt19937 random_generator(seed_random_generator1() ^ feature_id);

            // E(X) = (a + b) / 2 = (1 - range + 1 + range) / 2 = 1
            std::uniform_real_distribution<double> uniform_range(1.0f - range_, 1.0f + range_);

            for (std::size_t d = 0; d < divisions; ++d) {
                thresholds[d] = uniform_range(random_generator);
                sum += thresholds[d];
            }

            sum += uniform_range(random_generator);
        }

        std::size_t done = 0;
        const double length_step = static_cast<double>(object_count) / sum;

        // thresholds are converted from an arbitrary space into real values (via indices)
        // d - iterates over divisions (of a variable in a discretization)
        for (std::size_t d = 0; d < divisions; ++d) {
            done += std::lround(thresholds[d] * length_step);

            // Note: Check when will this happen, maybe could be skipped
            if (done >= object_count) {
                done = object_count - 1;
            }

            thresholds[d] = sorted_in_data[done];
        }
    }

    // o - iterates over objects
    for (std::size_t o = 0; o < object_count; ++o) {
        out_data[o] = 0;

        // out_data[o] (starting with 0) is incremented every time in_data[o] is above a threshold
        // divisions is a small number (<=15), no reason to use binsearch, hence linear
        // d - iterates over divisions (per object o)
        for (std::size_t d = 0; d < divisions; ++d) {
            out_data[o] += in_data[o] > thresholds[d];
        }
    }

    delete[] thresholds;
}


PYBIND11_MODULE(wrap_discretize, module) {
    module.doc() = "Wrapper of https://bitbucket.org/mdfs/mdfs/src/master/src/cpu/discretize.cpp";
    module.def("discretize", &discretize,
               "In place discretization.",
               "seed"_a,
               "discretization_index"_a,
               "feature_id"_a,
               "divisions"_a,
               "object_count"_a,
               "py_in_data"_a,
               "py_sorted_in_data"_a,
               "py_out_data"_a,
               "range_"_a
              );
}
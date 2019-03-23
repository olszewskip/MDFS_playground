/*
<%
cfg['compiler_args'] = ['-std=c++11', '-fopenmp']
cfg['linker_args'] = ['-fopenmp']
setup_pybind11(cfg)
%>
*/

#include <math.h>
#include <tuple>
//#include <vector>
#include <pybind11/pybind11.h>
//#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include <omp.h>

#include "matrix_class.h"

namespace py = pybind11;


std::tuple<double, double, double> work_3a(const int kData_dim,
                                           const std::tuple<int, int, int> py_Xdims,
                                           py::array_t<int> &py_X0,
                                           py::array_t<int> &py_X1,
                                           py::array_t<int> &py_X2,
                                           const int kN_classes,
                                           py::array_t<double> &py_pseudo_counts,
                                           py::array_t<int> &py_y) {
    
    const int kC_Xdim_0 = std::get<0>(py_Xdims);
    const int kC_Xdim_1 = std::get<1>(py_Xdims);
    const int kC_Xdim_2 = std::get<2>(py_Xdims);
    const int kC_ydim = kN_classes;
    
    py::buffer_info py_X0_buf = py_X0.request();
    auto *X0 = static_cast<int *>(py_X0_buf.ptr);
    py::buffer_info py_X1_buf = py_X1.request();
    auto *X1 = static_cast<int *>(py_X1_buf.ptr);
    py::buffer_info py_X2_buf = py_X2.request();
    auto *X2 = static_cast<int *>(py_X2_buf.ptr);
    
    py::buffer_info py_pseudo_counts_buf = py_pseudo_counts.request();
    auto *pseudo_counts = (double *) py_pseudo_counts_buf.ptr;
    
    py::buffer_info py_y_buf = py_y.request();
    auto *y = (int *) py_y_buf.ptr;
    
    Matrix<int> contingency_m(kC_Xdim_0, kC_Xdim_1, kC_Xdim_2, kC_ydim);
    for (int data_idx = 0; data_idx < kData_dim; data_idx++) {
        ++contingency_m(X0[data_idx], X1[data_idx], X2[data_idx], y[data_idx]);
    }
    
    double neg_H = 0., neg_H_X0 = 0., neg_H_X1 = 0., neg_H_X2 = 0.;
    
    #pragma omp parallel
    {
    #pragma omp sections
    {
        #pragma omp section
        {
        double local_neg_H = 0.;
        for (int C_Xidx_0 = 0; C_Xidx_0 < kC_Xdim_0; ++C_Xidx_0) {
            for (int C_Xidx_1 = 0; C_Xidx_1 < kC_Xdim_1; ++C_Xidx_1) {
                for (int C_Xidx_2 = 0; C_Xidx_2 < kC_Xdim_2; ++C_Xidx_2) {
                    double count_y = 0.;
                    for (int C_yidx = 0; C_yidx < kC_ydim; ++C_yidx) {
                        double count = pseudo_counts[C_yidx] + contingency_m(C_Xidx_0, C_Xidx_1, C_Xidx_2, C_yidx);
                        local_neg_H += count * log2(count);
                        count_y += count;
                    }
                    local_neg_H -= count_y * log2(count_y);
                }
            }
        }
        #pragma omp atomic
            neg_H += local_neg_H;
        }        
        
        #pragma omp section
        {
        double local_neg_H_X0 = 0.;
        for (int C_Xidx_1 = 0; C_Xidx_1 < kC_Xdim_1; C_Xidx_1++) {
            for (int C_Xidx_2 = 0; C_Xidx_2 < kC_Xdim_2; C_Xidx_2++) {
                double count_y_X0 = 0.;
                for (int C_yidx = 0; C_yidx < kC_ydim; C_yidx++) {
                    double count_X0 = 0.;
                    for (int C_Xidx_0 = 0; C_Xidx_0 < kC_Xdim_0; C_Xidx_0++) {
                        double count = pseudo_counts[C_yidx] + contingency_m(C_Xidx_0, C_Xidx_1, C_Xidx_2, C_yidx);
                        count_X0 += count;
                    }
                    local_neg_H_X0 += count_X0 * log2(count_X0);
                    count_y_X0 += count_X0;
                }
                local_neg_H_X0 -= count_y_X0 * log2(count_y_X0);
            }
        }
        #pragma omp atomic
            neg_H_X0 += local_neg_H_X0;
        }
        
        #pragma omp section
        {
        double local_neg_H_X1 = 0.;
        for (int C_Xidx_0 = 0; C_Xidx_0 < kC_Xdim_0; C_Xidx_0++) {
            for (int C_Xidx_2 = 0; C_Xidx_2 < kC_Xdim_2; C_Xidx_2++) {
                double count_y_X1 = 0.;
                for (int C_yidx = 0; C_yidx < kC_ydim; C_yidx++) {
                    double count_X1 = 0.;
                    for (int C_Xidx_1 = 0; C_Xidx_1 < kC_Xdim_1; C_Xidx_1++) {
                        double count = pseudo_counts[C_yidx] + contingency_m(C_Xidx_0, C_Xidx_1, C_Xidx_2, C_yidx);
                        count_X1 += count;
                    }
                    local_neg_H_X1 += count_X1 * log2(count_X1);
                    count_y_X1 += count_X1;
                }
                local_neg_H_X1 -= count_y_X1 * log2(count_y_X1);
            }
        }
        #pragma omp atomic
            neg_H_X1 += local_neg_H_X1;
        }
        
        #pragma omp section
        {
        double local_neg_H_X2 = 0.;
        for (int C_Xidx_0 = 0; C_Xidx_0 < kC_Xdim_0; C_Xidx_0++) {
            for (int C_Xidx_1 = 0; C_Xidx_1 < kC_Xdim_1; C_Xidx_1++) {
                double count_y_X2 = 0.;
                for (int C_yidx = 0; C_yidx < kC_ydim; C_yidx++) {
                    double count_X2 = 0.;
                    for (int C_Xidx_2 = 0; C_Xidx_2 < kC_Xdim_2; C_Xidx_2++) {
                        double count = pseudo_counts[C_yidx] + contingency_m(C_Xidx_0, C_Xidx_1, C_Xidx_2, C_yidx);
                        count_X2 += count;
                    }
                    local_neg_H_X2 += count_X2 * log2(count_X2);
                    count_y_X2 += count_X2;
                }
                local_neg_H_X2 -= count_y_X2 * log2(count_y_X2);
            }
        }
        #pragma omp atomic
            neg_H_X2 += local_neg_H_X2;
        }        
    }
    }
    
    return std::make_tuple(neg_H - neg_H_X0,
                           neg_H - neg_H_X1,
                           neg_H - neg_H_X2
                          ); 
}
       

PYBIND11_MODULE(fast, module) {
    module.def("work_3a", &work_3a);
}
/*
<%
setup_pybind11(cfg)
cfg['compiler_args'] = ['-std=c++17', '-fopenmp']
cfg['linker_args'] = ['-fopenmp']
%>
*/

#include <math.h>
#include <tuple>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

namespace py = pybind11;

std::tuple<double, double> work(const int kData_dim,
                                const int kDivisions,
                                py::array_t<int> &py_X0,
                                py::array_t<int> &py_X1,
                                const int kN_classes,
                                py::array_t<int> &py_y) {
    
    const int kC_Xdim = kDivisions + 1;
    const int kC_ydim = kN_classes;
    
    py::buffer_info py_X0_buf = py_X0.request();
    auto *X0 = (int *) py_X0_buf.ptr;
    py::buffer_info py_X1_buf = py_X1.request();
    auto *X1 = (int *) py_X1_buf.ptr;
    
    py::buffer_info py_y_buf = py_y.request();
    auto *y = (int *) py_y_buf.ptr;
    
    int contingency_m[kC_Xdim + 1][kC_Xdim + 1][kC_ydim + 1] = {};

    for (int data_idx = 0; data_idx < kData_dim; data_idx++) {
        py::print("data_idx=", data_idx);
        py::print("X0[data_idx]=", X0[data_idx]);
        py::print("X1[data_idx]=", X1[data_idx]);
        py::print("y[data_idx]=", y[data_idx]);
        contingency_m[X0[data_idx]][X1[data_idx]][y[data_idx]] += 1;
        contingency_m[X0[data_idx]][X1[data_idx]][kC_ydim] += 1;
        contingency_m[X0[data_idx]][kC_Xdim][y[data_idx]] += 1;
        contingency_m[kC_Xdim][X1[data_idx]][y[data_idx]] += 1;
        contingency_m[X0[data_idx]][kC_Xdim][kC_ydim] += 1;
        contingency_m[kC_Xdim][X1[data_idx]][kC_ydim] += 1;
        py::print("--");
        py::print(contingency_m[0][0][0], contingency_m[0][1][0], contingency_m[0][2][0]);
        py::print(contingency_m[1][0][0], contingency_m[1][1][0], contingency_m[1][2][0]);
        py::print(contingency_m[2][0][0], contingency_m[2][1][0], contingency_m[2][2][0]);
        py::print("--");
        py::print(contingency_m[0][0][1], contingency_m[0][1][1], contingency_m[0][2][1]);
        py::print(contingency_m[1][0][1], contingency_m[1][1][1], contingency_m[1][2][1]);
        py::print(contingency_m[2][0][1], contingency_m[2][1][1], contingency_m[2][2][1]);
        py::print("--");
        py::print(contingency_m[0][0][2], contingency_m[0][1][2], contingency_m[0][2][2]);
        py::print(contingency_m[1][0][2], contingency_m[1][1][2], contingency_m[1][2][2]);
        py::print(contingency_m[2][0][2], contingency_m[2][1][2], contingency_m[2][2][2]);
        py::print("----");

    }
    double neg_H = 0, neg_H_y = 0, neg_H_X0 = 0, neg_H_X1 = 0;
    

    for (int C_Xidx_0 = 0; C_Xidx_0 < kC_Xdim; C_Xidx_0++) {
        for (int C_Xidx_1 = 0; C_Xidx_1 < kC_Xdim; C_Xidx_1++) {
          for (int C_yidx = 0; C_yidx < kC_ydim; C_yidx++) {
            int &count = contingency_m[C_Xidx_0][C_Xidx_1][C_yidx];
            neg_H += count * log2(count);
              
          }
          int &count_y = contingency_m[C_Xidx_0][C_Xidx_1][kC_ydim];
          neg_H -= count_y * log2(count_y);
        }
    }


    for (int C_Xidx = 0; C_Xidx < kC_Xdim; C_Xidx++) {
        for (int C_yidx = 0; C_yidx < kC_ydim; C_yidx++) {
          int &count_X0 = contingency_m[kC_Xdim][C_Xidx][C_yidx];
          neg_H_X0 += count_X0 * log2(count_X0);
          int &count_X1 = contingency_m[C_Xidx][kC_Xdim][C_yidx];
          neg_H_X1 += count_X1 * log2(count_X1);
        }
        int &count_X0_y = contingency_m[kC_Xdim][C_Xidx][kC_ydim];
        neg_H_X0 -= count_X0_y * log2(count_X0_y);
        int &count_X1_y = contingency_m[C_Xidx][kC_Xdim][kC_ydim];
        neg_H_X1 -= count_X1_y * log2(count_X1_y);
    }


    return {neg_H - neg_H_X0, neg_H - neg_H_X1};
}
        
PYBIND11_MODULE(singleton_work, module) {
    module.def("work", &work);
}
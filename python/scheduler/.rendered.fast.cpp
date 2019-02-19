/*

*/

#include <math.h>
#include <tuple>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

namespace py = pybind11;

std::tuple<double, double> work_2(const int kData_dim,
                                  const int kDivisions,
                                  py::array_t<int> &py_X0,
                                  py::array_t<int> &py_X1,
                                  const int kN_classes,
                                  py::array_t<double> &py_pseudo_counts,
                                  py::array_t<int> &py_y) {
    
    const int kC_Xdim = kDivisions + 1;
    const int kC_ydim = kN_classes;
    
    py::buffer_info py_X0_buf = py_X0.request();
    auto *X0 = (int *) py_X0_buf.ptr;
    py::buffer_info py_X1_buf = py_X1.request();
    auto *X1 = (int *) py_X1_buf.ptr;
    
    py::buffer_info py_y_buf = py_y.request();
    auto *y = (int *) py_y_buf.ptr;
    
    py::buffer_info py_pseudo_counts_buf = py_pseudo_counts.request();
    auto *pseudo_counts = (double *) py_pseudo_counts_buf.ptr;
    
    double contingency_m[kC_Xdim + 1][kC_Xdim + 1][kC_ydim + 1] = {};

    for (int data_idx = 0; data_idx < kData_dim; data_idx++) {
        contingency_m[X0[data_idx]][X1[data_idx]][y[data_idx]] += 1.;
        contingency_m[X0[data_idx]][X1[data_idx]][kC_ydim] += 1.;
        contingency_m[X0[data_idx]][kC_Xdim][y[data_idx]] += 1.;
        contingency_m[kC_Xdim][X1[data_idx]][y[data_idx]] += 1.;
        contingency_m[X0[data_idx]][kC_Xdim][kC_ydim] += 1.;
        contingency_m[kC_Xdim][X1[data_idx]][kC_ydim] += 1.;
    }
    
    double neg_H = 0., neg_H_y = 0., neg_H_X0 = 0., neg_H_X1 = 0.;

    for (int C_Xidx_0 = 0; C_Xidx_0 < kC_Xdim; C_Xidx_0++) {
        for (int C_Xidx_1 = 0; C_Xidx_1 < kC_Xdim; C_Xidx_1++) {
          double count_y = contingency_m[C_Xidx_0][C_Xidx_1][kC_ydim];            
          for (int C_yidx = 0; C_yidx < kC_ydim; C_yidx++) {
            count_y += pseudo_counts[C_yidx];  
            double count = pseudo_counts[C_yidx] + contingency_m[C_Xidx_0][C_Xidx_1][C_yidx];
            neg_H += count * log2(count);  
          }
          neg_H -= count_y * log2(count_y);  
        }
    }

    for (int C_Xidx = 0; C_Xidx < kC_Xdim; C_Xidx++) {
        double count_X0_y = contingency_m[kC_Xdim][C_Xidx][kC_ydim];
        double count_X1_y = contingency_m[C_Xidx][kC_Xdim][kC_ydim];
        for (int C_yidx = 0; C_yidx < kC_ydim; C_yidx++) {
          count_X0_y += kC_Xdim * pseudo_counts[C_yidx];
          count_X1_y += kC_Xdim * pseudo_counts[C_yidx];
            
          double count_X0 = kC_Xdim * pseudo_counts[C_yidx] + contingency_m[kC_Xdim][C_Xidx][C_yidx];
          neg_H_X0 += count_X0 * log2(count_X0);
            
          double count_X1 = kC_Xdim * pseudo_counts[C_yidx] + contingency_m[C_Xidx][kC_Xdim][C_yidx];
          neg_H_X1 += count_X1 * log2(count_X1);
        }
        neg_H_X0 -= count_X0_y * log2(count_X0_y);
        neg_H_X1 -= count_X1_y * log2(count_X1_y);
    }

    return {neg_H - neg_H_X0, neg_H - neg_H_X1};
}
    
    
std::tuple<double, double, double> work_3(const int kData_dim,
                                          const int kDivisions,
                                          py::array_t<int> &py_X0,
                                          py::array_t<int> &py_X1,
                                          py::array_t<int> &py_X2,
                                          const int kN_classes,
                                          py::array_t<double> &py_pseudo_counts,
                                          py::array_t<int> &py_y) {
    
    const int kC_Xdim = kDivisions + 1;
    const int kC_ydim = kN_classes;
    
    py::buffer_info py_X0_buf = py_X0.request();
    auto *X0 = (int *) py_X0_buf.ptr;
    py::buffer_info py_X1_buf = py_X1.request();
    auto *X1 = (int *) py_X1_buf.ptr;
    py::buffer_info py_X2_buf = py_X2.request();
    auto *X2 = (int *) py_X2_buf.ptr;
    
    py::buffer_info py_pseudo_counts_buf = py_pseudo_counts.request();
    auto *pseudo_counts = (double *) py_pseudo_counts_buf.ptr;
    
    py::buffer_info py_y_buf = py_y.request();
    auto *y = (int *) py_y_buf.ptr;
    
    double contingency_m[kC_Xdim + 1][kC_Xdim + 1][kC_Xdim + 1][kC_ydim + 1] = {};

    for (int data_idx = 0; data_idx < kData_dim; data_idx++) {
        contingency_m[X0[data_idx]][X1[data_idx]][X2[data_idx]][y[data_idx]] += 1.;
        contingency_m[X0[data_idx]][X1[data_idx]][X2[data_idx]][kC_ydim] += 1.;
        contingency_m[X0[data_idx]][X1[data_idx]][kC_Xdim][y[data_idx]] += 1.;
        contingency_m[X0[data_idx]][kC_Xdim][X1[data_idx]][y[data_idx]] += 1.;
        contingency_m[kC_Xdim][X1[data_idx]][X2[data_idx]][y[data_idx]] += 1.;
        contingency_m[X0[data_idx]][X1[data_idx]][kC_Xdim][kC_ydim] += 1.;
        contingency_m[X0[data_idx]][kC_Xdim][X2[data_idx]][kC_ydim] += 1.;
        contingency_m[kC_Xdim][X0[data_idx]][X2[data_idx]][kC_ydim] += 1.;
    }
    
    double neg_H = 0., neg_H_y = 0., neg_H_X0 = 0., neg_H_X1 = 0., neg_H_X2 = 0.;
    
    for (int C_Xidx_0 = 0; C_Xidx_0 < kC_Xdim; C_Xidx_0++) {
        for (int C_Xidx_1 = 0; C_Xidx_1 < kC_Xdim; C_Xidx_1++) {
            for (int C_Xidx_2 = 0; C_Xidx_2 < kC_Xdim; C_Xidx_2++) {
              double count_y = contingency_m[C_Xidx_0][C_Xidx_1][C_Xidx_2][kC_ydim];            
              for (int C_yidx = 0; C_yidx < kC_ydim; C_yidx++) {
                count_y += pseudo_counts[C_yidx];  
                double count = pseudo_counts[C_yidx] + contingency_m[C_Xidx_0][C_Xidx_1][C_Xidx_2][C_yidx];
                neg_H += count * log2(count);  
              }
              neg_H -= count_y * log2(count_y);
            }
        }
    }

    for (int C_Xidx_i = 0; C_Xidx_i < kC_Xdim; C_Xidx_i++) {
        for (int C_Xidx_j = 0; C_Xidx_j < kC_Xdim; C_Xidx_j++) {
            double count_X0_y = contingency_m[kC_Xdim][C_Xidx_i][C_Xidx_j][kC_ydim];
            double count_X1_y = contingency_m[C_Xidx_i][kC_Xdim][C_Xidx_j][kC_ydim];
            double count_X2_y = contingency_m[C_Xidx_i][C_Xidx_j][kC_Xdim][kC_ydim];
            for (int C_yidx = 0; C_yidx < kC_ydim; C_yidx++) {
              count_X0_y += kC_Xdim * pseudo_counts[C_yidx];
              count_X1_y += kC_Xdim * pseudo_counts[C_yidx];
              count_X2_y += kC_Xdim * pseudo_counts[C_yidx];

              double count_X0 = kC_Xdim * pseudo_counts[C_yidx] + contingency_m[kC_Xdim][C_Xidx_i][C_Xidx_j][C_yidx];
              neg_H_X0 += count_X0 * log2(count_X0);

              double count_X1 = kC_Xdim * pseudo_counts[C_yidx] + contingency_m[C_Xidx_i][kC_Xdim][C_Xidx_j][C_yidx];
              neg_H_X1 += count_X1 * log2(count_X1);
                
              double count_X2 = kC_Xdim * pseudo_counts[C_yidx] + contingency_m[C_Xidx_i][C_Xidx_j][kC_Xdim][C_yidx];
              neg_H_X2 += count_X2 * log2(count_X2);
            }
            neg_H_X0 -= count_X0_y * log2(count_X0_y);
            neg_H_X1 -= count_X1_y * log2(count_X1_y);
            neg_H_X2 -= count_X2_y * log2(count_X2_y);
        }
    }

    return {neg_H - neg_H_X0, neg_H - neg_H_X1, neg_H - neg_H_X2};
}

        

PYBIND11_MODULE(fast, module) {
    module.def("work_2", &work_2);
    module.def("work_3", &work_3);
}
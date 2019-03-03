/*
<%
cfg['compiler_args'] = ['-std=c++11', '-fopenmp']
cfg['linker_args'] = ['-fopenmp']
setup_pybind11(cfg)
%>
*/

#include <math.h>
#include <tuple>
#include <vector>
#include <pybind11/pybind11.h>
//#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include <omp.h>

namespace py = pybind11;

std::tuple<double, double, double> work_3a(const int kData_dim,
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
    auto *X0 = static_cast<int *>(py_X0_buf.ptr);
    py::buffer_info py_X1_buf = py_X1.request();
    auto *X1 = static_cast<int *>(py_X1_buf.ptr);
    py::buffer_info py_X2_buf = py_X2.request();
    auto *X2 = static_cast<int *>(py_X2_buf.ptr);
    
    py::buffer_info py_y_buf = py_y.request();
    auto *y = static_cast<int *>(py_y_buf.ptr);
    
    py::buffer_info py_pseudo_counts_buf = py_pseudo_counts.request();
    auto *pseudo_counts = static_cast<double *>(py_pseudo_counts_buf.ptr);
    
    int contingency_m[kC_Xdim][kC_Xdim][kC_Xdim][kC_ydim] = {};
    double contingency_m_y[kC_Xdim][kC_Xdim][kC_Xdim] = {};

    for (int data_idx = 0; data_idx < kData_dim; data_idx++) {
        contingency_m[X0[data_idx]][X1[data_idx]][X2[data_idx]][y[data_idx]]++;
    }

    double neg_H = 0., neg_H_X0 = 0., neg_H_X1 = 0., neg_H_X2 = 0.;
    
    #pragma omp parallel
    {
    
    #pragma omp for reduction (+: neg_H, neg_H_X0, neg_H_X1, neg_H_X2)
    for (int C_Xidx_i = 0; C_Xidx_i < kC_Xdim; C_Xidx_i++) {
        for (int C_Xidx_j = 0; C_Xidx_j < kC_Xdim; C_Xidx_j++) {
            for (int C_yidx = 0; C_yidx < kC_ydim; C_yidx++) {
                double count_X0 = 0;
                double count_X1 = 0;
                double count_X2 = 0;
                for (int C_Xidx_k = 0; C_Xidx_k < kC_Xdim; C_Xidx_k++) {
                    double count_ijk = pseudo_counts[C_yidx] + contingency_m[C_Xidx_i][C_Xidx_j][C_Xidx_k][C_yidx];
                    count_X0 += pseudo_counts[C_yidx] + contingency_m[C_Xidx_k][C_Xidx_i][C_Xidx_j][C_yidx];
                    count_X1 += pseudo_counts[C_yidx] + contingency_m[C_Xidx_i][C_Xidx_k][C_Xidx_j][C_yidx];
                    count_X2 += count_ijk;
                    contingency_m_y[C_Xidx_i][C_Xidx_j][C_Xidx_k] += count_ijk;
                    neg_H += count_ijk * log2(count_ijk);
                }
                neg_H_X0 += count_X0 * log2(count_X0);
                neg_H_X1 += count_X1 * log2(count_X1);
                neg_H_X2 += count_X2 * log2(count_X2);
            }
        }
    }

    #pragma omp for reduction (+: neg_H, neg_H_X0, neg_H_X1, neg_H_X2)
    for (int C_Xidx_i = 0; C_Xidx_i < kC_Xdim; C_Xidx_i++) {
        for (int C_Xidx_j = 0; C_Xidx_j < kC_Xdim; C_Xidx_j++) {
            double count_X0_y = 0;
            double count_X1_y = 0;
            double count_X2_y = 0;
            for (int C_Xidx_k = 0; C_Xidx_k < kC_Xdim; C_Xidx_k++) {
                neg_H -= contingency_m_y[C_Xidx_i][C_Xidx_j][C_Xidx_k] * log2(contingency_m_y[C_Xidx_i][C_Xidx_j][C_Xidx_k]);
                count_X0_y += contingency_m_y[C_Xidx_k][C_Xidx_i][C_Xidx_j];
                count_X1_y += contingency_m_y[C_Xidx_i][C_Xidx_k][C_Xidx_j];
                count_X2_y += contingency_m_y[C_Xidx_i][C_Xidx_j][C_Xidx_k];
            }
            neg_H_X0 -= count_X0_y * log2(count_X0_y);
            neg_H_X1 -= count_X1_y * log2(count_X1_y);
            neg_H_X2 -= count_X2_y * log2(count_X2_y);
        }
    }
    
    }
        
    return std::make_tuple(neg_H - neg_H_X0,
                           neg_H - neg_H_X1,
                           neg_H - neg_H_X2
                          );
}

    

std::tuple<double, double, double> work_3b(const int kData_dim,
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
    auto *X0 = static_cast<int *>(py_X0_buf.ptr);
    py::buffer_info py_X1_buf = py_X1.request();
    auto *X1 = static_cast<int *>(py_X1_buf.ptr);
    py::buffer_info py_X2_buf = py_X2.request();
    auto *X2 = static_cast<int *>(py_X2_buf.ptr);
    
    py::buffer_info py_y_buf = py_y.request();
    auto *y = static_cast<int *>(py_y_buf.ptr);
    
    py::buffer_info py_pseudo_counts_buf = py_pseudo_counts.request();
    auto *pseudo_counts = static_cast<double *>(py_pseudo_counts_buf.ptr);
    
    int contingency_m[kC_Xdim][kC_Xdim][kC_Xdim][kC_ydim] = {};
    double contingency_m_y[kC_Xdim][kC_Xdim][kC_Xdim] = {};
        
    omp_lock_t contingency_m_lock[kC_Xdim][kC_Xdim][kC_Xdim][kC_ydim];
                                                                                                
    for (int C_Xidx_0 = 0; C_Xidx_0 < kC_Xdim; C_Xidx_0++) {
        for (int C_Xidx_1 = 0; C_Xidx_1 < kC_Xdim; C_Xidx_1++) {
            for (int C_Xidx_2 = 0; C_Xidx_2 < kC_Xdim; C_Xidx_2++) {
                for (int C_yidx = 0; C_yidx < kC_ydim; C_yidx++) {
                    omp_init_lock(&contingency_m_lock[C_Xidx_0][C_Xidx_1][C_Xidx_2][C_yidx]);
                }
            }
        }
    }
    
    double neg_H = 0., neg_H_X0 = 0., neg_H_X1 = 0., neg_H_X2 = 0.;
    
    #pragma omp parallel
    {
    
    #pragma omp for                                           
    for (int data_idx = 0; data_idx < kData_dim; data_idx++) {
        omp_set_lock(&contingency_m_lock[X0[data_idx]][X1[data_idx]][X2[data_idx]][y[data_idx]]);
        contingency_m[X0[data_idx]][X1[data_idx]][X2[data_idx]][y[data_idx]]++;
        omp_unset_lock(&contingency_m_lock[X0[data_idx]][X1[data_idx]][X2[data_idx]][y[data_idx]]);
    }
    
    #pragma omp for reduction (+: neg_H, neg_H_X0, neg_H_X1, neg_H_X2)
    for (int C_Xidx_i = 0; C_Xidx_i < kC_Xdim; C_Xidx_i++) {
        for (int C_Xidx_j = 0; C_Xidx_j < kC_Xdim; C_Xidx_j++) {
            for (int C_yidx = 0; C_yidx < kC_ydim; C_yidx++) {
                double count_X0 = 0;
                double count_X1 = 0;
                double count_X2 = 0;
                for (int C_Xidx_k = 0; C_Xidx_k < kC_Xdim; C_Xidx_k++) {
                    double count_ijk = pseudo_counts[C_yidx] + contingency_m[C_Xidx_i][C_Xidx_j][C_Xidx_k][C_yidx];
                    count_X0 += pseudo_counts[C_yidx] + contingency_m[C_Xidx_k][C_Xidx_i][C_Xidx_j][C_yidx];
                    count_X1 += pseudo_counts[C_yidx] + contingency_m[C_Xidx_i][C_Xidx_k][C_Xidx_j][C_yidx];
                    count_X2 += count_ijk;
                    contingency_m_y[C_Xidx_i][C_Xidx_j][C_Xidx_k] += count_ijk;
                    neg_H += count_ijk * log2(count_ijk);
                }
                neg_H_X0 += count_X0 * log2(count_X0);
                neg_H_X1 += count_X1 * log2(count_X1);
                neg_H_X2 += count_X2 * log2(count_X2);
            }
        }
    }

    #pragma omp for reduction (+: neg_H, neg_H_X0, neg_H_X1, neg_H_X2)
    for (int C_Xidx_i = 0; C_Xidx_i < kC_Xdim; C_Xidx_i++) {
        for (int C_Xidx_j = 0; C_Xidx_j < kC_Xdim; C_Xidx_j++) {
            double count_X0_y = 0;
            double count_X1_y = 0;
            double count_X2_y = 0;
            for (int C_Xidx_k = 0; C_Xidx_k < kC_Xdim; C_Xidx_k++) {
                neg_H -= contingency_m_y[C_Xidx_i][C_Xidx_j][C_Xidx_k] * log2(contingency_m_y[C_Xidx_i][C_Xidx_j][C_Xidx_k]);
                count_X0_y += contingency_m_y[C_Xidx_k][C_Xidx_i][C_Xidx_j];
                count_X1_y += contingency_m_y[C_Xidx_i][C_Xidx_k][C_Xidx_j];
                count_X2_y += contingency_m_y[C_Xidx_i][C_Xidx_j][C_Xidx_k];
            }
            neg_H_X0 -= count_X0_y * log2(count_X0_y);
            neg_H_X1 -= count_X1_y * log2(count_X1_y);
            neg_H_X2 -= count_X2_y * log2(count_X2_y);
        }
    }
    
    }
    
    for (int C_Xidx_0 = 0; C_Xidx_0 < kC_Xdim; C_Xidx_0++) {
        for (int C_Xidx_1 = 0; C_Xidx_1 < kC_Xdim; C_Xidx_1++) {
            for (int C_Xidx_2 = 0; C_Xidx_2 < kC_Xdim; C_Xidx_2++) {
                for (int C_yidx = 0; C_yidx < kC_ydim; C_yidx++) {
                    omp_destroy_lock(&contingency_m_lock[C_Xidx_0][C_Xidx_1][C_Xidx_2][C_yidx]);
                }
            }
        }
    }
        
    return std::make_tuple(neg_H - neg_H_X0,
                           neg_H - neg_H_X1,
                           neg_H - neg_H_X2
                          );
}

std::tuple<double, double, double> work_3c(const int kData_dim,
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
    auto *X0 = static_cast<int *>(py_X0_buf.ptr);
    py::buffer_info py_X1_buf = py_X1.request();
    auto *X1 = static_cast<int *>(py_X1_buf.ptr);
    py::buffer_info py_X2_buf = py_X2.request();
    auto *X2 = static_cast<int *>(py_X2_buf.ptr);
    
    py::buffer_info py_y_buf = py_y.request();
    auto *y = static_cast<int *>(py_y_buf.ptr);
    
    py::buffer_info py_pseudo_counts_buf = py_pseudo_counts.request();
    auto *pseudo_counts = static_cast<double *>(py_pseudo_counts_buf.ptr);
    
    int contingency_m[kC_Xdim][kC_Xdim][kC_Xdim][kC_ydim] = {};
    double contingency_m_y[kC_Xdim][kC_Xdim][kC_Xdim] = {};
    
    double neg_H = 0., neg_H_X0 = 0., neg_H_X1 = 0., neg_H_X2 = 0.;
    
    #pragma omp parallel
    {
    #pragma omp for
    for (int data_idx = 0; data_idx < kData_dim; data_idx++) {
        #pragma omp atomic
            contingency_m[X0[data_idx]][X1[data_idx]][X2[data_idx]][y[data_idx]]++;
    }
    
    #pragma omp for reduction (+: neg_H, neg_H_X0, neg_H_X1, neg_H_X2)
    for (int C_Xidx_i = 0; C_Xidx_i < kC_Xdim; C_Xidx_i++) {
        for (int C_Xidx_j = 0; C_Xidx_j < kC_Xdim; C_Xidx_j++) {
            for (int C_yidx = 0; C_yidx < kC_ydim; C_yidx++) {
                double count_X0 = 0;
                double count_X1 = 0;
                double count_X2 = 0;
                for (int C_Xidx_k = 0; C_Xidx_k < kC_Xdim; C_Xidx_k++) {
                    double count_ijk = pseudo_counts[C_yidx] + contingency_m[C_Xidx_i][C_Xidx_j][C_Xidx_k][C_yidx];
                    count_X0 += pseudo_counts[C_yidx] + contingency_m[C_Xidx_k][C_Xidx_i][C_Xidx_j][C_yidx];
                    count_X1 += pseudo_counts[C_yidx] + contingency_m[C_Xidx_i][C_Xidx_k][C_Xidx_j][C_yidx];
                    count_X2 += count_ijk;
                    contingency_m_y[C_Xidx_i][C_Xidx_j][C_Xidx_k] += count_ijk;
                    neg_H += count_ijk * log2(count_ijk);
                }
                neg_H_X0 += count_X0 * log2(count_X0);
                neg_H_X1 += count_X1 * log2(count_X1);
                neg_H_X2 += count_X2 * log2(count_X2);
            }
        }
    }

    #pragma omp for reduction (+: neg_H, neg_H_X0, neg_H_X1, neg_H_X2)
    for (int C_Xidx_i = 0; C_Xidx_i < kC_Xdim; C_Xidx_i++) {
        for (int C_Xidx_j = 0; C_Xidx_j < kC_Xdim; C_Xidx_j++) {
            double count_X0_y = 0;
            double count_X1_y = 0;
            double count_X2_y = 0;
            for (int C_Xidx_k = 0; C_Xidx_k < kC_Xdim; C_Xidx_k++) {
                neg_H -= contingency_m_y[C_Xidx_i][C_Xidx_j][C_Xidx_k] * log2(contingency_m_y[C_Xidx_i][C_Xidx_j][C_Xidx_k]);
                count_X0_y += contingency_m_y[C_Xidx_k][C_Xidx_i][C_Xidx_j];
                count_X1_y += contingency_m_y[C_Xidx_i][C_Xidx_k][C_Xidx_j];
                count_X2_y += contingency_m_y[C_Xidx_i][C_Xidx_j][C_Xidx_k];
            }
            neg_H_X0 -= count_X0_y * log2(count_X0_y);
            neg_H_X1 -= count_X1_y * log2(count_X1_y);
            neg_H_X2 -= count_X2_y * log2(count_X2_y);
        }
    }
        
    }

        
    return std::make_tuple(neg_H - neg_H_X0,
                           neg_H - neg_H_X1,
                           neg_H - neg_H_X2
                          );
}
        

PYBIND11_MODULE(fast11, module) {
    module.def("work_3a", &work_3a);
    module.def("work_3b", &work_3b);
    module.def("work_3c", &work_3c);
}
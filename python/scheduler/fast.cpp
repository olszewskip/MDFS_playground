/*
<%
setup_pybind11(cfg)
cfg['compiler_args'] = ['-std=c++17', '-fopenmp']
cfg['linker_args'] = ['-fopenmp']
%>
*/

#include <math.h>
#include <tuple>
#include <vector>
#include <pybind11/pybind11.h>
//#include <pybind11/stl.h>
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
    auto *X0 = static_cast<int *>(py_X0_buf.ptr);
    py::buffer_info py_X1_buf = py_X1.request();
    auto *X1 = static_cast<int *>(py_X1_buf.ptr);
    
    py::buffer_info py_y_buf = py_y.request();
    auto *y = static_cast<int *>(py_y_buf.ptr);
    
    py::buffer_info py_pseudo_counts_buf = py_pseudo_counts.request();
    auto *pseudo_counts = static_cast<double *>(py_pseudo_counts_buf.ptr);
    
    int contingency_m[kC_Xdim][kC_Xdim][kC_ydim] = {};
    double contingency_m_y[kC_Xdim][kC_Xdim] = {};
    
    for (int data_idx = 0; data_idx < kData_dim; data_idx++) {
        contingency_m[X0[data_idx]][X1[data_idx]][y[data_idx]]++;
    }
    
    double neg_H = 0., neg_H_X0 = 0., neg_H_X1 = 0.;
    double neg_H_y=0., neg_H_X0_y=0., neg_H_X1_y=0.;
    
    //#pragma omp parallel for reduction (+: neg_H, neg_H_X0, neg_H_X1)
    for (int C_Xidx_i = 0; C_Xidx_i < kC_Xdim; C_Xidx_i++) {
        for (int C_yidx = 0; C_yidx < kC_ydim; C_yidx++) {
            double count_X0 = 0;
            double count_X1 = 0;
            for (int C_Xidx_j = 0; C_Xidx_j < kC_Xdim; C_Xidx_j++) {
                double count_ij = pseudo_counts[C_yidx] + contingency_m[C_Xidx_i][C_Xidx_j][C_yidx];
                count_X0 += pseudo_counts[C_yidx] + contingency_m[C_Xidx_j][C_Xidx_i][C_yidx];
                count_X1 += count_ij;
                contingency_m_y[C_Xidx_i][C_Xidx_j] += count_ij;
                neg_H += count_ij * log2(count_ij);
            }
            neg_H_X0 += count_X0 * log2(count_X0);
            neg_H_X1 += count_X1 * log2(count_X1);
        }
    }
    
    //#pragma omp parallel for reduction (+: neg_H, neg_H_X0, neg_H_X1)
    for (int C_Xidx_i = 0; C_Xidx_i < kC_Xdim; C_Xidx_i++) {
        double count_X0_y = 0;
        double count_X1_y = 0;        
        for (int C_Xidx_j = 0; C_Xidx_j < kC_Xdim; C_Xidx_j++) {
            neg_H -= contingency_m_y[C_Xidx_i][C_Xidx_j] * log2(contingency_m_y[C_Xidx_i][C_Xidx_j]);
            count_X0_y += contingency_m_y[C_Xidx_j][C_Xidx_i];
            count_X1_y += contingency_m_y[C_Xidx_i][C_Xidx_j];
        }
        neg_H_X0 -= count_X0_y * log2(count_X0_y);
        neg_H_X1 -= count_X1_y * log2(count_X1_y);
    }
    
    return {neg_H - neg_H_X0,
            neg_H - neg_H_X1
           };
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
    
    //#pragma omp parallel for reduction (+: neg_H, neg_H_X0, neg_H_X1, neg_H_X2)
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

    //#pragma omp parallel for reduction (+: neg_H, neg_H_X0, neg_H_X1, neg_H_X2)
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
    
    return {neg_H - neg_H_X0,
            neg_H - neg_H_X1,
            neg_H - neg_H_X2
           };
}
    

std::tuple<double, double> work_2b(const int kData_dim,
                                  const int kDivisions,
                                  py::array_t<int> &py_X0,
                                  py::array_t<int> &py_X1,
                                  const int kN_classes,
                                  py::array_t<double> &py_pseudo_counts,
                                  py::array_t<int> &py_y) {
    
    const int kC_Xdim = kDivisions + 1;
    const int kC_ydim = kN_classes;
    
    py::buffer_info py_X0_buf = py_X0.request();
    auto *X0 = static_cast<int *>(py_X0_buf.ptr);
    py::buffer_info py_X1_buf = py_X1.request();
    auto *X1 = static_cast<int *>(py_X1_buf.ptr);
    
    py::buffer_info py_y_buf = py_y.request();
    auto *y = static_cast<int *>(py_y_buf.ptr);
    
    py::buffer_info py_pseudo_counts_buf = py_pseudo_counts.request();
    auto *pseudo_counts = static_cast<double *>(py_pseudo_counts_buf.ptr);
    
    int contingency_m[kC_Xdim + 1][kC_Xdim + 1][kC_ydim + 1] = {};
    

    for (int data_idx = 0; data_idx < kData_dim; data_idx++) {
        contingency_m[X0[data_idx]][X1[data_idx]][y[data_idx]]++;
        contingency_m[X0[data_idx]][X1[data_idx]][kC_ydim]++;
        contingency_m[X0[data_idx]][kC_Xdim][y[data_idx]]++;
        contingency_m[kC_Xdim][X1[data_idx]][y[data_idx]]++;
        contingency_m[X0[data_idx]][kC_Xdim][kC_ydim]++;
        contingency_m[kC_Xdim][X1[data_idx]][kC_ydim]++;
    }
    
    double neg_H = 0., neg_H_X0 = 0., neg_H_X1 = 0.;

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
    
    py::buffer_info py_pseudo_counts_buf = py_pseudo_counts.request();
    auto *pseudo_counts = (double *) py_pseudo_counts_buf.ptr;
    
    py::buffer_info py_y_buf = py_y.request();
    auto *y = (int *) py_y_buf.ptr;
    
    int contingency_m[kC_Xdim + 1][kC_Xdim + 1][kC_Xdim + 1][kC_ydim + 1] = {};

    for (int data_idx = 0; data_idx < kData_dim; data_idx++) {
        contingency_m[X0[data_idx]][X1[data_idx]][X2[data_idx]][y[data_idx]]++;
        contingency_m[X0[data_idx]][X1[data_idx]][X2[data_idx]][kC_ydim]++;
        contingency_m[X0[data_idx]][X1[data_idx]][kC_Xdim][y[data_idx]]++;
        contingency_m[X0[data_idx]][kC_Xdim][X2[data_idx]][y[data_idx]]++;
        contingency_m[kC_Xdim][X1[data_idx]][X2[data_idx]][y[data_idx]]++;
        contingency_m[X0[data_idx]][X1[data_idx]][kC_Xdim][kC_ydim]++;
        contingency_m[X0[data_idx]][kC_Xdim][X2[data_idx]][kC_ydim]++;
        contingency_m[kC_Xdim][X1[data_idx]][X2[data_idx]][kC_ydim]++;
    }
    
    double neg_H = 0., neg_H_X0 = 0., neg_H_X1 = 0., neg_H_X2 = 0.;
    
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
    
    py::buffer_info py_pseudo_counts_buf = py_pseudo_counts.request();
    auto *pseudo_counts = (double *) py_pseudo_counts_buf.ptr;
    
    py::buffer_info py_y_buf = py_y.request();
    auto *y = (int *) py_y_buf.ptr;
    
    int contingency_m[kC_Xdim][kC_Xdim][kC_Xdim][kC_ydim] = {};
    for (int data_idx = 0; data_idx < kData_dim; data_idx++) {
        contingency_m[X0[data_idx]][X1[data_idx]][X2[data_idx]][y[data_idx]]++;
    }
    
    double neg_H = 0., neg_H_X0 = 0., neg_H_X1 = 0., neg_H_X2 = 0.;
    
    #pragma omp parallel
    {
    #pragma omp sections
    {
        #pragma omp section
        {
        double local_neg_H = 0.;
        for (int C_Xidx_0 = 0; C_Xidx_0 < kC_Xdim; C_Xidx_0++) {
            for (int C_Xidx_1 = 0; C_Xidx_1 < kC_Xdim; C_Xidx_1++) {
                for (int C_Xidx_2 = 0; C_Xidx_2 < kC_Xdim; C_Xidx_2++) {
                    double count_y = 0.;
                    for (int C_yidx = 0; C_yidx < kC_ydim; C_yidx++) {
                        double count = pseudo_counts[C_yidx] + contingency_m[C_Xidx_0][C_Xidx_1][C_Xidx_2][C_yidx];
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
        for (int C_Xidx_1 = 0; C_Xidx_1 < kC_Xdim; C_Xidx_1++) {
            for (int C_Xidx_2 = 0; C_Xidx_2 < kC_Xdim; C_Xidx_2++) {
                double count_y_X0 = 0.;
                for (int C_yidx = 0; C_yidx < kC_ydim; C_yidx++) {
                    double count_X0 = 0.;
                    for (int C_Xidx_0 = 0; C_Xidx_0 < kC_Xdim; C_Xidx_0++) {
                        double count = pseudo_counts[C_yidx] + contingency_m[C_Xidx_0][C_Xidx_1][C_Xidx_2][C_yidx];
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
        for (int C_Xidx_0 = 0; C_Xidx_0 < kC_Xdim; C_Xidx_0++) {
            for (int C_Xidx_2 = 0; C_Xidx_2 < kC_Xdim; C_Xidx_2++) {
                double count_y_X1 = 0.;
                for (int C_yidx = 0; C_yidx < kC_ydim; C_yidx++) {
                    double count_X1 = 0.;
                    for (int C_Xidx_1 = 0; C_Xidx_1 < kC_Xdim; C_Xidx_1++) {
                        double count = pseudo_counts[C_yidx] + contingency_m[C_Xidx_0][C_Xidx_1][C_Xidx_2][C_yidx];
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
        for (int C_Xidx_0 = 0; C_Xidx_0 < kC_Xdim; C_Xidx_0++) {
            for (int C_Xidx_1 = 0; C_Xidx_1 < kC_Xdim; C_Xidx_1++) {
                double count_y_X2 = 0.;
                for (int C_yidx = 0; C_yidx < kC_ydim; C_yidx++) {
                    double count_X2 = 0.;
                    for (int C_Xidx_2 = 0; C_Xidx_2 < kC_Xdim; C_Xidx_2++) {
                        double count = pseudo_counts[C_yidx] + contingency_m[C_Xidx_0][C_Xidx_1][C_Xidx_2][C_yidx];
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
    
    return {neg_H - neg_H_X0, neg_H - neg_H_X1, neg_H - neg_H_X2};
}
        

PYBIND11_MODULE(fast, module) {
    module.def("work_2", &work_2);
    module.def("work_2b", &work_2b);
    module.def("work_3", &work_3);
    module.def("work_3b", &work_3b);
    module.def("work_3c", &work_3c);
}
/*

*/
#include <cmath>
#include <omp.h>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

#define NUM_THREADS 2
namespace py = pybind11;

int contingency(const int divisions, const int N, py::array_t<int> input1, py::array_t<int> input2){
    
    double l = std::log(3);
    py::print(l);
    
    int contingency_m[divisions][divisions] = {};
    
    for (int i = 0; i < divisions; i++)
        for (int j = 0; j < divisions; j++)
            py::print(contingency_m[i][j]);

    
    py::buffer_info buf1 = input1.request();
    py::buffer_info buf2 = input2.request();
    int *ptr1 = (int *) buf1.ptr;
    int *ptr2 = (int *) buf2.ptr;
    
    
    omp_set_num_threads(NUM_THREADS);
    
    int prod = 0;
    
    #pragma omp parallel for reduction(+:prod)
    for (int i = 0; i < N; i++) {
        prod += ptr1[i] * ptr2[i];
    }
    
    return prod;
}
        
PYBIND11_MODULE(np5, module) {
    module.def("contingency", &contingency);
}
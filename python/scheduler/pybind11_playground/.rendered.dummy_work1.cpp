/*

*/
#include <omp.h>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

namespace py = pybind11;

// define NUM_THREADS 2

double parallel_pi(int n) {
    
    double step = 1.0 / n;
    double pi = 0;

    #pragma omp parallel
    {
        double x;
        #pragma omp for reduction(+:pi)
        for(int i = 0 ; i < n ; i++){
            x = (i + 0.5) * step;
            pi += 4.0 / (1 + x*x);
        }
    }    
    
    pi *= step;
    return pi;
}


double scalar_prod(py::array_t<double> input1, py::array_t<double> input2){
    
    py::buffer_info buf1 = input1.request();
    py::buffer_info buf2 = input2.request();
    auto *ptr1 = (double *) buf1.ptr;
    auto *ptr2 = (double *) buf2.ptr;
    
    int element_count = 1;
    for (auto r: buf1.shape) {
      element_count *= r;
    }
    
    int nthreads;
    double prod = 0.;
    
    #pragma omp parallel
    {
        int id = omp_get_thread_num();
        if (id==0) nthreads = omp_get_num_threads();
    
        #pragma omp for reduction(+:prod)
        for (int i = 0; i < element_count; ++i) {
            prod += ptr1[i] * ptr2[i];
        }
    }
    
    py::print("I got", nthreads, "threads!");
    
    return prod;
}
        
PYBIND11_MODULE(dummy_work1, module) {
    module.def("parallel_pi", &parallel_pi);
    module.def("scalar_prod", &scalar_prod);
}


#include <omp.h>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>

namespace py = pybind11;

double calc_pi(int n) {
  /* Acquire GIL before calling Python code */
  py::gil_scoped_acquire acquire;

  int i;
  double step = 1.0/n;
  double s = 0;

  #pragma omp parallel
  {
    double x;
    #pragma omp for reduction(+:s)
    for (i=0; i<n; i++) {
      x = (i+0.5) * step;
      s += 4.0/(1 + x*x);
    }
  }
  return step * s;
};

PYBIND11_PLUGIN(openmpexample, m) {
    m.doc() = "";
  m.def("calc_pi", [](int n) {
      /* Release GIL before calling into C++ code */
      py::gil_scoped_release release;
      return calc_pi(n);
    });
}
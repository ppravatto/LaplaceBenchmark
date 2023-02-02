#include <pybind11/pybind11.h>
#include "laplace_benchmark.h"

namespace py = pybind11;

PYBIND11_MODULE(pylbench, module) {
    module.doc() = "Module to run CPUs benchmarks based on the Laplace grid algorithm";

    module.def(
        "laplace_2d_serial",
        &(laplace_benchmark::laplace_2d_serial),
        "Runs a serial laplace benchmark",
        py::arg("nrows"),
        py::arg("ncols"),
        py::arg("tol"),
        py::arg("max_iter"),
        py::arg("verbose"),
        py::return_value_policy::copy
    );

    module.def(
        "laplace_2d_openmp",
        &(laplace_benchmark::laplace_2d_openmp),
        "Runs a multi-threaded laplace benchmark using OpenMP",
        py::arg("nrows"),
        py::arg("ncols"),
        py::arg("tol"),
        py::arg("max_iter"),
        py::arg("verbose"),
        py::arg("nthreads"),
        py::return_value_policy::copy
    );
}



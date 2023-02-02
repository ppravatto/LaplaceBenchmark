#ifndef LAPLACE_BENCHMARK_H
#define LAPLACE_BENCHMARK_H

namespace laplace_benchmark{

    double laplace_2d_serial(int nrows, int ncols, double tol, int max_iter, bool verbose);
    double laplace_2d_openmp(int nrows, int ncols, double tol, int max_iter, bool verbose, int nthreads);
    
}

#endif
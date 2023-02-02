#include <iostream>
#include <omp.h>

#include "laplace_benchmark.h"

const int nrows = 1000;
const int ncols = 1000;
const int maxiter = 10000;
const double tol = 1e-8;

const int naverage = 10;

int main(){

    int maxthreads = omp_get_max_threads();

    for(int nthreads=1; nthreads<=maxthreads; nthreads++){

        std::cout << "Running benchmark on " << nthreads << " threads\n";

        double average_runtime = 0.;
        for (int i=0; i<naverage; i++){
            double runtime = laplace_benchmark::laplace_2d_openmp(nrows, ncols, tol, maxiter, false, nthreads);
            average_runtime += runtime;

            std::cout << "    -> Runtime " << i << ": " << runtime << " ms\n";
        }

        average_runtime /= double(naverage);
        std::cout << " -> Average runtime: " << average_runtime << " ms\n\n";

    }

    return 0;

}
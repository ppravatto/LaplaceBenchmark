#include <iostream>
#include <iomanip>
#include <algorithm>
#include <cmath>
#include <chrono>
#include <omp.h>

#include "laplace_benchmark.h"

namespace laplace_benchmark{

    double laplace_2d_serial(int nrows, int ncols, double tol, int max_iter, bool verbose){
        
        // Starting high resolution timer
        auto t1 = std::chrono::high_resolution_clock::now();
        
        //Allocaing memory for the matrices
        double * A = new double [nrows*ncols];
        double * A_new = new double [nrows*ncols];

        //Matrix initialization
        for(int r=0; r<nrows; r++){
            for(int c=0; c<ncols; c++){
                A[c+ncols*r] = (c==0)? 1. : 0.;
            }
        }

        // Iterate until either convergence is achieved or the maximum nuber of iterations is reached
        int iter = 0;
        double error = 1.;

        while(error > tol && iter < max_iter){
            
            error = 0.;

            // Compute average
            for(int r=1; r<nrows-1; r++){
                for(int c=1; c<ncols-1; c++){
                    A_new[c+ncols*r] = 0.25*(A[c+1+ncols*r] + A[c-1+ncols*r] + A[c+ncols*(r+1)] + A[c+ncols*(r-1)]);
                    error = std::max(error, std::abs(A_new[c+ncols*r] - A[c+ncols*r]));
                }
            }

            // Update A matrix
            for(int r=1; r<nrows-1; r++){
                for(int c=1; c<ncols-1; c++){
                    A[c+ncols*r] = A_new[c+ncols*r];
                }
            }

            // Print a message every 100 iterations if verbose is true 
            if(iter%100 == 0 && verbose == true){
                std::cout << std::scientific << std::setprecision(4) << "Iteration: " << iter << "\t" << "Error: " << error << std::endl;
            }

            iter++;
        }
        
        // Free the allocated memory
        delete[] A;
        delete[] A_new;

        // Stop the timer and compute runtime
        auto t2 = std::chrono::high_resolution_clock::now();
        double time = std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count();
        if(verbose == true) std::cout << "Runtime: " << time << "ms\n";
        
        return time;
    }

    double laplace_2d_openmp(int nrows, int ncols, double tol, int max_iter, bool verbose, int nthreads){

        // Set number of threads
        if(nthreads == -1){
            nthreads = omp_get_max_threads();
        }

        omp_set_num_threads(nthreads);
        
        // Starting high resolution timer
        auto t1 = std::chrono::high_resolution_clock::now();
        
        //Allocaing memory for the matrices
        double * A = new double [nrows*ncols];
        double * A_new = new double [nrows*ncols];

        //Matrix initialization
        #pragma omp parallel for
        for(int r=0; r<nrows; r++){
            for(int c=0; c<ncols; c++){
                A[c+ncols*r] = (c==0)? 1. : 0.;
            }
        }

        // Iterate until either convergence is achieved or the maximum nuber of iterations is reached
        int iter = 0;
        double error = 1.;

        while(error > tol && iter < max_iter){
            
            error = 0.;

            // Compute average
            #pragma omp parallel for reduction(max:error)
            for(int r=1; r<nrows-1; r++){
                for(int c=1; c<ncols-1; c++){
                    A_new[c+ncols*r] = 0.25*(A[c+1+ncols*r] + A[c-1+ncols*r] + A[c+ncols*(r+1)] + A[c+ncols*(r-1)]);
                    error = std::max(error, std::abs(A_new[c+ncols*r] - A[c+ncols*r]));
                }
            }

            // Update A matrix
            #pragma omp parallel for
            for(int r=1; r<nrows-1; r++){
                for(int c=1; c<ncols-1; c++){
                    A[c+ncols*r] = A_new[c+ncols*r];
                }
            }

            // Print a message every 100 iterations if verbose is true 
            if(iter%100 == 0 && verbose == true){
                std::cout << std::scientific << std::setprecision(4) << "Iteration: " << iter << "\t" << "Error: " << error << std::endl;
            }

            iter++;
        }
        
        // Free the allocated memory
        delete[] A;
        delete[] A_new;

        // Stop the timer and compute runtime
        auto t2 = std::chrono::high_resolution_clock::now();
        double time = std::chrono::duration_cast<std::chrono::milliseconds>(t2-t1).count();
        if(verbose == true) std::cout << "Runtime: " << time << "ms\n";
        
        return time;
    }

}
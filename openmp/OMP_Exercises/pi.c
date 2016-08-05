/*

This program will numerically compute the integral of

                  4/(1+x*x) 
				  
from 0 to 1.  The value of this integral is pi -- which 
is great since it gives us an easy way to check the answer.

The is the original sequential program.  It uses the timer
from the OpenMP runtime library

History: Written by Tim Mattson, 11/99.

*/
#include <stdio.h>
#include <omp.h>

static long num_steps = 1000000000;
double step;
int main ()
{
    int OMP_NUM_THREADS=atoi(getenv("OMP_NUM_THREADS"));
    printf("OMP_NUM_THREADS=%d\n", OMP_NUM_THREADS);
    int itr;
	double pi;
    double sum = 0.0;
	double start_time, run_time;
	step = 1.0/(double) num_steps;
    int load_size = num_steps / OMP_NUM_THREADS;
    
    printf("load_size=%d\n", load_size);
    
    double *psum = malloc(OMP_NUM_THREADS*8*8);
    memset(psum, 0, OMP_NUM_THREADS*8*8);

	start_time = omp_get_wtime();
    #pragma omp parallel
    {
        int i;
        int tid, start, end;
        double x;
        tid   = omp_get_thread_num();       
        start = tid * load_size;
        end   = (tid != OMP_NUM_THREADS-1) ? (tid+1)*load_size : num_steps;
            
        printf("tid=%d, start=%d, end=%d\n", tid, start, end);
        
	    for (i = start; i< end; i++)
        {
	        x = (i+0.5)*step;
	        *(psum+tid*8) = *(psum+tid*8) + 4.0/(1.0+x*x);
	    }

        #pragma omp atomic
        sum = sum + *(psum+tid*8);

        *(psum+tid*8) = 0;
    }
	run_time = omp_get_wtime() - start_time;

    for (itr=0; itr<OMP_NUM_THREADS; itr++)
    {
        printf("psum=%f\n", *(psum+itr*8));
        sum = sum + *(psum+itr*8);
    }
	pi = step * sum;

	printf("\n pi with %ld steps is %lf in %lf seconds\n ",num_steps,pi,run_time);
}	  

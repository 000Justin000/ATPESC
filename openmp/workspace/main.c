#include <stdio.h>
#include <omp.h>

void puttasks()
{
    for (int i=0; i<1000; i++)
        #pragma omp task
        printf("haha %d from %d\n", 1, omp_get_thread_num());
}

int main()
{
    #pragma omp parallel
    {
        #pragma omp master
        {
            for (int i=0; i<1000; i++)
                #pragma omp task
                printf("haha %d from %d\n", 0, omp_get_thread_num());
            
            #pragma omp task
            puttasks();

            for (int i=0; i<1000; i++)
                #pragma omp task
                printf("haha %d from %d\n", 2, omp_get_thread_num());
        }
        // there is a implicit #pragma omp taskwait, waiting all the tasks created in this region and their children to finish
        // this is a barrier?

        #pragma omp single
        {
            for (int i=0; i<1000; i++)
                #pragma omp task
                printf("haha %d from %d\n", 3, omp_get_thread_num());
        }            
    } 
    // the barrier in here wait until all the task in the parallel region to finish

    return 0;
}

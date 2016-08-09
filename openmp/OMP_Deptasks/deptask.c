#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include <unistd.h>

#define psize 10

int main()
{
    int acc[psize][psize];
    int mdf[psize][psize];

    for(int itr=0; itr<psize; itr++)
        for(int jtr=0; jtr<psize; jtr++)
        {
            acc[itr][jtr]=-1;
            mdf[itr][jtr]=-1;
        }

    #pragma omp parallel default(none) shared(acc,mdf)
    {
        #pragma omp single
        {
            for(int itr=0; itr<psize; itr++)
                for(int jtr=0; jtr<psize; jtr++)
                {
                    if(itr==0 && jtr==0)
                    {
                        #pragma omp task depend(out:acc[itr][jtr])
                        {
                            acc[itr][jtr] = 1;
                            mdf[itr][jtr] = omp_get_thread_num();
                            printf("(%d,%d) by %d\n", itr, jtr, omp_get_thread_num());
                            sleep(1);
                        }

                    }
                    else if(itr==0)
                    {
                        #pragma omp task depend(in:acc[itr][jtr-1]) depend(out:acc[itr][jtr])
                        {
                            acc[itr][jtr] = acc[itr][jtr-1];
                            mdf[itr][jtr] = omp_get_thread_num();
                            printf("(%d,%d) by %d\n", itr, jtr, omp_get_thread_num());
                            sleep(1);
                        }
                    }
                    else if(jtr==0)
                    {
                        #pragma omp task depend(in:acc[itr-1][jtr]) depend(out:acc[itr][jtr])
                        {
                            acc[itr][jtr] = acc[itr-1][jtr];
                            mdf[itr][jtr] = omp_get_thread_num();
                            printf("(%d,%d) by %d\n", itr, jtr, omp_get_thread_num());
                            sleep(1);
                        }
                    }
                    else if(itr>0 && jtr>0)
                    {
                        #pragma omp task depend(in:acc[itr][jtr-1],acc[itr-1][jtr]) depend(out:acc[itr][jtr])
                        {
                            acc[itr][jtr] = acc[itr][jtr-1]+acc[itr-1][jtr];
                            mdf[itr][jtr] = omp_get_thread_num();
                            printf("(%d,%d) by %d\n", itr, jtr, omp_get_thread_num());
                            sleep(1);
                        }
                    }
                }
        }
    }


    for(int itr=0; itr<psize; itr++)
    {
        for(int jtr=0; jtr<psize; jtr++)
            printf("%6d ", acc[itr][jtr]);
        printf("\n");
    }

    for(int itr=0; itr<psize; itr++)
    {
        for(int jtr=0; jtr<psize; jtr++)
            printf("%6d ", mdf[itr][jtr]);
        printf("\n");
    }

    return 0;
}

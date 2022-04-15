#include <stdlib.h>
#include <stdio.h>
#include <sys/time.h>

//const int NUM_PER_TASK = 262144;
//const int NUM_PER_TASK = 8192;
const int NUM_PER_TASK = 32768;
//const int NUM_PER_TASK = 256;

#pragma omp task in(array_data[0;NUM_PER_TASK]) out (result[0;1])
void array_sum (int array_data[NUM_PER_TASK], int result[1], int tid)
{
    int l_i;
    //printf("Execute a task %d\n", tid);
   
    result[0] = 0;
    for( l_i = 0; l_i < NUM_PER_TASK; l_i++ ) 
    {
        result[0] += array_data[l_i];
    }
}

int main (int argc, char* argv[])
{
    int **l_array_of_arrays;
    int *l_partial_sums;
    int l_num_tasks;
    int l_total;
    int l_i, l_j;
    struct timeval start;
    struct timeval stop;
    unsigned long elapsed;

    /******************************
     *      Check arguments       *
     ******************************/
    if(argc != 2)
    {
        printf("Usage: %s number_of_tasks\n", argv[0]);
        return 0;
    }
    l_num_tasks = atoi(argv[1]);

    /******************************
     *      Allocate array        *
     ******************************/
    l_partial_sums = (int*)malloc(l_num_tasks*sizeof(int));
    l_array_of_arrays = (int**)malloc(l_num_tasks*sizeof(int*));
    for(l_i = 0; l_i < l_num_tasks; l_i++)
    {
        l_array_of_arrays[l_i] = (int*)malloc(NUM_PER_TASK*sizeof(int));
        /******************************
        *      Initialize array      *
        ******************************/
       for(l_j = 0; l_j < NUM_PER_TASK; l_j++)
       {
           if((l_j % 2) == 0) l_array_of_arrays[l_i][l_j] = 1; //even
           else               l_array_of_arrays[l_i][l_j] = 0; //odd
       }
    }

    /******************************
     *         Call tasks         *
     ******************************/
    gettimeofday(&start,NULL);
    for(l_i = 0; l_i < l_num_tasks; l_i++)
    {
        array_sum(l_array_of_arrays[l_i], &(l_partial_sums[l_i]), l_i);
    }

    /******************************
     *    Accumulate results      *
     ******************************/
    l_total = 0;
    for(l_i = 0; l_i < l_num_tasks; l_i++)
    {
        #pragma omp taskwait on (l_partial_sums[l_i])
        //printf("%d -> %d\n", l_i, l_partial_sums[l_i]);
        l_total += l_partial_sums[l_i];
    }
    gettimeofday(&stop,NULL);

    printf("Result = %d\n",l_total);

    elapsed = 1000000 * (stop.tv_sec - start.tv_sec);
    elapsed += stop.tv_usec - start.tv_usec;

    printf("par_sec_time_us:%lu\n",elapsed);

    return 0;
}


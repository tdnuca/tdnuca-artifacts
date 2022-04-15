

#include <sys/time.h>
#include <malloc.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <assert.h>


/*
 * When the ratio of record movement across custers (led by centroids) hits
 * this low bound, the algorithms does not perform any more iterations.
 */
#define TERMINATION_THRESHOLD  0.01f

/*
 * Iteration hard limit.
 */
#define MAX_ITERATIONS         10

/*
 * Run the full number of iterations. We need this to get task timings.
 */
#define AVOID_THRESHOLD        0

/*
 * Size hint for the record block.
 */
#define BLOCK_SIZE             8192



/*********************************            *********************************/
/*********************************  SPU CODE  *********************************/
/*********************************            *********************************/


/*****************************************************************************/
/* SMP version: fake SPU ops. */
/*****************************************************************************/

typedef struct {
	float f1, f2, f3, f4;
} vector_float4;

#define spu_add(dst, a, b)			\
    do {					\
	dst.f1 = a.f1 + b.f1;			\
	dst.f2 = a.f2 + b.f2;			\
	dst.f3 = a.f3 + b.f3;			\
	dst.f4 = a.f4 + b.f4;			\
    } while(0)

#define spu_sub(dst, a, b)			\
    do {					\
	dst.f1 = a.f1 - b.f1;			\
	dst.f2 = a.f2 - b.f2;			\
	dst.f3 = a.f3 - b.f3;			\
	dst.f4 = a.f4 - b.f4;			\
    } while(0)

#define spu_mul(dst, a, b)			\
    do {					\
	dst.f1 = a.f1 * b.f1;			\
	dst.f2 = a.f2 * b.f2;			\
	dst.f3 = a.f3 * b.f3;			\
	dst.f4 = a.f4 * b.f4;			\
    } while(0)

#define spu_extract(a, idx)			\
    ( (idx == 0) ? a.f1				\
      : ( (idx == 1) ? a.f2			\
	  : ( (idx == 2) ? a.f3			\
	      : a.f4)))

/*****************************************************************************/
/*****************************************************************************/

/******************************************************************************
 *                                                                            *
 * Find the nearest centroid for a record.  Returns the center identifier     *
 * (index).                                                                   *
 *                                                                            *
 * This is a straightforward adaptation of the original function              *
 *                                                                            *
 ******************************************************************************/
int compare_to_centers (float  *record,
                        int     record_length,
                        int     number_of_centers,
                        float  *centers)
{
  short k;
  short i,j;
  short min_id=0;
  float min = 999999999;
  short blocksize=12;
  short loopcount = record_length / blocksize;
  short extra = record_length % blocksize;
  float *data2=centers;//+k*(cb.record_length);
  float *data1=record;

  float distance;
  //float distance_3 = 0.0;
  //float distance_4 = 0.0;
  //float distance_5= 0.0;
  int offsetX = record_length;
  float *data3=centers+offsetX;//+k*(record_length);
  float Cdistance;
  //float Cdistance_3 = 0.0;
  //float Cdistance_4 = 0.0;
  //float Cdistance_5= 0.0;

  offsetX +=record_length;
  float *data4=centers+offsetX;//+k*(cb.record_length);
  float Ddistance;
  //float Ddistance_3 = 0.0;
  //float Ddistance_4 = 0.0;
  //float Ddistance_5= 0.0;

  offsetX +=record_length;
  float *data5=centers+offsetX;//+k*(cb.record_length);
  float Edistance;
  //float Edistance_3 = 0.0;
  //float Edistance_4 = 0.0;
  //float Edistance_5= 0.0;

  for(k=0; k< number_of_centers-4; k+=4) { //for each center
	  //distance_3 = 0.0;
	  //distance_4 = 0.0;
	  //distance_5= 0.0;
	  vector_float4* varray1 = (vector_float4*)(data1);
	  vector_float4* varray2 = (vector_float4*)(data2);
	  vector_float4 result;
	  vector_float4 result2;
	  vector_float4 result3;
	  
	  //Cdistance_3 = 0.0;
	  //Cdistance_4 = 0.0;
	  //Cdistance_5= 0.0;
	  vector_float4* Cvarray2 = (vector_float4*)(data3);
	  vector_float4 Cresult;
	  vector_float4 Cresult2;
	  vector_float4 Cresult3;
	  
	  //Ddistance_3 = 0.0;
	  //Ddistance_4 = 0.0;
	  //Ddistance_5= 0.0;
	  vector_float4* Dvarray2 = (vector_float4*)(data4);
	  vector_float4 Dresult;
	  vector_float4 Dresult2;
	  vector_float4 Dresult3;
	  
	  
	  //Edistance_3 = 0.0;
	  //Edistance_4 = 0.0;
	  //Edistance_5= 0.0;
	  vector_float4* Evarray2 = (vector_float4*)(data5);
	  vector_float4 Eresult;
	  vector_float4 Eresult2;
	  vector_float4 Eresult3;
	  vector_float4 temp,Ctemp,Dtemp,Etemp;
	  int index=0;
	  
	  //accumulators
	  vector_float4 C1acc = {0,0,0,0};
	  vector_float4 C2acc = {0,0,0,0};
	  vector_float4 C3acc = {0,0,0,0};
	  vector_float4 C4acc = {0,0,0,0};
	  for(i=0;i<loopcount;i++) {
		  //prof_cp1();
		  spu_sub(result, varray1[index],varray2[index]);
		  spu_sub(result2, varray1[index+1],varray2[index+1]);
		  spu_sub(result3, varray1[index+2],varray2[index+2]);
		  spu_sub(Cresult, varray1[index],Cvarray2[index]);
		  spu_sub(Cresult2, varray1[index+1],Cvarray2[index+1]);
		  spu_sub(Cresult3, varray1[index+2],Cvarray2[index+2]);
		  spu_sub(Dresult, varray1[index],Dvarray2[index]);
		  spu_sub(Dresult2, varray1[index+1],Dvarray2[index+1]);
		  spu_sub(Dresult3, varray1[index+2],Dvarray2[index+2]);
		  spu_sub(Eresult, varray1[index],Evarray2[index]);
		  spu_sub(Eresult2, varray1[index+1],Evarray2[index+1]);
		  spu_sub(Eresult3, varray1[index+2],Evarray2[index+2]);
		  
		  spu_mul(result, result,result);
		  spu_mul(result2, result2,result2);
		  spu_mul(result3, result3,result3);
		  spu_mul(Cresult, Cresult,Cresult);
		  spu_mul(Cresult2, Cresult2,Cresult2);
		  spu_mul(Cresult3, Cresult3,Cresult3);
		  spu_mul(Dresult, Dresult,Dresult);
		  spu_mul(Dresult2, Dresult2,Dresult2);
		  spu_mul(Dresult3, Dresult3,Dresult3);
		  spu_mul(Eresult, Eresult,Eresult);
		  spu_mul(Eresult2, Eresult2,Eresult2);
		  spu_mul(Eresult3, Eresult3,Eresult3);
		  
		  spu_add(C1acc, result,result2);
		  spu_add(C2acc, Cresult,Cresult2);
		  spu_add(C3acc, Dresult,Dresult2);
		  spu_add(C4acc, Eresult,Eresult2);
		  spu_add(C1acc, C1acc,result3);
		  spu_add(C2acc, C2acc,Cresult3);
		  spu_add(C3acc, C3acc,Dresult3);
		  spu_add(C4acc, C4acc,Eresult3);
		  
		  index = index+3;
	  }
	  distance=spu_extract(C1acc,0)+spu_extract(C1acc,1)+spu_extract(C1acc,2)+spu_extract(C1acc,3);
	  Cdistance=spu_extract(C2acc,0)+spu_extract(C2acc,1)+spu_extract(C2acc,2)+spu_extract(C2acc,3);
	  Ddistance=spu_extract(C3acc,0)+spu_extract(C3acc,1)+spu_extract(C3acc,2)+spu_extract(C3acc,3);
	  Edistance=spu_extract(C4acc,0)+spu_extract(C4acc,1)+spu_extract(C4acc,2)+spu_extract(C4acc,3);
	  
	  
	  if(__builtin_expect((extra>0),0)) {
		  for(i=0;i<extra;i++) {
			  float distance1 =data1[index+i]-data2[index+i];
			  distance += distance1*distance1;
			  float Cdistance1 =data1[index+i]-data3[index+i];
			  Cdistance += Cdistance1*Cdistance1;
			  float Ddistance1 =data1[index+i]-data4[index+i];
			  Ddistance += Ddistance1*Ddistance1;
			  float Edistance1 =data1[index+i]-data5[index+i];
			  Edistance += Edistance1*Edistance1;
		  }
	  }
	  //vector float min1 =spu_promote(min,0);
	  //vector float dis1 =spu_promote(distance,0);
	  //vector unsigned int sel = spu_cmpgt(min1,dis1);
	  //min = spu_extract(spu_sel(dis1,min1,sel),0);
	  //min_id = spu_extract(spu_sel(k,min_id,sel),0);
	  //min = spu_extract(select);
	  
	  /* this is made spu_sel by the compiler */
	  if(__builtin_expect((min > distance),0)) {
		  min_id=k;
		  min=distance;
	  }
	  if(__builtin_expect((min > Cdistance),0)) {
		  min_id=k+1;
		  min=Cdistance;
	  }
	  if(__builtin_expect((min > Ddistance),0)) {
		  min_id=k+2;
		  min=Ddistance;
	  }
	  if(__builtin_expect((min > Edistance),0)) {
		  min_id=k+3;
		  min=Edistance;
	  }
	  data2+=record_length*4;
	  data3+=record_length*4;
	  data4+=record_length*4;
	  data5+=record_length*4;
  } //for each center mod 4
  
  return min_id;
}


/******************************************************************************
 *                                                                            *
 * Collect new centroid data from current centroids and the given chunk of    *
 * records.  The output data will be combined toghether from the other tasks  *
 * at the PPU.                                                                *
 *                                                                            *
 ******************************************************************************/
#pragma omp task in([number_of_records * DIMENSIONS]records, [CENTERS * DIMENSIONS]centers) \
                 inout([number_of_records]assigned_center) \
                 inout([CENTERS * DIMENSIONS]newcenters, [CENTERS + 1]newcenters_histogram)
void kmeans_calculate (int     DIMENSIONS,
                       int     CENTERS,
                       int     number_of_records,
                       float  *records,
                       float  *centers,
                       int    *assigned_center,
                       float  *newcenters,
                       int    *newcenters_histogram)
{
        int  i, j;
        int  min;

        for (i = 0;  i < number_of_records;  i++)
        {
                /* Find this record's nearest centroid */
                min     = compare_to_centers(&records[i * DIMENSIONS],
                                             DIMENSIONS, CENTERS, centers);
                if (assigned_center[i] != min)
                {
                        newcenters_histogram[CENTERS]++;
                        assigned_center[i] = min;
                }
                /* Update centroid's cluster information */
                newcenters_histogram[min]++;
                for (j = 0;  j < DIMENSIONS;  j++)
                        newcenters[min * DIMENSIONS + j] += records[i * DIMENSIONS + j];
        }
}



/*********************************            *********************************/
/*********************************  PPU CODE  *********************************/
/*********************************            *********************************/


/******************************************************************************
 *                                                                            *
 * Make a timer snapshot and, optionally, print the elapsed time since the    *
 * last time the function was called.                                         *
 *                                                                            *
 ******************************************************************************/
unsigned long time_int (int print)
{

        unsigned long elapsed_usec;
        static struct timeval  t1;  /* var for previous time stamp */
        static struct timeval  t2;  /* var of current time stamp */

        if (gettimeofday(&t2, NULL) == -1)
        {
                perror("gettimeofday");
                exit(9);
        }

	elapsed_usec = (t2.tv_sec - t1.tv_sec) * 1e6 +
		(t2.tv_usec - t1.tv_usec);
	
        if (print)
        {
                printf("Time spent [%.2fs] \n", ((double)elapsed_usec)*1e-6);
        }

        t1 = t2;
        return elapsed_usec;
}

/******************************************************************************
 * Element type in generic functions. 
 ******************************************************************************/

#define ELEM_INT    1
#define ELEM_FLOAT  2

/******************************************************************************
 *                                                                            *
 * Export memset into a task.
 *                                                                            *
 ******************************************************************************/

#pragma omp task out([bytes]buf1, [bytes]buf2, [bytes]buf3, [bytes]buf4) 
void memzero4_task(char* buf1, char* buf2, char* buf3, char* buf4, size_t bytes)
{
	int i;

	int size = bytes / sizeof(int);
	int* ibuf1 = (int*)buf1;
	int* ibuf2 = (int*)buf2;
	int* ibuf3 = (int*)buf3;
	int* ibuf4 = (int*)buf4;

	for(i=0; i<size; i++) {
		ibuf1[i] = ibuf2[i] = ibuf3[i] = ibuf4[i] = 0;
	}
}


/******************************************************************************
 *                                                                            *
 * Combine results produced by the tasks and calculate new centroids for the  *
 * next clustering iteration.                                                 *
 *                                                                            *
 ******************************************************************************/

#pragma omp task in([bytes]buf1, [bytes]buf2, [bytes]buf3)  \
	inout([bytes]buf0)
void sum4_task(int etype,
               int bytes,
               char* buf0,
               char* buf1,
               char* buf2,
               char* buf3)
{
	if(etype == ELEM_INT) {
		int size = bytes / sizeof(int);
		int* ibuf0 = (int*)buf0;
		int* ibuf1 = (int*)buf1;
		int* ibuf2 = (int*)buf2;
		int* ibuf3 = (int*)buf3;
		for(int i=0; i<size; i++) {
			ibuf0[i] += ibuf1[i] + ibuf2[i] + ibuf3[i];
		}
	} else if(etype == ELEM_FLOAT) {
		int size = bytes / sizeof(float);
		float* fbuf0 = (float*)buf0;
		float* fbuf1 = (float*)buf1;
		float* fbuf2 = (float*)buf2;
		float* fbuf3 = (float*)buf3;
		for(int i=0; i<size; i++) {
			fbuf0[i] += fbuf1[i] + fbuf2[i] + fbuf3[i];
		}
	}
}

void sum2Darray_i(int rows, int cols, int** array)
{
	// tree sum: first sum each group of 4 adjacent arrays, then sum groups
	// of arrays with distance 4, then distance 16, etc.
	// indices:
	// d: stride
	// base: first row in current sum
	for(int d=1; d <= (rows/4) ; d *= 4) {
		for(int base=0; base < rows; base += 4*d) {
			sum4_task(ELEM_INT, cols*sizeof(int), 
				  (char*)array[base], 
				  (char*)array[base + 1*d],
				  (char*)array[base + 2*d],
				  (char*)array[base + 3*d]);
		}
	}
}

void sum2Darray_f(int rows, int cols, float** array)
{
	// tree sum: first sum each group of 4 adjacent arrays, then sum groups
	// of arrays with distance 4, then distance 16, etc.
	// indices:
	// d: stride
	// base: first row in current sum
	for(int d=1; d <= (rows/4) ; d *= 4) {
		for(int base=0; base < rows; base += 4*d) {
			sum4_task(ELEM_FLOAT, cols*sizeof(float), 
				  (char*)array[base], 
				  (char*)array[base + 1*d],
				  (char*)array[base + 2*d],
				  (char*)array[base + 3*d]);
		}
	}
}


#pragma omp task in([CENTERS+1]newcenters_histograms0) \
	in([CENTERS*DIMENSIONS]newcenters0)	 \
	out([CENTERS*DIMENSIONS]centers)
void recalc_loop2(int DIMENSIONS,
		  int CENTERS,
		  int* newcenters_histograms0, 
		  float* newcenters0, 
		  float* centers)
{
	int i, j;

        /* Re-average centroids directly on the centroids array */
        for (i = 0;  i < CENTERS;  i++) {
                for (j = 0;  j < DIMENSIONS;  j++)
                {
                        if (newcenters_histograms0[i] > 0)
                                centers[i*DIMENSIONS+j] = newcenters0[i*DIMENSIONS+j] / newcenters_histograms0[i];
                        else
                                centers[i*DIMENSIONS+j] = -1;
                }
	}

}

void recalculate_centers (int      DIMENSIONS,
                          int      CENTERS,
                          int      num_spus,
                          float   *centers,
                          float  **newcenters,
                          int    **newcenters_histograms)
{
        int  i, j, k;

	
	sum2Darray_i(num_spus, CENTERS, newcenters_histograms);
	sum2Darray_f(num_spus, CENTERS*DIMENSIONS, newcenters);
#if 0
        /* Combine task results (co-ordinates and histograms) */
        for (k = 1;  k < num_spus;  k++) {
                for (i = 0;  i < CENTERS;  i++)
                {
                        newcenters_histograms[0][i] += newcenters_histograms[k][i];
                        for (j = 0;  j < DIMENSIONS;  j++)
                                newcenters[0][i*DIMENSIONS+j] += newcenters[k][i*DIMENSIONS+j];
                }
	}
#endif

	// export the loop into a task, which acts as a barrier because of
	// 'centers'.
	recalc_loop2(DIMENSIONS, CENTERS, 
		     newcenters_histograms[0], newcenters[0], centers);
}


/******************************************************************************
 *                                                                            *
 * Main function.                                                             *
 *                                                                            *
 ******************************************************************************/
int main (int argc, char **argv)
{
        /* Algorithm parameters */
        int      NUMBER_OF_RECORDS, DIMENSIONS, CENTERS, CSS_NUM_SPUS;
        int      number_of_records, block_size, block_records;
        int      localstore_usage;
        /* Auxiliary variables */
        int      k, i, j, z;
        int      iteration, changed;
        /* Other computations */
        char    *env_num_spus;
        int      sizeofdistance;
        double   gflops;
        unsigned long elapsed_usec;
        /* Data arrays */
        float   *records;
        float   *centers;
        int     *assigned_centers;
        float  **newcenters;
        int    **newcenters_histograms;
	int newcenters_elms;
	int newcenters_histograms_elms;

	srand(17);

        if (argc < 4)
        {
                printf("%s <data_points> <dim> <num_centers>\n", argv[0]);
                printf("eg: %s 200000 60 64\n", argv[0]);
                return 1;
        }

        NUMBER_OF_RECORDS  = atoi(argv[1]);
        DIMENSIONS         = atoi(argv[2]);
        CENTERS            = atoi(argv[3]);

        /* Check number of centers (to avoid uchar overflow) */
        if (CENTERS > 256)
        {
                printf("Current program only supports 256 centers\n.");
                return 2;
        }

        env_num_spus = getenv("CSS_NUM_SPUS");
        if (env_num_spus != NULL)
                CSS_NUM_SPUS = atoi(env_num_spus);
        else
                CSS_NUM_SPUS = 8;
	CSS_NUM_SPUS = 64;


	// When summing the 2D arrays, we sum 4 rows at a time, so we need to
	// make sure the number of rows is a power of 4.
	assert(((CSS_NUM_SPUS - 1) & CSS_NUM_SPUS) == 0); // power of 2
	// TODO: make sure the power of 2 is even


        block_size    = DIMENSIONS * sizeof(float);
        block_records = 1;
        while (block_size < BLOCK_SIZE)
        {
                block_size    *= 2;
                block_records *= 2;
        }

        localstore_usage = block_size + block_records * sizeof(int) +
                           2 * (CENTERS * DIMENSIONS * sizeof(float)) +
                           (CENTERS + 1) * sizeof(int);

        /* Print precalculated information */
        printf("***********************************************\n");
        printf("***********************************************\n");
        printf("******                                   ******\n");
        printf("******          CellSs K-Means           ******\n");
        printf("******  BARCELONA SUPERCOMPUTING CENTER  ******\n");
        printf("******                                   ******\n");
        printf("******     ORIGINAL CBE VERSION FROM     ******\n");
        printf("******       OHIO STATE UNIVERSITY       ******\n");
        printf("******                                   ******\n");
        printf("***********************************************\n");
        printf("***********************************************\n\n");
        printf("Number of data points is %d\n",     NUMBER_OF_RECORDS);
        printf("Number of centers is %d\n",         CENTERS);
        printf("Number of dimensions %d\n\n",       DIMENSIONS);

        printf("Record array chunk size is %d\n",     block_size);
        printf("Number of records per chunk is %d\n", block_records);
        printf("Predicted localstore usage is %d bytes\n\n", localstore_usage);

	printf("kmeans_calculate parallelism: %d\n", 
	       NUMBER_OF_RECORDS/block_records);

        /* Allocate space for data */
        records               = (float*)memalign(128, NUMBER_OF_RECORDS * DIMENSIONS * sizeof(float));
        centers               = (float*)memalign(128, CENTERS * DIMENSIONS * sizeof(float));
        assigned_centers      = (int*)memalign(128, NUMBER_OF_RECORDS * sizeof(int));
        newcenters            = (float**)malloc(CSS_NUM_SPUS * sizeof(float *));
	newcenters_elms       = CENTERS * DIMENSIONS;
        newcenters_histograms = (int**)malloc(CSS_NUM_SPUS * sizeof(int *));
	newcenters_histograms_elms = (CENTERS + 1);
        for (k = 0;  k < CSS_NUM_SPUS;  k++)
        {
                newcenters[k]            = 
			(float*)memalign(128, 
					 newcenters_elms * sizeof(float));
                newcenters_histograms[k] = 
			(int*)memalign(128, 
				       newcenters_histograms_elms*sizeof(int));
        }

        /* time_int(0); */

        /* Make record data */
        for (i = 0;  i < NUMBER_OF_RECORDS;  i++) {
                for (j = 0;  j < DIMENSIONS;  j++) {
                        records[i * DIMENSIONS + j] = rand() / 100000.0f + j;
		}
	}

        /* Initialize the centers */
        for (i = 0;  i < CENTERS;  i++) {
                memcpy(&centers[i * DIMENSIONS], &records[i * DIMENSIONS], 
		       DIMENSIONS * sizeof(float));
	}

        /* Clear assigned centers */
        memset(assigned_centers, -1, NUMBER_OF_RECORDS * sizeof(int));

        /* time_int(1); */

        changed        = 0;
        iteration      = 0;
        gflops         = 0.0;
        sizeofdistance = DIMENSIONS * 3 + 1;

	fprintf(stderr, "Start parallel section\n");

        time_int(0);

        /* Main loop */
        do {

                for (k = 0;  k < CSS_NUM_SPUS;  k+=4)
                {
                        memzero4_task((char*)newcenters[k],
				      (char*)newcenters[k+1],
				      (char*)newcenters[k+2],
				      (char*)newcenters[k+3],
				      newcenters_elms * sizeof(float));

                        memzero4_task((char*)newcenters_histograms[k], 
				      (char*)newcenters_histograms[k+1], 
				      (char*)newcenters_histograms[k+2], 
				      (char*)newcenters_histograms[k+3], 
				      newcenters_histograms_elms*sizeof(int));
                }

                k = 0;
		z = 0;
                for (i = 0;  i < NUMBER_OF_RECORDS;  i += block_records)
                {
                        number_of_records = (i + block_records > NUMBER_OF_RECORDS ? NUMBER_OF_RECORDS - i : block_records);

                        kmeans_calculate(DIMENSIONS, CENTERS, number_of_records,
                                         &records[i * DIMENSIONS], centers,
                                         &assigned_centers[i], newcenters[k],
                                         newcenters_histograms[k]);

                        k = (k + 1) % CSS_NUM_SPUS;
			z++;
                }

                /* Delay accounting calculations in order to launch tasks as
                 * soon as possible */
                //if (changed)
                //        printf("TOTAL CHANGED IS %d\n\n", changed);
                gflops += (double) NUMBER_OF_RECORDS * CENTERS * sizeofdistance;
                iteration++;
                printf("EXECUTING ITERATION #%d, calculations %d\n", iteration,z);

#pragma omp taskwait
                changed = 0;
                for (k = 0;  k < CSS_NUM_SPUS;  k++)
                        changed += newcenters_histograms[k][CENTERS];

		fprintf(stderr, "change: %f, threshold: %f\n",
			(float)changed/NUMBER_OF_RECORDS, (float)TERMINATION_THRESHOLD);

                if (!(AVOID_THRESHOLD) && 
		    ((float)changed/NUMBER_OF_RECORDS < TERMINATION_THRESHOLD))
                {
                        //printf("TOTAL CHANGED IS %d\n\n", changed);
                        //printf("Number of changed records has fallen below the threshold of %f\n\n", TERMINATION_THRESHOLD);
                        break;
                }
                if (iteration >= MAX_ITERATIONS)
                {
                        //printf("TOTAL CHANGED IS %d\n\n", changed);
                        //printf("We are at the maximum number of iterations ...\n\n");
                        break;
                }

                recalculate_centers(DIMENSIONS, CENTERS, CSS_NUM_SPUS,
                                    centers, newcenters, newcenters_histograms);
        } while (1);
#pragma omp taskwait
        elapsed_usec = time_int(1);
	fprintf(stderr, "Finished parallel section\n");
	
        gflops = gflops / (1024.0 * 1024.0 * 1024.0);
        printf(" ROUGH ESTIMATED GFLOPS=%lf\n\n\n", gflops / (((double)elapsed_usec) * 1e-6) );

        printf("par_sec_time_us:%lu\n",elapsed_usec);
        return 0;
}


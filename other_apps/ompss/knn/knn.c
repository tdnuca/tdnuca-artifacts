#include <sys/time.h>
#include <malloc.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

/*
 * Size threshold for array chunks.  Training and test arrays will not fit in
 * the localsotre, training information but it generally won't fit in the local
 * store; therefore, they will have to be splitted in smaller pieces.  This
 * macro controls the size threshold of the array blocks.
 *
 * Use only power of two values.
 */
#define BLOCK_SIZE  32768
//#define BLOCK_SIZE  8192 /* SMP version */



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
		dst.f1 = a.f1 + b.f1;		\
		dst.f2 = a.f2 + b.f2;		\
		dst.f3 = a.f3 + b.f3;		\
		dst.f4 = a.f4 + b.f4;		\
	} while(0)

#define spu_sub(dst, a, b)			\
	do {					\
		dst.f1 = a.f1 - b.f1;		\
		dst.f2 = a.f2 - b.f2;		\
		dst.f3 = a.f3 - b.f3;		\
		dst.f4 = a.f4 - b.f4;		\
	} while(0)

#define spu_mul(dst, a, b)			\
	do {					\
		dst.f1 = a.f1 * b.f1;		\
		dst.f2 = a.f2 * b.f2;		\
		dst.f3 = a.f3 * b.f3;		\
		dst.f4 = a.f4 * b.f4;		\
	} while(0)

#define spu_extract(a, idx)			\
	( (idx == 0) ? a.f1			\
	  : ( (idx == 1) ? a.f2			\
	      : ( (idx == 2) ? a.f3		\
		  : a.f4)))

/******************************************************************************
 *                                                                            *
 * Find the classes of the K nearest neighbours.  Trivial adaptation from the *
 * original one.  Which basically means no readability improvements were      *
 * carried out.                                                               *
 *                                                                            *
 ******************************************************************************/
#define RECORD_SIZE  DIMENSIONS
#pragma omp task in([train_size * (DIMENSIONS + 1)]Train)   \
	inout([data_size * (DIMENSIONS + K)]Data)	\
	inout([data_size * K]local_distance_data)	\
	inout([data_size]ASSIGNED_CLASSES)
void compare (int     K,
              int     DIMENSIONS,
              int     CLASSES,
              int     train_size,
              int     data_size,
              float  *Train,
              float  *Data,
              float  *local_distance_data,
              int    *ASSIGNED_CLASSES)
{
	//  prof_cp30();
	//prof_cp0();
	//loop data
	short loopsize = RECORD_SIZE / 12;
	short extra = RECORD_SIZE % 12;
	float *record = Data;
	float *train = Train;
	float distance=0;
	float Xdistance=0;
	short i,j,k;
	record=Data;
	vector_float4 result;
	vector_float4 Bresult;
	vector_float4 Cresult;
	vector_float4 Xresult;
	vector_float4 XBresult;
	vector_float4 XCresult;
	vector_float4* varray1;
	vector_float4* varray2;
	vector_float4* Xvarray2;
	float Bdis;
	float Cdis;
	float XBdis;
	short len = RECORD_SIZE/4;
	float XCdis;
	short TRAIN_STEP=2*(RECORD_SIZE+1);
	float record_max;
	//printf("train size is %d, data size is %d\n",train_size,data_size);
	for(i=0;i<data_size;i++) //for each data point
	{
		//printf("looking at record %d\n",i);
		varray1 = (vector_float4*)record;
		//reset training data
		train=Train;
		for(j=0;j<train_size;j+=2) //for each training point
		{
			distance=0;
			Xdistance=0;
			varray2 = (vector_float4*)train;
			Xvarray2 = (vector_float4*)(train+RECORD_SIZE+1);
			for(k=0;k<len;k+=3) {
				//compare data to first training point
				spu_sub(result, varray1[k],varray2[k]);
				spu_sub(Bresult, varray1[k+1],varray2[k+1]);
				spu_sub(Cresult, varray1[k+2],varray2[k+2]);
				spu_sub(Xresult, varray1[k],Xvarray2[k]);
				spu_sub(XBresult, varray1[k+1],Xvarray2[k+1]);
				spu_sub(XCresult, varray1[k+2],Xvarray2[k+2]);
				spu_mul(result, result,result);
				spu_mul(Bresult, Bresult,Bresult);
				spu_mul(Cresult, Cresult,Cresult);
				spu_mul(Xresult, Xresult,Xresult);
				spu_mul(XBresult, XBresult,XBresult);
				spu_mul(XCresult, XCresult,XCresult);
				distance += spu_extract(result,0)+spu_extract(result,1)+spu_extract(result,2)+spu_extract(result,3);
				Bdis = spu_extract(Bresult,0)+spu_extract(Bresult,1)+spu_extract(Bresult,2)+spu_extract(Bresult,3);
				Cdis = spu_extract(Cresult,0)+spu_extract(Cresult,1)+spu_extract(Cresult,2)+spu_extract(Cresult,3);
				Xdistance += spu_extract(Xresult,0)+spu_extract(Xresult,1)+spu_extract(Xresult,2)+spu_extract(Xresult,3);
				XBdis = spu_extract(XBresult,0)+spu_extract(XBresult,1)+spu_extract(XBresult,2)+spu_extract(XBresult,3);
				XCdis = spu_extract(XCresult,0)+spu_extract(XCresult,1)+spu_extract(XCresult,2)+spu_extract(XCresult,3);
				distance +=Bdis+Cdis;
				Xdistance +=XBdis+XCdis;

			}
			if(__builtin_expect((extra > 0),0))
			{
				//printf("here\n");
				short w;
				for(w=0;w<extra; w++)
				{
					float dis = record[12*loopsize+w]-train[12*loopsize+w];
					distance +=dis*dis;
					float Xdis = record[12*loopsize+w]-train[12*loopsize+RECORD_SIZE+1+w];
					Xdistance +=Xdis*Xdis;
				}
			}
	  
			short rr=0;
			float temp2;
			float temp3;
			int temp4;
			//if we havent already added K items
			if(__builtin_expect(ASSIGNED_CLASSES[i]<K,0))
			{
				//printf("assigned is %d, distance is %.2f, data point is %d\n",ASSIGNED_CLASSES[i],distance,i);
				rr=0;
				//print_clases(record);
				while(rr<K)
				{
					if(local_distance_data[i*K+rr] > distance)
					{
						rr++;
					}
					else
					{
						temp4=(int)train[RECORD_SIZE];
						//printf("label is %d\n",temp4);
						while(rr<K)
						{
							temp2 = local_distance_data[i*K+rr];
							local_distance_data[i*K+rr]=distance;
							distance = temp2;
							//set class too
							temp3 = record[RECORD_SIZE+rr];
							record[RECORD_SIZE+rr]=temp4;
							temp4=temp3;
							rr++;
						}
					}
				}
				//print_classes(record);
				//getchar();
				ASSIGNED_CLASSES[i]++;
			}
			else
			{
				if(__builtin_expect(local_distance_data[i*K] > distance,0))//its smaller
				{
					rr=K-1;
					while(rr>=0)
					{
						if(local_distance_data[i*K+rr] < distance)
						{
							rr--;
						}
						else
						{
							temp4=(int)train[RECORD_SIZE];
							//printf("label is %d\n",temp4);
							while(rr>=0)
							{
								temp2 = local_distance_data[i*K+rr];
								local_distance_data[i*K+rr]=distance;
								distance = temp2;
								//set class too
								temp3 = record[RECORD_SIZE+rr];
								record[RECORD_SIZE+rr]=temp4;
								temp4=temp3;
								rr--;
							}
						}
					}
				}
			}
			//now do the second point
			if(__builtin_expect(ASSIGNED_CLASSES[i]<K,0))
			{
				//printf("assigned is %d, distance is %.2f, data point is %d\n",ASSIGNED_CLASSES[i],distance,i);
				rr=0;
				//print_classes(record);
				while(rr<K)
				{
					if(local_distance_data[i*K+rr] > Xdistance)
					{
						rr++;
					}
					else
					{
						//get the class label
						temp4=(int)train[2*RECORD_SIZE+1];
						//printf("label is %d\n",temp4);
				  
						while(rr<K)
						{
							temp2 = local_distance_data[i*K+rr];
							local_distance_data[i*K+rr]=Xdistance;
							Xdistance = temp2;
							//set class too
							temp3 = record[RECORD_SIZE+rr];
							record[RECORD_SIZE+rr]=temp4;
							temp4=temp3;
							rr++;
						}
					}
				}
				//print_classes(record);
				//getchar();
				ASSIGNED_CLASSES[i]++;
			}
			else
			{
				if(__builtin_expect(local_distance_data[i*K] > Xdistance,0))//its smaller
				{
					rr=K-1;
					while(rr>=0)
					{
						if(local_distance_data[i*K+rr] < Xdistance)
						{
							rr--;
						}
						else
						{
							temp4=(int)train[2*RECORD_SIZE+1];
							//printf("label is %d\n",temp4);
					  
							while(rr>=0)
							{
								temp2 = local_distance_data[i*K+rr];
								local_distance_data[i*K+rr]=Xdistance;
								Xdistance = temp2;
								//set class too
								temp3 = record[RECORD_SIZE+rr];
								record[RECORD_SIZE+rr]=temp4;
								temp4=temp3;
								rr--;
							}
						}
					}
				}
			}
			train+=TRAIN_STEP;
			//print_classes(record);
		}//for each training point
		//we can now put the proper
		//print_classes(record);
		//getchar();
		record += DIMENSIONS + K;//full length, including k assignments
	}
	//prof_cp1();
}
#undef RECORD_SIZE


/******************************************************************************
 *                                                                            *
 * Assign the class of each test point as the most frequent among the K       *
 * nearest neighbours.                                                        *
 *                                                                            *
 ******************************************************************************/
#pragma omp task in([records * (K + DIMENSIONS)]data)	\
	out([records]class_labels)
void get_best_classlabels (int     K,
                           int     DIMENSIONS,
                           int     CLASSES,
                           int     records,
                           float  *data,
                           int    *class_labels)
{
        int  i, j;
        int  label, max_class_id, max_classcount;
        int  classcounts[CLASSES];

        for (i = 0;  i < records;  i++)
        {
                memset(classcounts, 0, sizeof(int) * CLASSES);

                /* Count nearest neighbour classes */
                for (j = DIMENSIONS;  j < DIMENSIONS + K;  j++)
                {
                        label = data[i * (DIMENSIONS + K) + j];
                        classcounts[label]++;
                }

                /* Select the most popular one */
                max_class_id   = -1;
                max_classcount =  0;
                for (j = 0;  j < CLASSES;  j++)
                        if (classcounts[j] > max_classcount)
                        {
                                max_class_id   = j;
                                max_classcount = classcounts[j];
                        }
                class_labels[i] = max_class_id;
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
        unsigned long elapsed_usec = 0;
        static struct timeval  t1;  /* previous time stamp */
        static struct timeval  t2;  /* current time stamp */

        if (gettimeofday(&t2, NULL) == -1)
        {
                perror("gettimeofday");
                exit(9);
        }

        if (print)
        {
                elapsed_usec = (t2.tv_sec  - t1.tv_sec) * 1e6 +
			(t2.tv_usec - t1.tv_usec);
                printf("Time spent [%.2fs] \n", ((double)elapsed_usec)*1e-6);
        }

        t1 = t2;
        return elapsed_usec;
}


/******************************************************************************
 *                                                                            *
 * Main function.                                                             *
 *                                                                            *
 ******************************************************************************/
int main (int argc, char **argv)
{
        /* Algorithm parameters */
        int     K, DIMENSIONS, CLASSES, TRAINING_RECORDS, TEST_RECORDS;
        int     test_block_records, train_block_records;
        int     number_of_records,  training_records;
        size_t  train_record_size, test_record_size, localstore_usage;
        size_t  train_block_size,  test_block_size;
        /* Data arrays */
        float  *training;
        float  *test;
        float  *local_distance_data;
        int    *classlabels;
        int    *assigned_classes;
        /* Loop indices */
        int     i, j;

        unsigned long elapsed_usec;

        if (argc < 6)
        {
                printf("usage: %s <labeled_pts> <points_to_label> <dim> <k_neighbors> <N_number_of_classes>\n",
                       argv[0]);
                printf("eg: %s 200000 5000 48 30 20\n", argv[0]);
                return 1;
        }

        TRAINING_RECORDS = atoi(argv[1]);
        TEST_RECORDS     = atoi(argv[2]);
        DIMENSIONS       = atoi(argv[3]);
        K                = atoi(argv[4]);
        CLASSES          = atoi(argv[5]);

        /* Compute training and test block sizes to friendly values  */
        test_record_size   = (DIMENSIONS + K) * sizeof(float);
        test_block_records = 1;
        test_block_size    = test_record_size;
        while (test_block_size < BLOCK_SIZE)
        {
                test_block_records *= 2;
                test_block_size    *= 2;
        }
        train_record_size   = (DIMENSIONS + 1) * sizeof(float);
        train_block_records = 1;
        train_block_size    = train_record_size;
        while (train_block_size < BLOCK_SIZE)
        {
                train_block_records *= 2;
                train_block_size    *= 2;
        }

        localstore_usage  = test_block_records * sizeof(int) +       /* class_labels */
		test_block_records * sizeof(int) +       /* assigned_classes */
		test_block_records * K * sizeof(float) + /* local_distance_data */
		train_block_size +                       /* train */
		test_block_size;                         /* data  */

        /* Check alignment of test slices */
        if (test_block_size % 16 != 0)
        {
                printf("That combination of <dim = %d> and <k_neighbors = %d> could lead to\n",
                       DIMENSIONS, K);
                printf("memory alignment problems.\n\nPoint size is %zd bytes, which would make a task block of %zd bytes.\n",
                       test_record_size, test_block_size);
                return 2;
        }

        printf("***********************************************\n");
        printf("***********************************************\n");
        printf("******                                   ******\n");
        printf("******         CellSs KNN DRIVER         ******\n");
        printf("******  BARCELONA SUPERCOMPUTING CENTER  ******\n");
        printf("******                                   ******\n");
        printf("******     ORIGINAL CBE VERSION FROM     ******\n");
        printf("******       OHIO STATE UNIVERSITY       ******\n");
        printf("******                                   ******\n");
        printf("***********************************************\n");
        printf("***********************************************\n\n");
        printf("BLOCK_SIZE is %d\n",                   BLOCK_SIZE);
        printf("Number of labeled samples is %d\n",    TRAINING_RECORDS);
        printf("Number of points to label is %d\n",    TEST_RECORDS);
        printf("Number of dimensions is %d\n",         DIMENSIONS);
        printf("Number of neighbors (k) is %d\n",      K);
        printf("Number of classes is %d\n\n",          CLASSES);
        printf("Training record size is %zd\n",        train_record_size);
        printf("Training block size is %zd\n",         train_block_size);
        printf("Training block contains %d records\n", train_block_records);
        printf("Test record size is %zd\n",            test_record_size);
        printf("Test block size is %zd\n",             test_block_size);
        printf("Test block contains %d records\n\n",   test_block_records);

        printf("Predicted localstore usage is %zd bytes\n\n", localstore_usage);

        /* Allocate space for data */
        training            = memalign(128, TRAINING_RECORDS * (DIMENSIONS + 1) * sizeof(float));
        test                = memalign(128, TEST_RECORDS * (DIMENSIONS + K) * sizeof(float));
        local_distance_data = memalign(128, TEST_RECORDS * K * sizeof(float));
        classlabels         = memalign(128, TEST_RECORDS * sizeof(int));
        assigned_classes    = memalign(128, TEST_RECORDS * sizeof(int));

        memset(local_distance_data, 0, TEST_RECORDS * K * sizeof(float));
        memset(assigned_classes,    0, TEST_RECORDS * sizeof(int));

        time_int(0);
        /* Make training data */
        for (i = 0;  i < TRAINING_RECORDS;  i++)
        {
                float x, y;

                for (j = 0;  j < DIMENSIONS;  j++)
                {
                        x = rand() / 10000000.0;
                        training[i * (DIMENSIONS + 1) + j] = x;
                }
                y = rand() % CLASSES;
                training[i * (DIMENSIONS + 1) + DIMENSIONS] = y; /* Set class */
        }
        /* Make unknown data */
        for (i = 0;  i < TEST_RECORDS;  i++)
        {
                for (j = 0;  j < DIMENSIONS;  j++)
                        test[i * (DIMENSIONS + K) + j]= rand() / 10000000.0;
                /* Set classes */
                for (j = 0;  j < K;  j++)
                        test[i * (DIMENSIONS + K) + DIMENSIONS + j] = -1;
        }
        time_int(1);

        time_int(0);

        /* Launch the distance calculation tasks */
        for (j = 0;  j < TRAINING_RECORDS;  j += train_block_records)
        {
                training_records = 
			((j + train_block_records > TRAINING_RECORDS)
			 ? (TRAINING_RECORDS - j)
			 : train_block_records);

                for (i = 0;  i < TEST_RECORDS;  i += test_block_records)
                {
                        number_of_records = 
				(i + test_block_records > TEST_RECORDS 
				 ? TEST_RECORDS - i 
				 : test_block_records);

                        compare(K, DIMENSIONS, CLASSES, training_records, 
				number_of_records, 
				&training[j * (DIMENSIONS + 1)], 
				&test[i * (DIMENSIONS + K)],
                                &local_distance_data[i * K],
				&assigned_classes[i]);
                }
        }
        /* And now the class selection tasks */
        for (i = 0;  i < TEST_RECORDS;  i += test_block_records)
        {
                number_of_records = 
			((i + test_block_records > TEST_RECORDS) 
			 ? (TEST_RECORDS - i) 
			 : test_block_records);

                get_best_classlabels(K, DIMENSIONS, CLASSES, number_of_records,
				     &test[i * (DIMENSIONS + K)], 
				     &classlabels[i]);
        }
#pragma omp taskwait
        elapsed_usec = time_int(1);

        printf("sample of the computed class labels ... \n");
        for (i = 1;  i < TEST_RECORDS;  i++)
        {
                if (i % 50 == 0)
                        printf("rec %6d = cls %3d    ", i, classlabels[i]);
                if (i % 200 == 0)
                        printf("\n");
        }
        printf("\n\n");

        free(training);
        free(test);
        free(classlabels);

        printf("par_sec_time_us:%lu\n",elapsed_usec);

        return 0;
}


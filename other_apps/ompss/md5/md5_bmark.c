/*
*   MD5 Benchmark
*   -------------
*   File: md5_bmark.c
*
*   This is the main file for the md5 benchmark kernel. This benchmark was
*   written as part of the StarBENCH benchmark suite at TU Berlin. It performs
*   MD5 computation on a number of self-generated input buffers in parallel,
*   automatically measuring execution time.
*
*   Copyright (C) 2011 Michael Andersch
*
*   This program is free software: you can redistribute it and/or modify
*   it under the terms of the GNU General Public License as published by
*   the Free Software Foundation, either version 3 of the License, or
*   (at your option) any later version.
*
*   This program is distributed in the hope that it will be useful,
*   but WITHOUT ANY WARRANTY; without even the implied warranty of
*   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*   GNU General Public License for more details.
*
*   You should have received a copy of the GNU General Public License
*   along with this program.  If not, see <http://www.gnu.org/licenses/>.
*
*/

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <unistd.h>
#include <string.h>
#include <sys/time.h>

#include "md5.h"
#include "md5_bmark.h"

typedef struct timeval timer;

#define TIME(x) gettimeofday(&x, NULL)

/* Function declarations */
int initialize(md5bench_t* args);
int finalize(md5bench_t* args);
void run(md5bench_t* args);
void process(uint8_t* in, uint8_t* out, int bufsize);
void listInputs();
long timediff(timer* starttime, timer* finishtime);

static char* usage = "Usage: %s <options>\n"
                     "-i inputnum          choose input data set\n"
                     "-c iterations        specify benchmark iterations\n"
                     "-s                   show available input sets\n"
                     "-h                   this help text\n";

// Input configurations
static data_t datasets[] = {
    {64, 512, 0},
    {64, 1024, 0},
    {64, 2048, 0},
    {64, 4096, 0},
    {128, 1024*512, 1},
    {128, 1024*1024, 1},
    {128, 1024*2048, 1},
    {128, 1024*4096, 1},
	// Lluc: Ad-hoc input set with 1000 tasks of 32KB each
    {1000, (32*1024) / (sizeof(uint8_t) * 2), 1},
    {2000, (32*1024) / (sizeof(uint8_t) * 2), 1},
    {6400, (32*1024) / (sizeof(uint8_t) * 2), 1},
    {1000, (64*1024) / (sizeof(uint8_t) * 2), 1},
    {1000, (128*1024) / (sizeof(uint8_t) * 2), 1},
    {1000, (256*1024) / (sizeof(uint8_t) * 2), 1},
    {1000, (512*1024) / (sizeof(uint8_t) * 2), 1},
    {1000, (1024*1024) / (sizeof(uint8_t) * 2), 1},
    {1000, (2048*1024) / (sizeof(uint8_t) * 2), 1},
};

/*
*   Function: initialize
*   --------------------
*   To initialize the benchmark parameters. Generates the input buffers from random data.
*/
int initialize(md5bench_t* args) {
    int index = args->input_set;
    if(index < 0 || index >= sizeof(datasets)/sizeof(datasets[0])) {
        fprintf(stderr, "Invalid input set specified! Clamping to set 0\n");
        index = 0;
    }

    args->numinputs = datasets[index].numbufs;
    args->size = datasets[index].bufsize;
    args->inputs = (uint8_t**)calloc(args->numinputs, sizeof(uint8_t*));
    args->out = (uint8_t*)calloc(args->numinputs, DIGEST_SIZE);
    if(args->inputs == NULL || args->out == NULL) {
        fprintf(stderr, "Memory Allocation Error\n");
        return -1;
    }

    fprintf(stderr, "Reading input set: %d buffers, %d bytes per buffer\n", datasets[index].numbufs, datasets[index].bufsize);

    // Now the input buffers need to be generated, for replicability, use same seed
    srand(datasets[index].rseed);

    for(int i = 0; i < args->numinputs; i++) {
        args->inputs[i] = (uint8_t*)malloc(sizeof(uint8_t)*datasets[index].bufsize);
        uint8_t* p = args->inputs[i];
        if(p == NULL) {
            fprintf(stderr, "Memory Allocation Error\n");
            return -1;
        }
        for(int j = 0; j < datasets[index].bufsize; j++)
            *p++ = rand() % 255;
    }

    return 0;
}

/*
*   Function: process
*   -----------------
*   Processes one input buffer, delivering the digest into out.
*/
#pragma omp task in([bufsize] in) out([bufsize] out)
void process(uint8_t* in, uint8_t* out, int bufsize) {
    MD5_CTX context;
    uint8_t digest[16];
    
    MD5_Init(&context);
    MD5_Update(&context, in, bufsize);
    MD5_Final(digest, &context);

    memcpy(out, digest, DIGEST_SIZE);
}

/*
*   Function: run
*   -------------
*   Main benchmarking function. If called, processes buffers with MD5
*   until no more buffers available. The resulting message digests
*   are written into consecutive locations in the preallocated output
*   buffer.
*/
void run(md5bench_t* args) {
	for(int i = 0; i < args->iterations; i++) {

        int buffers_to_process = args->numinputs;
        int next = 0;
        uint8_t** in = args->inputs;
        uint8_t* out = args->out;

        while(buffers_to_process > 0) {
            process(in[next], out+next*DIGEST_SIZE, args->size);
            next++;
            buffers_to_process--;
        }
        #pragma omp taskwait
	}
}

/*
*   Function: finalize
*   ------------------
*   Cleans up memory used by the benchmark for input and output buffers.
*/
int finalize(md5bench_t* args) {

    char outname[] = "output.txt";
    char buffer[64];
    FILE* fp;

    if(args->outflag) {
        fp = fopen(outname, "w");

        for(int i = 0; i < args->numinputs; i++) {
            sprintf(buffer, "Buffer %d has checksum ", i);
            fwrite(buffer, sizeof(char), strlen(buffer)+1, fp);
            for(int j = 0; j < DIGEST_SIZE*2; j+=2) {
                sprintf(buffer+j, "%x", args->out[DIGEST_SIZE*i+j/2] & 0xf);
                sprintf(buffer+j+1, "%x", args->out[DIGEST_SIZE*i+j/2] & 0xf0);
            }
            buffer[32] = '\0';
            fwrite(buffer, sizeof(char), 32, fp);
            fputc('\n', fp);
        }
        
        fclose(fp);
    }

    if(args->inputs) {
        for(int i = 0; i < args->numinputs; i++) {
            if(args->inputs[i])
                free(args->inputs[i]);
        }

        free(args->inputs);
    }

    if(args->out)
        free(args->out);

    return 0;
}

/*
*   Function: listInputs
*   --------------------
*   Lists all available input configurations of the benchmark on stderr.
*/
void listInputs() {
    double datasize_mb;
    fprintf(stderr, "Available input configurations:\n");

    for(int i = 0; i < sizeof(datasets)/sizeof(datasets[0]); i++) {
        datasize_mb = (double)(datasets[i].numbufs * datasets[i].bufsize)/1000000.0;
        fprintf(stderr, "Index %d: %d buffers of %d bytes, totalling %.1fMB\n", i, datasets[i].numbufs, datasets[i].bufsize, datasize_mb);
    }
}

/*
*   Function: timediff
*   ------------------
*   Compute the difference between timers starttime and finishtime in msecs.
*/
long timediff(timer* starttime, timer* finishtime)
{
    long msec;
    msec=(finishtime->tv_sec-starttime->tv_sec)*1000;
    msec+=(finishtime->tv_usec-starttime->tv_usec)/1000;
    return msec;
}

/** MAIN **/
int main(int argc, char** argv) {

    int opt;
    extern char* optarg;
    extern int optind;

    timer io_start, b_start, b_end;

    md5bench_t args;
    args.input_set = 0;
    args.iterations = 1;
    args.outflag = 0;

    // Who we are
    fprintf(stderr, "StarBench - MD5 Kernel\n");

    // Parse command line options
    while ( (opt=getopt(argc,argv,"i:c:ohs")) != EOF) {
        switch (opt) {
            case 'i':
                args.input_set = atoi(optarg);
                break;
            case 'c':
                args.iterations = atoi(optarg);
                break;
            case 'h': 
                fprintf(stderr, usage, argv[0]);
                return -1;
                break;
            case 's':
                listInputs();
                return -1;
                break;
            case 'o':
                args.outflag = 1;
                break;
            default:
				fprintf(stderr, usage, argv[0]);
                return -1;
                break;
        }
    }

	if(args.iterations < 1) {
		fprintf(stderr, "Bad argument given\n");
		exit(EXIT_FAILURE);
	}

    TIME(io_start);

    // Parameter initialization
    if(initialize(&args)) {
        fprintf(stderr, "Initialization Error\n");
        exit(EXIT_FAILURE);
    }

    TIME(b_start);
    
    run(&args);

    TIME(b_end);

    // Free memory
    if(finalize(&args)) {
		fprintf(stderr, "Finalization Error\n");
		exit(EXIT_FAILURE);
	}

    double io_time = (double)timediff(&io_start, &b_start)/1000;
    double b_time = (double)timediff(&b_start, &b_end)/1000;

    printf("\nI/O time: %.3f\nTime: %.3f\n", io_time, b_time);

    return 0;
}

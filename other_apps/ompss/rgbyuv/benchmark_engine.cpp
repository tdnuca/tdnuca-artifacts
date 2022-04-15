/*
*   File: benchmark_engine.cpp
*   --------------------------
*   Implementations of the benchmarkengine class functions.
*/

#include "benchmark_engine.h"

/*
 *  Intitialization of the benchmarking object. Creates the control objects for the two kernels.
 */
bool BenchmarkEngine::init(string srcname, string destname, unsigned int angle) {
    re = new RotateEngine();
    ce = new ConvertEngine();
    if(!re->init(srcname, destname, angle))
        return false;
    if(!ce->init(re->getOutput()))
        return false;
    return true;
};

/* Wrappers for OmpSs tasking */

#pragma omp task out(*dep)
void computeRowTask(int row, int height, int width, int xot, int yot, int xos, int yos, int ra, char* dep, RotateEngine* re) {
	re->computeRow(row, height, width, xot, yot, xos, yos, ra, dep);
}

#pragma omp task in(*dep)
void convertLineTask(int line, char* dep, ConvertEngine *ce) {
	ce->convertLine(line, dep);
}

/*
 *  Manages the adding of tasks to the runtime taskgraph.
 */
void BenchmarkEngine::run() {
    Image* input = re->getInput();

    unsigned int height = input->getHeight();
    unsigned int width = input->getWidth();
    //unsigned int depth = input->getDepth();
    float x_offset_source = (float)width / 2.0;
    float y_offset_source = (float)height / 2.0;
    unsigned int rev_angle = 360 - re->getangle();
    float x_offset_target = (float)re->gettargetw()/2.0;
    float y_offset_target = (float)re->gettargeth()/2.0;

    char* dependencies = new char[re->gettargeth()];
    for(int i = 0; i < re->gettargeth(); i++) {
        computeRowTask(i, re->gettargeth(), re->gettargetw(), x_offset_target, y_offset_target, x_offset_source, y_offset_source, rev_angle, &dependencies[i], re);
        convertLineTask(i, &dependencies[i], ce);
    }
    #pragma omp barrier
    delete[] dependencies;
}

/*
 *  Cleans up all memory used during the benchmark.
 */
void BenchmarkEngine::finish() {
    re->finish();
    ce->finish();
    delete re;
    delete ce;
}  

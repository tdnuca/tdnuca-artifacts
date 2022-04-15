/*
*   File: benchmark_engine.h
*   ------------------------
*   Header file containing the definition of the main benchmarking class.
*/

/**********************************************************************************
                INCLUDES & DEFINES
***********************************************************************************/
#ifndef B_ENGINE_H
#define B_ENGINE_H

#include "rotation_engine.h"
#include "convert_engine.h"
#include <string>

using namespace std;

class BenchmarkEngine {
public:
    bool init(string srcname, string destname, unsigned int angle);
    void run();
    void finish();
private:
    RotateEngine* re;
    ConvertEngine* ce;
};

#endif // B_ENGINE_H
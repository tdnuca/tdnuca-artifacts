/*
*   File: convert_engine.h
*   ----------------------
*   Contains the definition of the class encapsulating the color conversion kernel.
*/

/**********************************************************************************
                INCLUDES & DEFINES
***********************************************************************************/
#ifndef C_ENGINE_H
#define C_ENGINE_H

#include <stdint.h>
#include <stdio.h>
#include "image.h"
#include <sys/time.h>
#include <cmath>

class ConvertEngine {
public:
    ConvertEngine();
    bool init(Image* in);
    void run();
    void finish();
    void convertLine(int line, char* dep);
private:
    Image* input, *output;
    int width, height;
    bool initialized, done;
};

#endif // C_ENGINE_H
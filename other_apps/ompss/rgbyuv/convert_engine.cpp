/*
*   File: convert_engine.cpp
*   ------------------------
*   Contains implementations of the convertengine class functions.
*/

#include "convert_engine.h"

typedef struct timeval timer;

ConvertEngine::ConvertEngine() {
    done = false;
    initialized = false;
};

bool ConvertEngine::init(Image* in) {
    input = in;
    width = input->getWidth();
    height = input->getHeight();
    output = new Image();
    output->createImageFromTemplate(width, height, 3);
    initialized = true;
    return true;
}

//#pragma omp task input(*dep)
void ConvertEngine::convertLine(int line, char* dep) {
    //char d = *dep;

    uint8_t R, G, B;
    Pixel p; // Output, components used for YUV

    for(int j = 0; j < width; j++) {
        R=input->getPixelAt(j,line).r;
        G=input->getPixelAt(j,line).g;
        B=input->getPixelAt(j,line).b;

        p.r = round(0.256788*R+0.504129*G+0.097906*B) + 16;
        p.g = round(-0.148223*R-0.290993*G+0.439216*B) + 128;
        p.b = round(0.439216*R-0.367788*G-0.071427*B) + 128;
            
        output->setPixelAt(j,line,&p);
    }

    if(line == height-1)
        done = true;
}

void ConvertEngine::run() {
    for(int i = 0; i < height; i++) {
        convertLine(i, NULL);
    }
    done = true;
    return;
}

void ConvertEngine::finish() {
    char uvheader[256];
    snprintf(uvheader, (size_t)255, "P6\n%d %d\n%d\n", width, height, RGB_MAX_COLOR);
    
    // Write the full YUV image
    FILE* fp = fopen("yuvout.ppm", "w");
    if (fp == NULL) {
        fprintf(stderr, "Cannot open output file\n");
        output->clean();
        delete output;
    }
        
    fprintf(fp, "%s", uvheader);
    for(int i = 0; i < height; i++) {
        for(int j = 0; j < width; j++) {
            fputc(output->getPixelAt(j,i).r, fp);
            fputc(output->getPixelAt(j,i).g, fp);
            fputc(output->getPixelAt(j,i).b, fp);
        }
    }
    fclose(fp);

    // Input image is cleaned up by the rotation engine, so no need to delete it
    output->clean();
    delete output;
}

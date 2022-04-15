Application: TinyJPEG

This application benchmark performs decoding of JPEG images with fixed encoding of 2x2 MCU size and YUV color, producing RGB .tga output files. Images can contain up to one RST marker per MCU line.

Installation:

To install the benchmark, navigate to the directory this file is located in, open up a terminal and simply type 'make'. For certain architectures 
or special compilation options, you might need to change compilation parameters in the makefile.

Usage:

You may execute the benchmark by navigating to this directory after compilation and typing

./tinyjpeg 

The benchmark will then list a number of parameters you have to enter in order to execute. Specifying --benchmark on the commandline will automatically produce a measurement of the time it took to decode the picture.

Benchmark versions:

Serial
POSIX Threads
OpenMP SuperScalar

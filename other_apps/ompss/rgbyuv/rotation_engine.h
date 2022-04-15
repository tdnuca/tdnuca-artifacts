/*
*	"Image Rotate", a program to rotate images by a user-specifyable angle.
*	
*	Copyright (C) 2010 Michael Andersch
*	
*	This program is free software; you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation; either version 3 of the License, or (at your option) any later version.
*	
*	This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
*	
*	You should have received a copy of the GNU General Public License along with this program; if not, see <http://www.gnu.org/licenses/>.
*/

/*
*	File: rotation_engine.h
*	-----------------------
*	Header file for the rotation core. Contains the rotateengine object
*	definition.
*/

/**********************************************************************************
				INCLUDES & DEFINES
***********************************************************************************/
#ifndef R_ENGINE_H
#define R_ENGINE_H

#include <iostream>
#include <fstream>
#include <stdio.h>
#include <cmath>
#include <float.h>
#include "image.h"

using namespace std;

#define BLOCK 16

/*
*	Class: RotateEngine
*	-------------------
*   Container object for the benchmark. Contains benchmark state, i.e.
*   input and output images as well as additional kernel data.
*/
class RotateEngine {
    public:
		RotateEngine();
		void run();
		void finish();
		bool init(string srcname, string destname, unsigned int angle);
        void printRotationState();
        Image* getOutput();
        Image* getInput();
        int gettargeth();
        int gettargetw();
        unsigned int getangle();
        void computeRow(int row, int maxrow, int rowwidth, float xot, float yot, float xos, float yos, unsigned rev_angle, char* dep);
    private:
        /* Variables */
        string srcname, destname;
		Image input, output;
		unsigned int angle;
        int target_h, target_w;
        bool initialized, done;
		Coord ul, ur, ll, lr, c1, c2, c3, c4;
        /* Functions */
        bool writeOutImage();
		void rotatePoint(float *pt, float *target, unsigned int angle);
		double round(double num, int digits);
		int computeTargetHeight();
		int computeTargetWidth();
		float findMax(float* seq);
		float findMin(float* seq);
		void filter(Pixel* colors, Pixel* dest, Coord* sample_pos);
		void interpolateLinear(Pixel* a, Pixel* b, Pixel* dest, float weight);
};

#endif

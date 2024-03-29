#include "misc.hh"
#include <cstdlib>
#include <limits>
#include <sys/time.h>
#include <string>
#include <sstream>

double wtime()
{
  struct timeval tv;
  gettimeofday(&tv, 0);
  
  return tv.tv_sec+1e-6*tv.tv_usec;
}

algoparam_t::algoparam_t()
  : m(NULL),
    uvis(NULL),
    heatsrcs(NULL)
{}

algoparam_t::~algoparam_t()
{
  if (m != NULL)
    delete m;
  
  if (heatsrcs != NULL)
    delete[] heatsrcs;
  
  if (uvis != NULL)
    delete[] uvis;
}

int algoparam_t::read_input( std::istream& infile )
{
  algoparam_t * param = this;

  unsigned i;
  std::string buf;

  std::getline(infile, buf);
  std::istringstream ss(buf);
  ss >> maxiter;
    
  buf.clear();
  std::getline(infile, buf);
  ss.str(buf);
  ss >> resolution;

  param->visres = param->resolution;
  //param->visres = 256;

  buf.clear();
  std::getline(infile, buf);
  ss.str(buf);
  ss >> algorithm;

  buf.clear();
  std::getline(infile, buf);
  ss.str(buf);
  ss >> numsrcs;

  //(param->heatsrcs) = 
  //  (heatsrc_t*) malloc( sizeof(heatsrc_t) * (param->numsrcs) );
  heatsrcs = new heatsrc_t[numsrcs];
  
  for( i=0; i<param->numsrcs; i++ ) {
    buf.clear();
    std::getline(infile, buf);
    ss.str(buf);
    ss >> heatsrcs[i].posx >> heatsrcs[i].posy >> heatsrcs[i].range >> heatsrcs[i].temp;
  }
  return 1;
}

int algoparam_t::init(unsigned long NB)
{
  if (resolution%NB != 0) {
    std::cerr << "The number of blocks must divide the resolution!" << std::endl;
    exit(1);
  }
  m = new Matrix(NB, resolution/NB, algorithm == 0);
  uvis = new double[visres*visres];
  std::fill(&uvis[0], &uvis[visres*visres - 1], 0.0);

  for (unsigned long i = 0; i < numsrcs; ++i) {
    m->setBorder(heatsrcs[i]);
  }

  return 1;
}

void algoparam_t::print_params() const
{
  fprintf(stdout, "Iterations        : %u\n", maxiter);
  fprintf(stdout, "Resolution        : %u\n", resolution);
  fprintf(stdout, "Algorithm         : %d (%s)\n",
	  algorithm,
	  (algorithm == 0) ? "Jacobi" : (algorithm ==1) ? "Gauss-Seidel" : "Red-Black" );
  fprintf(stdout, "Num. Heat sources : %u\n", numsrcs);

  for(unsigned i = 0; i < numsrcs; i++){
    fprintf(stdout, "  %2d: (%2.2f, %2.2f) %2.2f %2.2f \n",
	    i+1,
	    heatsrcs[i].posx,
	    heatsrcs[i].posy,
	    heatsrcs[i].range,
	    heatsrcs[i].temp );
  }

  std::fflush(stdout);
}

void algoparam_t::write_image( std::FILE * f )
{
  double *u = uvis;
  int sizex = visres;
  int sizey = visres;

  double min = std::numeric_limits<double>::max();
  double max = -std::numeric_limits<double>::max();

  // find minimum and maximum 
  for (int i=0; i < sizey; i++) {
    for (int j=0; j < sizex; j++) {
      if (u[i*sizex+j] > max)
        max = u[i*sizex+j];
      if (u[i*sizex+j] < min)
        min = u[i*sizex+j];
    }
  }
  

  fprintf(f, "P3\n");
  fprintf(f, "%u %u\n", sizex, sizey);
  fprintf(f, "%u\n", 255);

  for (int i = 0; i < sizey; i++) {
    for (int j = 0; j < sizex; j++) {
      int k = (int)(255.0*(u[i*sizex + j]-min)/(max-min));
      fprintf(f, "%d %d %d  ", magma[k].r, magma[k].g, magma[k].b);
    }
    fprintf(f, "\n");
  }
}

int algoparam_t::coarsen()
{
  unsigned long NB = m->NB;
  unsigned long BS = m->BS;
  for (unsigned long i = 0; i < m->NB; ++i) {
    for (unsigned long j = 0; j < m->NB; ++j) {
      unsigned long extra = 0;
      for (unsigned ii = 0; ii < m->BS; ++ii) {
	for (unsigned jj = 0; jj < m->BS; ++jj) {
	  uvis[i*NB*BS*BS + j*BS + extra + jj] = m->block(i, j)[ii*m->BS+jj];
	}
	extra += resolution;
      }
    }
  }

  return 1;
}

colormap_t magma = {{0, 0, 4},
                    {1, 0, 5},
                    {1, 1, 6},
                    {1, 1, 8},
                    {2, 1, 10},
                    {2, 2, 12},
                    {2, 2, 14},
                    {3, 3, 16},
                    {4, 3, 18},
                    {4, 4, 20},
                    {5, 4, 22},
                    {6, 5, 24},
                    {6, 5, 26},
                    {7, 6, 28},
                    {8, 7, 30},
                    {9, 7, 32},
                    {10, 8, 34},
                    {11, 9, 36},
                    {12, 9, 38},
                    {13, 10, 41},
                    {14, 11, 43},
                    {16, 11, 45},
                    {17, 12, 47},
                    {18, 13, 50},
                    {19, 13, 52},
                    {20, 14, 54},
                    {21, 14, 57},
                    {23, 15, 59},
                    {24, 15, 61},
                    {25, 16, 64},
                    {26, 16, 66},
                    {28, 16, 68},
                    {29, 17, 71},
                    {30, 17, 73},
                    {32, 17, 76},
                    {33, 17, 78},
                    {35, 18, 81},
                    {36, 18, 83},
                    {38, 18, 86},
                    {39, 18, 88},
                    {41, 17, 90},
                    {42, 17, 93},
                    {44, 17, 95},
                    {46, 17, 97},
                    {47, 17, 99},
                    {49, 17, 101},
                    {51, 16, 103},
                    {52, 16, 105},
                    {54, 16, 107},
                    {56, 16, 109},
                    {58, 15, 111},
                    {59, 15, 112},
                    {61, 15, 113},
                    {63, 15, 115},
                    {65, 15, 116},
                    {66, 15, 117},
                    {68, 15, 118},
                    {70, 16, 119},
                    {71, 16, 120},
                    {73, 16, 121},
                    {75, 17, 122},
                    {76, 17, 122},
                    {78, 17, 123},
                    {79, 18, 124},
                    {81, 18, 124},
                    {83, 19, 125},
                    {84, 19, 125},
                    {86, 20, 126},
                    {87, 21, 126},
                    {89, 21, 126},
                    {91, 22, 127},
                    {92, 22, 127},
                    {94, 23, 127},
                    {95, 24, 128},
                    {97, 24, 128},
                    {98, 25, 128},
                    {100, 26, 128},
                    {101, 26, 129},
                    {103, 27, 129},
                    {105, 28, 129},
                    {106, 28, 129},
                    {108, 29, 129},
                    {109, 30, 129},
                    {111, 30, 130},
                    {112, 31, 130},
                    {114, 31, 130},
                    {116, 32, 130},
                    {117, 33, 130},
                    {119, 33, 130},
                    {120, 34, 130},
                    {122, 34, 130},
                    {123, 35, 130},
                    {125, 36, 130},
                    {127, 36, 130},
                    {128, 37, 130},
                    {130, 37, 130},
                    {131, 38, 130},
                    {133, 38, 130},
                    {134, 39, 130},
                    {136, 40, 130},
                    {138, 40, 130},
                    {139, 41, 130},
                    {141, 41, 129},
                    {142, 42, 129},
                    {144, 42, 129},
                    {146, 43, 129},
                    {147, 43, 129},
                    {149, 44, 129},
                    {151, 44, 129},
                    {152, 45, 128},
                    {154, 46, 128},
                    {155, 46, 128},
                    {157, 47, 128},
                    {159, 47, 127},
                    {160, 48, 127},
                    {162, 48, 127},
                    {164, 49, 127},
                    {165, 49, 126},
                    {167, 50, 126},
                    {169, 50, 126},
                    {170, 51, 125},
                    {172, 51, 125},
                    {174, 52, 124},
                    {175, 52, 124},
                    {177, 53, 124},
                    {178, 53, 123},
                    {180, 54, 123},
                    {182, 54, 122},
                    {183, 55, 122},
                    {185, 56, 121},
                    {187, 56, 121},
                    {188, 57, 120},
                    {190, 57, 120},
                    {192, 58, 119},
                    {193, 59, 118},
                    {195, 59, 118},
                    {196, 60, 117},
                    {198, 60, 117},
                    {200, 61, 116},
                    {201, 62, 115},
                    {203, 63, 115},
                    {204, 63, 114},
                    {206, 64, 113},
                    {208, 65, 112},
                    {209, 66, 112},
                    {211, 66, 111},
                    {212, 67, 110},
                    {214, 68, 109},
                    {215, 69, 109},
                    {217, 70, 108},
                    {218, 71, 107},
                    {220, 72, 106},
                    {221, 73, 106},
                    {222, 74, 105},
                    {224, 75, 104},
                    {225, 76, 103},
                    {226, 77, 102},
                    {228, 78, 102},
                    {229, 79, 101},
                    {230, 81, 100},
                    {231, 82, 99},
                    {233, 83, 99},
                    {234, 84, 98},
                    {235, 86, 97},
                    {236, 87, 97},
                    {237, 89, 96},
                    {238, 90, 95},
                    {239, 92, 95},
                    {240, 93, 94},
                    {241, 95, 94},
                    {242, 97, 93},
                    {242, 98, 93},
                    {243, 100, 93},
                    {244, 102, 93},
                    {245, 104, 92},
                    {245, 105, 92},
                    {246, 107, 92},
                    {247, 109, 92},
                    {247, 111, 92},
                    {248, 113, 92},
                    {248, 114, 92},
                    {249, 116, 92},
                    {249, 118, 93},
                    {250, 120, 93},
                    {250, 122, 93},
                    {250, 124, 94},
                    {251, 126, 94},
                    {251, 128, 95},
                    {251, 129, 95},
                    {252, 131, 96},
                    {252, 133, 96},
                    {252, 135, 97},
                    {253, 137, 98},
                    {253, 139, 99},
                    {253, 141, 99},
                    {253, 143, 100},
                    {253, 145, 101},
                    {254, 147, 102},
                    {254, 149, 103},
                    {254, 150, 104},
                    {254, 152, 105},
                    {254, 154, 106},
                    {254, 156, 107},
                    {255, 158, 108},
                    {255, 160, 109},
                    {255, 162, 111},
                    {255, 164, 112},
                    {255, 165, 113},
                    {255, 167, 114},
                    {255, 169, 115},
                    {255, 171, 117},
                    {255, 173, 118},
                    {255, 175, 119},
                    {255, 177, 121},
                    {255, 179, 122},
                    {255, 180, 124},
                    {255, 182, 125},
                    {255, 184, 126},
                    {255, 186, 128},
                    {255, 188, 129},
                    {255, 190, 131},
                    {255, 191, 132},
                    {255, 193, 134},
                    {255, 195, 135},
                    {255, 197, 137},
                    {255, 199, 139},
                    {255, 201, 140},
                    {255, 203, 142},
                    {255, 204, 143},
                    {255, 206, 145},
                    {255, 208, 147},
                    {255, 210, 148},
                    {255, 212, 150},
                    {255, 214, 152},
                    {255, 215, 153},
                    {255, 217, 155},
                    {254, 219, 157},
                    {254, 221, 159},
                    {254, 223, 160},
                    {254, 225, 162},
                    {254, 226, 164},
                    {254, 228, 166},
                    {254, 230, 167},
                    {254, 232, 169},
                    {254, 234, 171},
                    {254, 236, 173},
                    {253, 237, 175},
                    {253, 239, 177},
                    {253, 241, 179},
                    {253, 243, 180},
                    {253, 245, 182},
                    {253, 246, 184},
                    {253, 248, 186},
                    {253, 250, 188},
                    {253, 252, 190},
                    {253, 254, 192}};

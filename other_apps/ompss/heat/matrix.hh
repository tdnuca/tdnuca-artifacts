#ifndef _MATRIX_HH
#define _MATRIX_HH

struct heatsrc_t
{
  double posx;
  double posy;
  double range;
  double temp;
};

struct Matrix
{
  unsigned long NB;
  unsigned long BS;
  int num_sockets;
  int socket_change;

private:

  double **u[2]; // u matrix
  double **vertical[2];   // left and right halos (containers)
  double **horizontal[2]; // top and bottom halos (containers)
  const bool dbuffer;
  int current;

public:
  
  Matrix(unsigned long NB, unsigned long BS, bool dbuffer=false);

  ~Matrix();
  
  void setBorder(const heatsrc_t &src);

  inline double * block(unsigned long i, unsigned long j)
  {
    return u[dbuffer ? current : 0][i*NB + j];
  }

  inline double * oldblock(unsigned long i, unsigned long j)
  {
    return u[dbuffer ? 1-current : 0][i*NB + j];
  }

  inline double * top(unsigned long i, unsigned long j)
  {
    return horizontal[dbuffer ? 1-current : 0][i*NB + j];
  }

  inline double * bottom(unsigned long i, unsigned long j)
  {
    return horizontal[dbuffer ? current : 0][(i+1)*NB + j];
  }

  inline double * left(unsigned long i, unsigned long j)
  {
    return vertical[dbuffer ? 1-current : 0][i*(NB+1) + j];
  }

  inline double * right(unsigned long i, unsigned long j)
  {
    return vertical[dbuffer ? current : 0][i*(NB+1) + j+1];
  }

  inline Matrix & operator++()
  {
    current = (current + 1)%2;

    return *this;
  }
};
#endif

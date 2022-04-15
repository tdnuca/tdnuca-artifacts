#include "misc.hh"

unsigned long wtime()
{
  struct timeval tv;
  gettimeofday(&tv, 0);

  return tv.tv_sec*1000000+tv.tv_usec;
}

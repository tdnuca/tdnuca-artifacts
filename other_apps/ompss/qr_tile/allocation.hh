#ifndef __ALLOCATION_H
#define __ALLOCATION_H

#include <unistd.h>

#ifndef _PAGE_SIZE
#define _PAGE_SIZE sysconf(_SC_PAGESIZE)
#endif

#include <cstdlib>

template<typename T>
void allocate(T** a, size_t size)
{
  *a = static_cast<T*>(malloc(size*sizeof(T)));
}

template<typename T>
void allocate_aligned(T** a, size_t size)
{
  if (posix_memalign(reinterpret_cast<void **>(a), _PAGE_SIZE, size*sizeof(T)) != 0) {
    *a = 0;
  }
}

#endif

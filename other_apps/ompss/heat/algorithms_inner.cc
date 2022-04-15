#include "algorithms_inner.hh"
#include <algorithm>
double inner_heat_step(const unsigned long BS,
		       double * mu,
		       double * top, const bool wr_top,
		       double * bottom, const bool wr_bottom,
		       double * left, const bool wr_left,
		       double * right, const bool wr_right,
		       const bool check)
{
  typedef double rm[BS][BS];
  rm & u = *reinterpret_cast<rm *>(mu);
  double my_sum = 0.0;
  double diff;
	
  // TOP LEFT
  double un = 0.25*( ( left[0] +
		       u[0][1] ) +
		     ( top[0] +
		       u[1][0] ) );
  if (check) {
    diff = un - u[0][0];
    my_sum += diff*diff;
  }
  u[0][0] = un;
	
  // TOP ROW, COL jj=1..BS-2
  for (unsigned long jj = 1; jj < BS-1; ++jj) {
    un = 0.25*( ( u[0][jj-1] +
		  u[0][jj+1] ) +
		( top[jj] +
		  u[1][jj] ) );
	  
    if (check) {
      diff = un - u[0][jj];
      my_sum += diff*diff;
    }
    u[0][jj] = un;
  }
	
  // TOP RIGHT
  un = 0.25*( ( u[0][BS-2] +
		right[0] ) +
	      ( top[BS-1] +
		u[1][BS-1] ) );
  if (check) {
    diff = un - u[0][BS-1];
    my_sum += diff*diff;
  }
  u[0][BS-1] = un;
	
  // ROW ii=1..BS-2
  for (unsigned long ii = 1; ii < BS-1; ++ii) {
	  
    // ROW ii, COL 0
    un = 0.25*( ( left[ii] +
		  u[ii][1] ) +
		( u[ii-1][0] +
		  u[ii+1][0] ) );
    if (check) {
      diff = un - u[ii][0];
      my_sum += diff*diff;
    }
    u[ii][0] = un;

    // ROW ii, COL jj=1..BS-2
    for (unsigned long jj = 1; jj < BS-1; ++jj) {
      // ROW ii, COL jj
      un = 0.25*( ( u[ii][jj-1] +
		    u[ii][jj+1] ) +
		  ( u[ii-1][jj] +
		    u[ii+1][jj] ) );
      if (check) {
	diff = un - u[ii][jj];
	my_sum += diff*diff;
      }
	    
      u[ii][jj] = un;
    }

    // ROW ii, COL BS-1
    un = 0.25*( ( u[ii][BS-2] +
		  right[ii] ) +
		( u[ii-1][BS-1] +
		  u[ii+1][BS-1] ) );
    if (check) {
      diff = un - u[ii][BS-1];
      my_sum += diff*diff;
    }
    u[ii][BS-1] = un;
  }
	
  // BOTTOM LEFT
  un = 0.25*( ( left[BS-1] +
		u[BS-1][1] ) +
	      ( u[BS-2][0] +
		bottom[0] ) );
  if (check) {
    diff = un - u[BS-1][0];
    my_sum += diff*diff;
  }
  u[BS-1][0] = un;

  // BOTTOM ROW, COL jj=1..BS-2
  for (unsigned long jj = 1; jj < BS-1; ++jj) {
    un = 0.25*( ( u[BS-1][jj-1] +
		  u[BS-1][jj+1] ) +
		( u[BS-2][jj] +
		  bottom[jj] ) );
	  
    if (check) {
      diff = un - u[BS-1][jj];
      my_sum += diff*diff;
    }
    u[BS-1][jj] = un;
  }
	
  // BOTTOM RIGHT
  un = 0.25*( ( u[BS-1][BS-2] +
		right[BS-1] ) +
	      ( u[BS-2][BS-1] +
		bottom[BS-1] ) );
  if (check) {
    diff = un - u[BS-1][BS-1];
    my_sum += diff*diff;
  }
  u[BS-1][BS-1] = un;

  // WRITE TO THE HALOS IF NEEDED
  if (wr_top)
    std::copy(&u[0][0], &u[0][BS], top);

  if (wr_bottom)
    std::copy(&u[BS-1][0], &u[BS-1][BS], bottom);

  if (wr_left)
    for (unsigned long i = 0; i < BS; ++i)
      left[i] = u[i][0];

  if (wr_right)
    for (unsigned long i = 0; i < BS; ++i)
      right[i] = u[i][BS-1];
      

  return my_sum;
}


double inner_heat_step(const unsigned long BS,
		       double * oldmu,
		       double * mu,
		       double * top, const bool wr_top,
		       double * bottom, const bool wr_bottom,
		       double * left, const bool wr_left,
		       double * right, const bool wr_right,
		       const bool check)
{
  typedef double (*rm)[BS];
  rm u = reinterpret_cast<rm>(mu);
  rm old_u = reinterpret_cast<rm>(oldmu);
  double my_sum = 0.0;
  double diff;
	
  // TOP LEFT
  double un = 0.25*( ( left[0] +
		       old_u[0][1] ) +
		     ( top[0] +
		       old_u[1][0] ) );
  if (check) {
    diff = un - old_u[0][0];
    my_sum += diff*diff;
  }
  u[0][0] = un;
	
  // TOP ROW, COL jj=1..BS-2
  for (unsigned long jj = 1; jj < BS-1; ++jj) {
    un = 0.25*( ( old_u[0][jj-1] +
		  old_u[0][jj+1] ) +
		( top[jj] +
		  old_u[1][jj] ) );
	  
    if (check) {
      diff = un - old_u[0][jj];
      my_sum += diff*diff;
    }
    u[0][jj] = un;
  }
	
  // TOP RIGHT
  un = 0.25*( ( old_u[0][BS-2] +
		right[0] ) +
	      ( top[BS-1] +
		old_u[1][BS-1] ) );
  if (check) {
    diff = un - old_u[0][BS-1];
    my_sum += diff*diff;
  }
  u[0][BS-1] = un;
	
  // ROW ii=1..BS-2
  for (unsigned long ii = 1; ii < BS-1; ++ii) {
	  
    // ROW ii, COL 0
    un = 0.25*( ( left[ii] +
		  old_u[ii][1] ) +
		( old_u[ii-1][0] +
		  old_u[ii+1][0] ) );
    if (check) {
      diff = un - old_u[ii][0];
      my_sum += diff*diff;
    }
    u[ii][0] = un;

    // ROW ii, COL jj=1..BS-2
    for (unsigned long jj = 1; jj < BS-1; ++jj) {
      // ROW ii, COL jj
      un = 0.25*( ( old_u[ii][jj-1] +
		    old_u[ii][jj+1] ) +
		  ( old_u[ii-1][jj] +
		    old_u[ii+1][jj] ) );
      if (check) {
	diff = un - old_u[ii][jj];
	my_sum += diff*diff;
      }
	    
      u[ii][jj] = un;
    }

    // ROW ii, COL BS-1
    un = 0.25*( ( old_u[ii][BS-2] +
		  right[ii] ) +
		( old_u[ii-1][BS-1] +
		  old_u[ii+1][BS-1] ) );
    if (check) {
      diff = un - old_u[ii][BS-1];
      my_sum += diff*diff;
    }
    u[ii][BS-1] = un;
  }
	
  // BOTTOM LEFT
  un = 0.25*( ( left[BS-1] +
		old_u[BS-1][1] ) +
	      ( old_u[BS-2][0] +
		bottom[0] ) );
  if (check) {
    diff = un - old_u[BS-1][0];
    my_sum += diff*diff;
  }
  u[BS-1][0] = un;

  // BOTTOM ROW, COL jj=1..BS-2
  for (unsigned long jj = 1; jj < BS-1; ++jj) {
    un = 0.25*( ( old_u[BS-1][jj-1] +
		  old_u[BS-1][jj+1] ) +
		( old_u[BS-2][jj] +
		  bottom[jj] ) );
	  
    if (check) {
      diff = un - old_u[BS-1][jj];
      my_sum += diff*diff;
    }
    u[BS-1][jj] = un;
  }
	
  // BOTTOM RIGHT
  un = 0.25*( ( old_u[BS-1][BS-2] +
		right[BS-1] ) +
	      ( old_u[BS-2][BS-1] +
		bottom[BS-1] ) );
  if (check) {
    diff = un - old_u[BS-1][BS-1];
    my_sum += diff*diff;
  }
  u[BS-1][BS-1] = un;

  // WRITE TO THE HALOS IF NEEDED
  if (wr_top)
    std::copy(&u[0][0], &u[0][BS], top);

  if (wr_bottom)
    std::copy(&u[BS-1][0], &u[BS-1][BS], bottom);

  if (wr_left)
    for (unsigned long i = 0; i < BS; ++i)
      left[i] = u[i][0];

  if (wr_right)
    for (unsigned long i = 0; i < BS; ++i)
      right[i] = u[i][BS-1];
      

  return my_sum;
}

#ifndef _HEAT_INNER_HH
#define _HEAT_INNER_HH

double inner_heat_step(const unsigned long BS,
		       double *mu,
		       double *top, const bool wr_top,
		       double *bottom, const bool wr_bottom,
		       double *left, const bool wr_left,
		       double *right, const bool wr_right,
		       const bool check);

double inner_heat_step(const unsigned long BS,
		       double *oldmu,
		       double *mu,
		       double *top, const bool wr_top,
		       double *bottom, const bool wr_bottom,
		       double *left, const bool wr_left,
		       double *right, const bool wr_right,
		       const bool check);

#endif

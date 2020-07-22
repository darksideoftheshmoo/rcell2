#ifndef _complex_
   #define _complex_
   #include "complex.h"
#endif
#ifndef _point_
   #define _point_
   #include "point.h"
#endif
struct complex *FFT_of_r_vs_theta(struct point *);
double FFT_ratio(struct point *);

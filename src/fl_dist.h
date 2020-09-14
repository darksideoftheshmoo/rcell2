#ifndef _point_
   #define _point_
   #include "point.h"
#endif
 
void maximum_pixels_within_fixed_radius(float *,int,int,
					struct point *,float,
					int **,int **, int *,struct point **);
void maximum_contiguous_pixels(float *,int,int,
			       struct point *,
			       int **,int **, int *);
void get_gauss_2d_parameters(float *,int,int,float,
			     struct point *,
			     float *,float *,float *,
			     int **,int **, int *);
void find_vacuole(struct point *, float *,int,int, int **,int **,int *);

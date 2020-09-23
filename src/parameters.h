//V1.2a The double variable were float in the previous version. 
//Using double for compatibility with GOptions library

extern double max_d_over_s_cut;
extern int max_pixels_per_cell;
extern int min_pixels_per_cell;
extern double background_reject_factor;
extern int nucleus_radii[6]; //V1.4.5 radii size initialization


extern double max_split_d_over_minor;
extern double I_over_U_for_match;

#ifndef _parameters_
#define _parameters_
double max_d_over_s_cut;
int max_pixels_per_cell;
int min_pixels_per_cell;
int nucleus_radii[6]; //V1.4.5 radii size initialization

double background_reject_factor;
double max_split_d_over_minor;
double I_over_U_for_match;
#endif

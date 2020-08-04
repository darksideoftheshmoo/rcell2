/*

Cell-ID is intended to identify cells in images and to calculate a
number of statistics, including statistics derived from corresponding
fluorescence images.
Copyright (C) 2005 Andrew Gordon

This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.

This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with this library; if not, write to the Free Software
Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

Andrew Gordon can be contacted at
agordon@molsci.org
Molecular Sciences Institute
2168 Shattuck Ave, 2nd Floor
Berkeley, CA  94704

**********************************************************************
Start-copyright-notice-for-libtiff
Libtiff software is used by Cell-ID for some of the reading in of
TIF image file data and also for creating new TIF files. Libtiff is
available at http://www.remotesensing.org/libtiff/. The libtiff software
was written by Sam Leffler while working for Silicon Graphics and
contains the following copyright notice:

   "Copyright (c) 1988-1997 Sam Leffler
    Copyright (c) 1991-1997 Silicon Graphics, Inc.

    Permission to use, copy, modify, distribute, and sell this software and
    its documentation for any purpose is hereby granted without fee, provided
    that (i) the above copyright notices and this permission notice appear in
    all copies of the software and related documentation, and (ii) the names
    of Sam Leffler and Silicon Graphics may not be used in any advertising or
    publicity relating to the software without the specific, prior written
    permission of Sam Leffler and Silicon Graphics.

    THE SOFTWARE IS PROVIDED "AS-IS" AND WITHOUT WARRANTY OF ANY KIND,
    EXPRESS, IMPLIED OR OTHERWISE, INCLUDING WITHOUT LIMITATION, ANY
    WARRANTY OF MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE.

    IN NO EVENT SHALL SAM LEFFLER OR SILICON GRAPHICS BE LIABLE FOR
    ANY SPECIAL, INCIDENTAL, INDIRECT OR CONSEQUENTIAL DAMAGES OF ANY KIND,
    OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS,
    WHETHER OR NOT ADVISED OF THE POSSIBILITY OF DAMAGE, AND ON ANY THEORY OF
    LIABILITY, ARISING OUT OF OR IN CONNECTION WITH THE USE OR PERFORMANCE
    OF THIS SOFTWARE."

End-copyright-notice-for-Libtiff
*********************************************




*/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "fit.h"
#include "tif_routines.h"
#include "segment.h"
#include "nums.h"
#include "fft_stats.h"
#include "image_type.h"
#include "split_and_overlap.h"
#include "contiguous.h"
#include "fft.h"
#include "fl_dist.h"
#include "parameters.h"
#include "flatten.h"

#define pi 3.1415926535
#define twopi 2.0*3.1415926535
#define nbins_for_cut_calculation 1000
// #define I_over_U_for_match 0.2--> now in parameters.txt
#define isect_max 1000
//#define min_pixels_per_cell 75
#define min_pixels_per_cell_hard 100

#define ccd_floor 128.0
#define ccd_floorp1 129.0

int xmax,ymax;
int xmax_ymax;

// declare locally thsese "extern" variables from header files, because the compiler failed
int nucleus_radii[6];              // will this help with: cellMagick.so: undefined symbol: nucleus_radii ??
int image_type;                    // it seems to help and be harmless: https://stackoverflow.com/questions/36040861/multiple-definitions-of-a-global-variable
double max_d_over_s_cut;
int max_pixels_per_cell;
int min_pixels_per_cell;
double max_split_d_over_minor;
double background_reject_factor;
int recalculate_internal;
double I_over_U_for_match;
int overall_id_offset;
int third_image_type;
// end declare locally

//Global Arrays
float *c=NULL;    //Bright field image array
float *fl=NULL;   //Fluorescence image array
float *third_image=NULL;
int *d=NULL;

int *work_array=NULL;

void sort_list_by_fl(struct point *);

void check_mem(void);

int total_overlap_of_all_cells(int,int);
struct point *point_malloc(void);
void point_free(struct point *);
void point_list_free(struct point *);
void point_list_adjust(struct point *,int,int);
struct point *fix_spirals(struct point *);
struct point *clean_up_tails(struct point *);
struct point *find_interior_points(struct point *);
float area_of_cell(struct point *);
float circumference_of_cell(struct point *);
struct point *circularize_points(struct point *, float);
void calculate_cut(void);
float cut,cut_high,cut_low,mid; //Cuts from gaussian fit to histogram
//of pixels.
float max_d_over_s(struct point *, struct point **, struct point **);
struct point *copy_cell_for_split_regions(struct point *,int);
void neighbor_interpolation(int,int,float,float,float *,float *,float *);
int get_median(int,int,int,float *,float *,float *,float *);
int get_median_deviations(int,int,float,float *,float *);
int combine_cells_in_cs_array(int,int,int);
int mark_distance_from_boundary(struct point *, int);

#define n_points_r_vs_theta 128
float r_vs_theta[n_points_r_vs_theta];
void r_vs_theta_from_boundary(struct point*);
struct point *boundary_from_r_vs_theta(void);
void fill_cos_sin_arrays(void);
double dtheta;
double *cos_theta=NULL;
double *sin_theta=NULL;
float centroid_x,centroid_y;
void statistics_from_r_vs_theta(int,
				float *,float *,float *,float *,
				float *,float *);

float smallest_circumference;

//Calculate the mode of the background for each of the time points
#define max_time 50000
float back_cur=0.0;
float back_pixels[max_time];
float back_pixels_1[max_time];
float back_pixels_2[max_time];
float total_fluorescence[max_time];
float total_area[max_time];

//Pointers to save the pixels which constitute the largest background.
//These pixels will be used to calculate the background level in the
//fluorescence image, and they'll be calculated from the BF image.
unsigned short int *pixels_for_background_x=NULL;
unsigned short int *pixels_for_background_y=NULL;
int n_pixels_for_background;

//Use over[] array to record where the cells are for overlap calculation.
//Put it up here so don't have to keep re-allocating the space.  Will
//calloc() space for it below, so it must be initialized to NULL otherwise
//space won't be allocated for it.
unsigned char *over=NULL;
int over_array_max=0;
void update_overlap_value(void);
void fill_overlap_array_with_point_list(struct point *);
void fill_overlap_array_with_point_list_offset(struct point *,int,int);
unsigned char overlap_value=0;
//Make maximum overlap_value the maximum for a signed short (ie,
// 7fff, ie, assume 2-bytes.  This should be by far large enough
//since at most it needs to be larger than the number of cells found
//in any given image (0x7fff=32767). (On the other hand, an
//unsigned char (max=255) would probably be too small sometimes.)
#define overlap_value_max 254

//We define a "cell" as a closed set of segments.
#define max_cells 200000
//struct point *vacuole[max_cells];
struct point *boundary[max_cells];  // "declare boundary as array 20000 of pointer to struct point" https://cdecl.org/?q=struct+point+*boundary%5B20000%5D%3B
struct point *interior[max_cells];
struct point *boundary_p1[max_cells];
struct point *interior_p1[max_cells];
struct point *boundary_m1[max_cells];
struct point *interior_m1[max_cells];
struct point *boundary_m2[max_cells];
struct point *interior_m2[max_cells];
struct point *boundary_m3[max_cells];
struct point *interior_m3[max_cells];

int save_realign_offset_x[max_cells];
int save_realign_offset_y[max_cells];
void re_align_cell(struct point *,int *,int *);

float calculate_volume_cone_old(struct point *);
int calculate_volume(struct point *,struct point *,
		     float *,float *,
		     float *,float *,
		     float *,float *,
		     float *,float *);
float integrate_gaussian_from_0_to_x(float,float);

double *gaussian_integral=NULL;
int gaussian_integral_nbins=5000;
double max_gaussian_integral=5.0;
float gaussian_integral_bin_width;

//For "local background correction"
struct point *boundary_p5[max_cells];
float back_p5[max_cells];
float pixels_back_p5[max_cells];
float pixels_total_p5[max_cells];
struct point *boundary_phalfminor[max_cells];
float back_phalfminor[max_cells];
float pixels_back_phalfminor[max_cells];
float pixels_total_phalfminor[max_cells];

float mean_x[max_cells],mean_y[max_cells];

float fft_stat[max_cells];
float vol_rotation[max_cells];
float vol_cone[max_cells];
float vol_sphere[max_cells];
float surface_area[max_cells];
float vol_eff_1[max_cells];
float vol_eff_2[max_cells];
float vol_eff_3[max_cells];
float vol_eff_4[max_cells];
float vol_eff_5[max_cells];
float vol_eff_6[max_cells];
float circumference[max_cells];
float major_axis1[max_cells],minor_axis1[max_cells];
float major_axis2[max_cells],minor_axis2[max_cells];
int n_points[max_cells];
float fluorescence[max_cells];
float fluorescence_p1[max_cells];
float fluorescence_m1[max_cells];
float fluorescence_m2[max_cells];
float fluorescence_m3[max_cells];
float med_fl[max_cells];
float mad_fl[max_cells];
float mcd_fl[max_cells];
float mcd_rms_fl[max_cells];
float vacuole_area[max_cells];
float vacuole_fl[max_cells];
float cell_area[max_cells];
float cell_area_p1[max_cells];
float cell_area_m1[max_cells];
float cell_area_m2[max_cells];
float cell_area_m3[max_cells];
int x_nucl[max_cells];
int y_nucl[max_cells];

#define nucleus_distribution_types 6 //V1.4.5 modified from 8 to 6
float fl_nucleus[max_cells][nucleus_distribution_types];
float area_nucleus[max_cells][nucleus_distribution_types];
float fl_nucleus_from_search[max_cells][nucleus_distribution_types];
struct point *p_nuclear_center[max_cells];
int fret_copy_type[max_cells];

//Some variables for the pixel-by-pixel comparison of current
//image to the previous
float *flprev=NULL;
float *cur_prev=NULL;;
int first_cur_prev_comparison=0; //A flag for first time
float diff_fl_prev(int,int); //For calculating offset

void calculate_global_stats_from_interior_and_boundary();

float pos_sig_mean_x[max_cells];
float pos_sig_mean_y[max_cells];
float pos_sig_mean_r[max_cells];
float pos_sig_rms_r[max_cells];
float neg_sig_mean_x[max_cells];
float neg_sig_mean_y[max_cells];
float neg_sig_mean_r[max_cells];
float neg_sig_rms_r[max_cells];

int location_in_cs_array[max_cells];
int pixel_list[max_cells];
int n_found=0;
int n_before_fret_copy=0;
//For connection between bottom and top of split images
float fret_mx,fret_my,fret_bx,fret_by;
int *fret_labels=NULL;
#define fret_offset 1000
#define lower_fret_region 1
#define higher_fret_region 2
int fret_region_use;
int have_fret_image=0;


#define n_fft_max 100

#define max_offsets 1000
int nucleus_offset_x[max_offsets];
int nucleus_offset_y[max_offsets];
int nucleus_n_offset=0;


struct blob {
  int index; //A unique number to label this cell list
  int flag; //for whatever use
  float x,y; //mean of interior points to cell
  float a; //Area of interior of cell
  int n; //Number of interior pixels to cell
  float fluor; //Sum of fluorescence from fluorescence picture
  float fluor_p1; //Sum of fluorescence from fluorescene picture
  float fluor_m1; //Sum of fluorescence from fluorescence picture
  float fluor_m2; //Sum of fluorescence from fluorescence picture
  float fluor_m3; //Sum of fluorescence from fluorescence picture
  float area_p1;
  float area_m1;
  float area_m2;
  float area_m3;
  float vacuole_area;
  float vacuole_fl;

  //Some local background stuff
  float back_p5;
  float pixels_back_p5;
  float pixels_total_p5;
  float back_phalfminor;
  float pixels_back_phalfminor;
  float pixels_total_phalfminor;

  int i_time; //Which loop was this found in
  int secs; //Time in seconds, starting at first image.  Time stamp is
  //take from metamorph's tag in the tiff file
  float i_over_u; //Intersection over Union with previous time point
  float fft_stat;
	int x_nucleus, y_nucleus; //V1.4.5
  float area_nucleus1;
  float area_nucleus2;
  float area_nucleus3;
  float area_nucleus4;
  float area_nucleus5;
  float area_nucleus6;
  float area_nucleus7;
  float area_nucleus8;
  float fl_nucleus1;
  float fl_nucleus2;
  float fl_nucleus3;
  float fl_nucleus4;
  float fl_nucleus5;
  float fl_nucleus6;
  float fl_nucleus7;
  float fl_nucleus8;
  float fl_nucleus_from_search1;
  float fl_nucleus_from_search2;
  float fl_nucleus_from_search3;
  float fl_nucleus_from_search4;
  float fl_nucleus_from_search5;
  float fl_nucleus_from_search6;
  float fl_nucleus_from_search7;
  float fl_nucleus_from_search8;

  //Some variables for the measure of the comparison of the current
  //fluorescence with the previous image.
  float pos_sig_mean_x;
  float pos_sig_mean_y;
  float pos_sig_mean_r;
  float pos_sig_rms_r;
  float neg_sig_mean_x;
  float neg_sig_mean_y;
  float neg_sig_mean_r;
  float neg_sig_rms_r;

  float vol_rotation;
  float vol_cone;
  float vol_sphere;
  float surface_area;
  float vol_eff_1;
  float vol_eff_2;
  float vol_eff_3;
  float vol_eff_4;
  float vol_eff_5;
  float vol_eff_6;
  float circumference;
  float major_axis_length1;
  float minor_axis_length1;
  float major_axis_length2;
  float minor_axis_length2;

  struct point *boundary;
  struct point *interior;
  struct blob *next;
  struct blob *prev;
};
struct blob *cs[max_cells];
int n_known=0;
int next_index=0;
int total_time=0;

struct point *make_boundary_list(int,int,int);
int neighboring_points(int, int, int);

//For my own data allocation system
//#define n_mem_blocks 1000
//Not using anymore (as of 5/22/03) so just make it 1
#define n_mem_blocks 1
struct point_mem{
  struct point p;
  struct point_mem *next;
};
struct point_mem *mem_blocks[n_mem_blocks];
struct point_mem *pmem_cur=NULL;
struct point_mem *pmem_start=NULL;
#define mem_size 10000
int next_block=0;

/*************************************************************/
int find_cells(struct point ***boundary_out,struct point ***interior_out){

  int i,j,k,l;

  float minor_s0,major_s0,minor_s1,major_s1;
  int di_split,dj_split;
  float tmp1,tmp2,r,r_max;
  int ix,iy;
  int isection,uon;
  float I_over_U,I_over_U_max;
  //double theta,ctheta,stheta;
  int n_small,n_norm,n_small_max;
  //int n_norm_max;

  int n_p;
  float tmp_x,tmp_y,tmp;
  float sum1_xy,sum_xx,sum_x,sum1_y,sum_n;
  float sum2_xy,sum2_y;
  int recalc_fret_offsets_loop;
  struct point *p;
  struct point *ptmp;
  struct point *ptmp2;
  struct point *s0;
  struct point *s1;
  int n_cells=0;

  int max_size;
  float cut_step;

  int n_fft_cur;
  int n_fft;
  int i_fft[n_fft_max];
  int j_fft[n_fft_max];
  double re_fft[n_fft_max];
  double im_fft[n_fft_max];

  struct complex *fft_in=NULL;
  struct complex *fft_out=NULL;
  double dtmp;
  double max1,cur_max;
  double min1;

  int smallest_white_area=75;

  struct contiguous_search *csearch;        /* Declare csearch pointer of type contiguous_search */
  struct contiguous_search csearch_memory;  /* Declare csearch_memory of type contiguous_search */
  //Some arrays for contiguous searches
  int *clist_x;
  int *clist_y;
  int n_contiguous=0;
  int list_cur=0;

  //For connection between bottom and top of split images
  int fret_offx,fret_offy;

  static int first=1;

  //float *ftmp;
  //float laplace[]=
  //  {
  //    0.0, 0.0, -1.0, 0.0, 0.0,
  //    0.0, -1.0, -2.0, -1.0, 0.0,
  //    -1.0, -2.0, 16.0, -2.0, -1.0,
  //    0.0, -1.0, -2.0, -1.0, 0.0,
  //    0.0, 0.0, -1.0, 0.0, 0.0
  //  };
  //int laplace_lims;
  //int laplace_size=5;
  //int ux,uy;
  //int jlim;
  //double r_dot_x,r_dot_y;
  //double rmax,rmin;

  int u;

  if (first==1){
    first=0;
    //Do some initialiations here even though it's not
    //strictly necessary. (Note that in C external and
    //static variables are supposed to be initialized to
    //0, but 1) 0 isn't necessarily NULL for pointers, although
    //in reality it probably always is, and 2) I may want
    //to make these boundary, interior, etc variables non-global
    //at some point, and weird crashes might ensue.
    for(u=0;u<max_cells;u++){
      boundary[u]=NULL;
      interior[u]=NULL;
      boundary_p1[u]=NULL;
      interior_p1[u]=NULL;
      boundary_m1[u]=NULL;
      interior_m1[u]=NULL;
      boundary_m2[u]=NULL;
      interior_m2[u]=NULL;
      boundary_m3[u]=NULL;
      interior_m3[u]=NULL;
      boundary_p5[u]=NULL;
      boundary_phalfminor[u]=NULL;
    }
  }

  *boundary_out=NULL; //to default to
  *interior_out=NULL;

  csearch=(&csearch_memory);  // csearch now points to csearch_memory, this coding is weird...
  //Define some arrays for the contiguous searching.
  (csearch->data_array)=c;
  (csearch->label_array)=d;
  (csearch->xmax)=xmax;      // assign xmax value to xmas property of the csearch struct
  (csearch->ymax)=ymax;      // assign ymax value to xmas property of the csearch struct

  n_pixels_for_background=0;


  if (image_type==bright_field){
    flatten_image(c,xmax,ymax,overwrite_image,x_and_y,linear);
  }
  if ((image_type==fret_bf_bottom_only)||
      (image_type==fret_bf_top_only)||
      (image_type==fret_bf_bottom_and_top)){

    have_fret_image=1;
    //x-projection only
    flatten_image(c,xmax,ymax,overwrite_image,x_only,nonlinear);
  }

  //Which region of BF to use. (If we're doing bottom and top, then
  //start with bottom.)
  fret_region_use=lower_fret_region;
  if (image_type==fret_bf_top_only){
    fret_region_use=higher_fret_region;
  }


  if (have_fret_image==1){

    if(fret_labels==NULL){
      printf("Haven't calculated upper and lower regions from");
      printf("fluorescence image (while removing edge points).\n");
      fflush(stdout);
      perror(0); //exit(0); //http://r-pkgs.had.co.nz/src.html
    }
    //Calculate offset from bottom to top of image
    fret_mx=0.0; //Starting values passed in
    fret_my=0.0;
    fret_bx=0.0;
    fret_by=-(float)(ymax/2);
    if (image_type!=fret_bf_bottom_and_top){
      calculate_split_offset(fret_labels,
			     fl,
			     &fret_mx,
			     &fret_bx,
			     &fret_my,
			     &fret_by,
			     xmax,ymax);
    }else{
      calculate_split_offset(fret_labels,
			     //c,
			     fl,
			     &fret_mx,
			     &fret_bx,
			     &fret_my,
			     &fret_by,
			     xmax,ymax);
    }

  }

  n_found=0;

  for(i=0;i<max_cells;i++){
    //fret_copy_type[] will be used to flag whether the cells we
    //end up with were found in the image, or were copied from the
    //lower or upper region to the opposite region.
    //0==>original, 1==>a copy.
    //Set to -999 here just to make sure that it gets set for every id.
    //(We'll check it below, not being set is bad.)
    fret_copy_type[i]=-999;
  }

  n_before_fret_copy=0;

  //Set a marker for doing fret_region_use=lower_ and then higher_fret_region
  //(Note n_found will continue to be incremented on second pass.
 fret_loop:

  if (have_fret_image==1){
    if (fret_region_use==lower_fret_region){
      printf("Searching lower part of split BF image for cells.\n");
    }else{
      printf("Searching higher part of split BF image for cells.\n");
    }
  }

  //Calculate cut to put on the data for the searching
  if((image_type==bright_field)||(image_type==confocal_transmission)
                               ||(have_fret_image==1)){
    calculate_cut();
    smallest_white_area=75;
    smallest_circumference=(float)sqrt(((double)smallest_white_area)/3.14159)
      *2.0*3.14159;
  }else if(image_type==metamorph_deconvolution){
    //Deconvolution case, calculate cut differently
    //We're going to vary the low cut and keep track of the number of very
    //small contiguous regions and the number of normal ones.  We're going to
    //find the cut values that maximize the number of very small ones and then
    //maximizes  the ratio of the small to normal ones.  (The number of small
    //grows as  the giant background regions start to break up.)  The number
    //of normal ones starts to decline as we cut into the good regions.
    cut_step=10.0;
    cut_low=0.5;
    cut_high=1.0e30;
    n_small_max=0;
    r_max=0.0;
    //The background are the points between cut_low and cut_high
    csearch->cut_behavior=cut_between;
    (csearch->p)=NULL; //No point list means loop over entire array

  decon_start:

    csearch->cut_low=cut_low;
    csearch->cut_high=cut_high;
    do_contiguous_search(csearch);
    //Now search results
    n_small=0;
    n_norm=0;
    for(i=0;i<(csearch->n_lists_found);i++){
      n_contiguous=(csearch->npoints_in_list)[i];
      if((n_contiguous>150)&&(n_contiguous<1000)){
	      n_norm++;
      }else if((n_contiguous>0)&&(n_contiguous<20)){
	      n_small++;
      }
    }
    printf("For cut_low=%e nsmall=%i, nnorm=%i.\n",
	   cut_low,n_small,n_norm);fflush(stdout);
    if (n_small>n_small_max){
      tmp1=cut_low;
      n_small_max=n_small;
      r_max=0.0; //don't consider max ratio until after we find the n_small max
      if (image_type==hexagonal_grid){
	      cut_step=1.0;
      }
    } else if (n_small>0){
      r=((float)n_norm)/((float)n_small);
      if (r>r_max){
	      tmp2=cut_low;
	      r_max=r;
      }
    }else{
      r=-1.0;
    }

    if ((r>=r_max)||(n_small>=n_small_max)){
      cut_low+=cut_step;
      goto decon_start;
    }

    if(image_type==metamorph_deconvolution){
      cut_high=(tmp1+tmp2)/2.0;
      cut_low=cut_high+5.0;
      printf("Decon cut values are low=%e and high=%e\n",cut_low,cut_high);
    }else{
      cut_high=tmp1;
      cut_low=cut_high+5.0;
      printf("Array data cut vals are low=%e and high=%e\n",cut_low,cut_high);
    }

    //V1.2a TODO hard coded!
    smallest_white_area=75;
    smallest_circumference=(float)sqrt(((double)smallest_white_area)/3.14159)
      *2.0*3.14159;

  }else if(image_type==hexagonal_grid){
    //We have a giant hexagonal grid.  Use a fourier transform to pick
    //it out.

    free(fft_in);
    free(fft_out);
    fft_in=(struct complex *)malloc(xmax_ymax*sizeof(struct complex));
    for(k=0;i<xmax_ymax;k++){
      fft_in[k].r=(double)(c[k]);
      fft_in[k].i=0.0;
    }

    fft_out=FFT_2d(fft_in,xmax,1);

    //Find top six locations
    max1=1.0e50;
    cur_max=-1.0;
    n_fft=6;
    n_fft_cur=0;
    l=0;
    while(l<n_fft){
      for(i=0;i<xmax;i++){
	      for(j=0;j<ymax;j++){
	        if((i==0)&&(j==0)) continue;  //Skip the 0-frequency part.
	        k=(j*xmax+i);
	        dtmp=(fft_out[k].r)*(fft_out[k].r)+(fft_out[k].i)*(fft_out[k].i);
	        if(dtmp>=max1) continue;
	        if(dtmp>cur_max){
	          cur_max=dtmp;
	          i_fft[n_fft_cur]=i;
	          j_fft[n_fft_cur]=j;
	          re_fft[n_fft_cur]=fft_out[k].r;
	          im_fft[n_fft_cur]=fft_out[k].i;
	        }
	      }
      }
      max1=cur_max;
      cur_max=-1.0;
      n_fft_cur++;

      l++;
    }
    n_fft=n_fft_cur;

    //Zero everything but the top n_fft points
    for(k=0;k<(xmax_ymax);k++){
      fft_in[k].r=0.0;
      fft_in[k].i=0.0;
    }
    for(l=0;l<n_fft;l++){
      i=i_fft[l];
      j=j_fft[l];
      k=(j*xmax+i);
      fft_in[k].r=re_fft[l];
      fft_in[k].i=im_fft[l];
      printf("(%i,%i)=(%e,%e)\n",i,j,fft_in[k].r,fft_in[k].i);
    }

    free(fft_out); //To prevent memory leak
    fft_out=FFT_2d(fft_in,xmax,-1);

    //Now reset c[] to be the real part of the fixed image.
    max1=-1.0e50;
    min1=1.0e50;
    for(i=0;i<xmax;i++){
      for(j=0;j<ymax;j++){
	      k=(j*xmax+i);
	      dtmp=fft_out[k].r;
	      c[k]=(float)(dtmp);
	      if(dtmp>max1) max1=dtmp;
	      if(dtmp<min1) min1=dtmp;
      }
    }
    //Subtract off minimum so c[] is always >=0
    tmp1=(float)min1;
    for(i=0;i<xmax_ymax;i++){
	    c[i]-=tmp1;
    }

    cut_high=0.4*((float)(max1-min1));
    cut_low=cut_high+5.0;
    smallest_white_area=0;
    smallest_circumference=0.0;
    max_d_over_s_cut=6.0e30; //Remove cut basically
    printf("Cut=%e\n",cut_high);

  }else if(image_type==membrane_tagged_fl){
    //Nothing yet....
    //TODO
  } //end of cut calculation for different image types

  printf("BF-cut-low=%e, BF-cut-high=%e\n",cut_low,cut_high);fflush(stdout);
  //Now, the cells are surrounded by a dark region and inside the cells is
  //higher than the background.  Divide the pixels into high pixels, low
  //pixels and middle pixels (background).
  for(i=0;i<xmax_ymax;i++){
    if(c[i]>cut_high){      //inside cells
      d[i]=cell_in;
    }else if(c[i]>cut_low){ //background
      d[i]=cell_out;
    }else{
      d[i]=cell_border;
    }
  }

  //Change the behavior of the contiguous search to use the
  //integer array d[]

  (csearch->label_array)=d;
  (csearch->cut_behavior)=equal_to_labels;
  (csearch->p)=NULL; //No point list means use entire image

  //If we have a deconvolution image then the border regions might be set
  //to 0.  If this is the case, then we don't want "half-cells" that get
  //cut off by the border region.  To prevent the program from picking them
  //up, we set the border regions to be "cell_in" and then find all the
  //contiguous pixels starting at a point in the border, and then removing
  //all of those from the image.  This will pick up the cells that are partly
  //in the border region and then remove them.
  if((image_type==bright_field)||(image_type==confocal_transmission)
                               ||(have_fret_image==1)){
    //V1.2a TODO hard coded!
    max_size=1000; //Normal case, inside cells has some noise, use this
    //to set some noise regions to be "inside cells"
  }else if(image_type==metamorph_deconvolution){ //deconvolution case
    max_size=1000;

    for(k=0;k<xmax;k++){
      j=xmax-1-k;
      for(i=k;i<=j;i++){ //Proceed around border in decreasing squares
	      l=0;
	      u=(i*xmax+k);
	      if (c[u]<0.5){
	        d[u]=cell_in;
	      }else{
	        l++; //count how many of the four cases failed.
	      }
	      u=(i*xmax+j);
	      if (c[u]<0.5){
	        d[u]=cell_in;
	      }else{
	        l++; //count how many of the four cases failed.
	      }
	      u=(k*xmax+i);
	      if (c[u]<0.5){
	        d[u]=cell_in;
	      }else{
	        l++; //count how many of the four cases failed.
	      }
	      u=(j*xmax+i);
	      if (c[u]<0.5){
	        d[u]=cell_in;
	      }else{
	        l++; //count how many of the four cases failed.
	      }
	      if (l==4){ //All four are non-zero
	        goto done_square_loop; //Stop at first zero pixel
	      }
      }
    }

    done_square_loop:

    if (k>0) {
      //Find all contiguous "cell_in" pixels starting in this border region.
      (csearch->p)=point_malloc();
      (csearch->p)->i=0;
      (csearch->p)->j=0;
      (csearch->p)->prev=NULL;
      (csearch->p)->next=NULL;
      (csearch->label_value)=cell_in;
      (csearch->label_value)=cell_in;
      do_contiguous_search(csearch);
      //Remove all of found pixels since they overlap the border.
      //printf("n,k=%i,%i\n",list_cur,n_contiguous);
      if ((csearch->n_lists_found)>0){
	     list_cur=(csearch->list_start)[0];
	     n_contiguous=(csearch->npoints_in_list)[0];
	     clist_x=(csearch->list_found_x);
	     clist_y=(csearch->list_found_y);
	     for(k=list_cur;k<(list_cur+n_contiguous);k++){
	       d[(clist_y[k]*xmax)+(clist_x[k])]=pixel_removed;
	     }
      }
      point_free(csearch->p);
    }

  }else{
    max_size=100000;
  }


  if ((image_type==bright_field)||(image_type==confocal_transmission)
                                ||(have_fret_image==1)){
    //Change background groups that are too small to cell_in

    (csearch->p)=NULL; //Use entire array
    (csearch->cut_behavior)=equal_to_labels;
    (csearch->label_array)=d;
    (csearch->label_value)=cell_out;
    do_contiguous_search(csearch);
    (clist_x)=(csearch->list_found_x);
    (clist_y)=(csearch->list_found_y);
    for(i=0;i<(csearch->n_lists_found);i++){
      n_contiguous=(csearch->npoints_in_list)[i];
      list_cur=(csearch->list_start)[i];
      if(n_contiguous>max_size){ //These guys are solidly background, so
	                               //remove them from image
	    //First check if they're so big that we should average their
	    //fluorescence level to get the background.
	    //if (n_contiguous>n_pixels_for_background){ //Save for below
	    //  free(pixels_for_background_x);
	    //  free(pixels_for_background_y);
	    //  pixels_for_background_x=(unsigned short int *)
	    //    malloc(n_contiguous*sizeof(unsigned short int));
	    //  pixels_for_background_y=(unsigned short int *)
	    //    malloc(n_contiguous*sizeof(unsigned short int));
	    //  n_pixels_for_background=0;
	    //  for(k=list_cur;k<(list_cur+n_contiguous);k++){
	    //    pixels_for_background_x[n_pixels_for_background]=clist_x[k];
	    //    pixels_for_background_y[n_pixels_for_background]=clist_y[k];
	    //    n_pixels_for_background++;
	    //  }
        //}
	      //Now remove these pixels
	      for(k=list_cur;k<(list_cur+n_contiguous);k++){
	        d[(clist_y[k]*xmax)+(clist_x[k])]=pixel_removed;
	      }
      }else{ //Otherwise set them to cell_in to correct interior of
	           //cells for background
        for(k=list_cur;k<(list_cur+n_contiguous);k++){
	        d[(clist_y[k]*xmax)+(clist_x[k])]=cell_in;
        }
      }
    }
  }

  //Change cell_border groups that are too small to cell_in
  if(image_type!=hexagonal_grid){

    (csearch->p)=NULL; //Use entire array
    (csearch->cut_behavior)=equal_to_labels;
    (csearch->label_array)=d;
    (csearch->label_value)=cell_border;
    do_contiguous_search(csearch);
    (clist_x)=(csearch->list_found_x);
    (clist_y)=(csearch->list_found_y);
    for(i=0;i<(csearch->n_lists_found);i++){
      n_contiguous=(csearch->npoints_in_list)[i];
      list_cur=(csearch->list_start)[i];
      if(n_contiguous<50){ //Set them to cell_in
	      for(k=list_cur;k<(list_cur+n_contiguous);k++){
	        d[(clist_y[k]*xmax)+(clist_x[k])]=cell_in;
	      }
      }
    }

  }


  //Remove groups of white that are too small, and give the others
  //their own label
  (csearch->p)=NULL; //Use entire array
  (csearch->cut_behavior)=equal_to_labels;
  (csearch->label_array)=d;
  (csearch->label_value)=cell_in;
  do_contiguous_search(csearch);
  (clist_x)=(csearch->list_found_x);
  (clist_y)=(csearch->list_found_y);
  for(i=0;i<(csearch->n_lists_found);i++){
    n_contiguous=(csearch->npoints_in_list)[i];
    list_cur=(csearch->list_start)[i];
    if(n_contiguous<smallest_white_area){ //Remove points
      for(k=list_cur;k<(list_cur+n_contiguous);k++){
	      d[(clist_y[k]*xmax)+(clist_x[k])]=pixel_removed;
      }
    }else{ //Save these guys with their own label
      n_p=100+n_cells;
      for(k=list_cur;k<(list_cur+n_contiguous);k++){
	      d[(clist_y[k]*xmax)+(clist_x[k])]=n_p;
      }
      pixel_list[n_cells]=list_cur;
      n_cells++;
    }
  }

  //Treat each contiguous list of white cells as a potential cell.
  //We have the start of each location in the clist_x,clist_y arrays
  //in pixel_list[n_cells].  Now take each blob and make a border.
  //Set work_array[] to 0 since make_boundary_list() uses it.
  memset(work_array,0,(xmax_ymax)*sizeof(int));

  for(i=0;i<n_cells;i++){
    j=pixel_list[i];
    p=make_boundary_list(clist_x[j],clist_y[j],100+i);

    //Remove cells that have a border point on edge of image
    j=0;
    for(ptmp=p;ptmp!=NULL;ptmp=ptmp->next){
      if((ptmp->i==0)||(ptmp->i==(xmax-1))||(ptmp->j==0)||(ptmp->j==(ymax-1))){
	      j=1;
	      break;
      }
      //Do same for the split image edges
      if(have_fret_image==1){
	      u=((ptmp->j)+1)*xmax+(ptmp->i);
	      if (u>=xmax_ymax){
	        j=1;
	        break;
	      }else{
	        if (fret_labels[u]!=fret_region_use){
	          j=1;
	          break;
	        }
	      }
	      u=((ptmp->j)-1)*xmax+(ptmp->i);
	      if (u<0){
	        j=1;
	        break;
	      }else{
	        if (fret_labels[u]!=fret_region_use){
	          j=1;
	          break;
	        }
	      }
	      u=(ptmp->j)*xmax+((ptmp->i)+1);
	      if (u>=xmax_ymax){
	        j=1;
	        break;
	      }else{
	        if (fret_labels[u]!=fret_region_use){
	          j=1;
	          break;
	        }
	      }
	      u=(ptmp->j)*xmax+((ptmp->i)-1);
	      if (u<0){
	        j=1;
	        break;
	      }else{
	        if (fret_labels[u]!=fret_region_use){
	          j=1;
	          break;
	        }
	      }
      }
    }

    if(j==0){
      if(p!=NULL){
	      boundary[n_found]=p;
	      n_found++;
      }
    }else{ //free points
      for(ptmp=p;ptmp!=NULL;ptmp=ptmp2){
	      ptmp2=ptmp->next;
	      point_free(ptmp);
      }
    }

  }

  //Calculate maximum of s/d where s is length along circumference
  //and d is distance between any two points.  If this is very high, it
  //means cell is pinched and we probably should split it into two cells.
  for(i=0;i<n_found;i++){
    tmp=max_d_over_s(boundary[i],&s0,&s1);
    if(tmp<=0.0){ //Remove the points entirely (this happens if
      //circumference is too small)
      point_list_free(boundary[i]);
      boundary[i]=boundary[n_found-1]; //Put last in list here
      n_found--;
      i--; //Check this point again
    }else if(tmp>max_d_over_s_cut){
      //s0 and s1 point to the two points for whom d/s is maximum
      //Split the boundary at this point and make two new cells
      //s0 is the earlier in the list and s1 the later
      //To make the split, first remove the NULL at the end of the list
      //and the NULL at the beginning of the list
      for(p=boundary[i];(p->next)!=NULL;p=p->next) ;
      p->next=boundary[i]; //p is end of the list
      boundary[i]->prev=p; //boundary[i] is beginning

      //Put NULL's at the end of the new lists
      (s0->prev)->next=NULL;
      (s1->prev)->next=NULL;
      //And the beginnings also
      s0->prev=NULL;
      s1->prev=NULL;

      //Pause here to check that we really want to split this cell.
      //compare the distance between the two points we're going to
      //split it at, with the minor axis of the two new cells.
      r_vs_theta_from_boundary(s0);
      statistics_from_r_vs_theta(0,&tmp,&major_s0,&minor_s0,
				 &tmp,&tmp,&tmp);
      r_vs_theta_from_boundary(s1);
      statistics_from_r_vs_theta(0,&tmp,&major_s1,&minor_s1,
				 &tmp,&tmp,&tmp);
      di_split=((s0->i)-(s1->i));
      dj_split=((s0->j)-(s1->j));
      tmp=((float)sqrt((double)(di_split*di_split+dj_split*dj_split)));
      if (
	      ((tmp/minor_s0)<max_split_d_over_minor)&&
	      ((tmp/minor_s1)<max_split_d_over_minor)){

        //Go ahead with split
	      boundary[i]=s0; //Overwrite current position with first cell
	      boundary[n_found]=s1; //Put second at the end of the list
	      n_found++;
	      i--; //Repeat the check on the first cell
      }else{
	      printf("Undoing split: %i (%e,%e,%e)\n",i,tmp,minor_s0,minor_s1);
	      //Undo the mess we made above. (Recombine the two cells)
	      for(p=s0;(p->next)!=NULL;p=(p->next)) ; //End of s0
	      //Heal the split here
	      (p->next)=s1;
	      (s1->prev)=p;
	      boundary[i]=s0; //Doesn't matter where we start it
      }
    }
  }

  //Check for bad fft's. (Note that I calculate fft again below.
  //The reason is that when I remove the cells I have to remove all
  //the global arrays associated with them. (I really should have all
  //of those compacted into a structure.) But if I remove them here,
  //I don't have to do that since no data has yet been calculated.
  for(i=0;i<n_found;i++){
    tmp=FFT_ratio(boundary[i]);
    if(tmp<=0.0){ //Remove the points entirely (this happens
      //if can't get a good calculation of fft. It implies
      //cells are convex I think. In the data it looks to happen
      //for spurious cells that look like "halos" around other cells.)
      point_list_free(boundary[i]);
      boundary[i]=boundary[n_found-1]; //Put last in list here
      n_found--;
      i--; //Check this point again
    }
  }


  //Calculate interior of each of these and remove cells that
  //overlap too much (We have to calculate interior[] lists up here
  //since remove_overlaps() uses them. It also uses the number of points
  //in the array n_points[], so we have to calculate that here also.
  for(i=n_before_fret_copy;i<n_found;i++){
    interior[i]=find_interior_points(boundary[i]);
    n_p=0;
    for(p=interior[i];p!=NULL;p=(p->next))n_p++;
    //Remove cells that are too large or too small
    if ((n_p>max_pixels_per_cell)||(n_p<min_pixels_per_cell)){
      //Free all points
      point_list_free(interior[i]);
      point_list_free(boundary[i]);
      boundary[i]=boundary[n_found-1];
      n_found--;
      i--; //Check i again
    }else{
      n_points[i]=n_p;
    }
  }

  remove_overlaps();

  //Take care of FRET stuff comparing bottom and top, etc.
  if ((image_type==fret_bf_bottom_only)||
      (image_type==fret_bf_top_only)){
    //Only half the image has data
    //For this case, we simply copy all the cells we found
    //in the bottom or top part of the image to the opposte part.

    //Run over all the found boundaries and make a copy of
    //them in the upper image, and increment n_found.
    //Note that boundary[i] is such that the first n_before_fret_copy
    //are for id<fret_offset, and the fret_offset+ cells are all
    //labeled _after_ the id<fret_offset cells.
    n_before_fret_copy=n_found;
    for(i=0;i<n_before_fret_copy;i++){
      if (image_type==fret_bf_bottom_only){ //bottom to top
	boundary[n_found]=copy_cell_for_split_regions(boundary[i],0);
      }else{ //top to bottom
	boundary[n_found]=copy_cell_for_split_regions(boundary[i],1);
      }
      interior[n_found]=find_interior_points(boundary[n_found]);
      //Number of interior points might be slightly different than
      //n_points[i].
      n_p=0;
      for(p=interior[n_found];p!=NULL;p=(p->next))n_p++;
      n_points[n_found]=n_p;

      fret_copy_type[i]=0; //0 means not a copy
      fret_copy_type[n_found]=1; //1 means is a copy

      n_found++; //Just made a copy
    }
  }else if (image_type==fret_bf_bottom_and_top){
    //Continue fret split image stuff if we're searching both top and
    //bottom separately for cells
    if (fret_region_use==lower_fret_region){
      fret_region_use=higher_fret_region;
      n_before_fret_copy=n_found; //To mark how many already found
      printf("Number new cells found in lower region: %i\n",n_found);
      goto fret_loop;
    }
    //If we're here then we've already done two passes. We re-sort the
    //cells so that the bottom and top cells match. Ie, for the cells
    //0<=i<n_before_copy, the match should be cell i+n_before_fret_copy.
    //This will mark the cells as matches for later stuff.

    printf("Number new cells found in upper region: %i\n",
	   n_found-n_before_fret_copy);

    for(recalc_fret_offsets_loop=0;recalc_fret_offsets_loop<1;
	recalc_fret_offsets_loop++){
      //First find all the matches and then use those to re-define
      //fret_bx,fret_mx, etc. The do it again.
      //(added 10/22/03)
      sum1_xy=0.0;sum_xx=0.0;sum_x=0.0;sum1_y=0.0;sum_n=0.0;
      sum2_xy=0.0;sum2_y=0.0;
      //Loop over all found lower cells
      for(i=0;i<n_before_fret_copy;i++){
	//Fill the overlap array in its upper region with cell data from
	//the lower region so can check for lower/upper overlaps.
	update_overlap_value();
	for(p=interior[i];p!=NULL;p=(p->next)){
	  tmp=((float)(p->i));
	  fret_offx=((int)(fret_bx+fret_mx*tmp+0.5));
	  fret_offy=((int)(fret_by+fret_my*tmp+0.5));
	  ix=(p->i)+fret_offx; //Note it may not be in image!
	  iy=(p->j)+fret_offy;
	  if ((ix>=0)&&(ix<xmax)&&(iy>=0)&&(iy<ymax)){
	    u=(iy*xmax+ix);
	    over[u]=overlap_value;
	  }
	}
	//Now find best overlap of all the higher cells (note above we can't
	//just do overlap(interior[j],0,joffset) since the offset from
	//bottom to the top isn't just a constant.
	k=-1;
	I_over_U_max=-9999.0;
	for(j=(n_before_fret_copy+i);j<n_found;j++){
	  isection=overlap(interior[j],0,0);//Compare jth with ith point
	  uon=n_points[i]+n_points[j]-isection; //union of points
	  if (uon>0){ //just in case
	    I_over_U=((float)isection)/((float)uon);
	  }else{
	    I_over_U=1.0;
	  }
	  if (I_over_U>I_over_U_max){
	    I_over_U_max=I_over_U;
	    k=j;
	  }
	}
	if(I_over_U_max>=I_over_U_for_match){ //Found a match
	  //Make dx and dy from the means of the two cells.
	  //(Unfortunately haven't calculated that yet, so waste time
	  //calculating again here)
	  tmp1=0.0;
	  tmp2=0.0;
	  tmp_x=0.0;
	  tmp_y=0.0;
	  tmp=0.0;
	  for(p=interior[i];p!=NULL;p=(p->next)){
	    tmp1+=((float)(p->i));
	    tmp2+=((float)(p->j));
	    tmp+=1.0;
	  }
	  if (tmp<=2.0){
	    goto skip_new_offset_function;
	  }
	  tmp1/=tmp;
	  tmp2/=tmp;
	  tmp=0.0;
	  for(p=interior[k];p!=NULL;p=(p->next)){
	    tmp_x+=((float)(p->i));
	    tmp_y+=((float)(p->j));
	    tmp+=1.0;
	  }
	  if (tmp<=2.0){
	    goto skip_new_offset_function;
	  }
	  tmp_x/=tmp;
	  tmp_y/=tmp;

	  //Accumulate dx,dy statistics to calculate new slope
	  //and intercept from these.
	  tmp_x-=tmp1;
	  tmp_y-=tmp2; //Vector from lower to upper region
	  //For fitting del(x)=mx+b and del(y)=mx+b
	  sum1_xy+=(tmp_x*tmp1);
	  sum1_y+=(tmp_x);
	  sum2_xy+=(tmp_y*tmp1);
	  sum2_y+=(tmp_y);
	  sum_xx+=(tmp1*tmp1);
	  sum_x+=(tmp1);
	  sum_n+=1.0;
	}
      }
      //Now use the sum's to calculate a new linear fit to the matching
      //found cells.
      if (sum_n<=3.0){
	goto skip_new_offset_function;
      }
      sum1_xy/=sum_n;
      sum1_y/=sum_n;
      sum2_xy/=sum_n;
      sum2_y/=sum_n;
      sum_xx/=sum_n;
      sum_x/=sum_n;

      tmp=sum_xx-sum_x*sum_x;
      if (tmp>0.0){
	printf("----->%e %e ",fret_mx,fret_bx);
	fret_mx=(sum1_xy-sum_x*sum1_y)/tmp;
	fret_bx=(-sum1_xy*sum_x+sum1_y*sum_xx)/tmp;
	printf("----->%e %e\n",fret_mx,fret_bx);
	printf("----->%e %e ",fret_my,fret_by);
	fret_my=(sum2_xy-sum_x*sum2_y)/tmp;
	fret_by=(-sum2_xy*sum_x+sum2_y*sum_xx)/tmp;
	printf("----->%e %e\n",fret_my,fret_by);
      }
    }
    fflush(stdout);

    skip_new_offset_function:



    //Loop over all found lower cells
    for(i=0;i<n_before_fret_copy;i++){
      //Fill the overlap array in its upper region with cell data from
      //the lower region so can check for lower/upper overlaps.
      update_overlap_value();
      for(p=interior[i];p!=NULL;p=(p->next)){
	tmp=((float)(p->i));
	fret_offx=((int)(fret_bx+fret_mx*tmp+0.5));
	fret_offy=((int)(fret_by+fret_my*tmp+0.5));
	ix=(p->i)+fret_offx; //Note it may not be in image!
	iy=(p->j)+fret_offy;
	if ((ix>=0)&&(ix<xmax)&&(iy>=0)&&(iy<ymax)){
	  u=(iy*xmax+ix);
	  over[u]=overlap_value;
	}
      }
      //Now find best overlap of all the higher cells (note above we can't
      //just do overlap(interior[j],0,joffset) since the offset from
      //bottom to the top isn't just a constant.
      k=-1;
      I_over_U_max=-9999.0;
      for(j=(n_before_fret_copy+i);j<n_found;j++){
	isection=overlap(interior[j],0,0);//Compare jth with ith point
	uon=n_points[i]+n_points[j]-isection; //union of points
	if (uon>0){ //just in case
	  I_over_U=((float)isection)/((float)uon);
	}else{
	  I_over_U=1.0;
	}
	if (I_over_U>I_over_U_max){
	  I_over_U_max=I_over_U;
	  k=j;
	}
      }
      if (I_over_U_max<I_over_U_for_match){ //No good match

	if (I_over_U_max>0.0){
	  //We're going to delete the smaller of the two cells, and then
	  //copy the better one onto it with copy_cell_...()
	  //j will be the cell we delete and then overwrite below:
	  if (n_points[i]<n_points[k]){
	    j=i;
	  }else{
	    j=k;
	  }
	  //Free the points for cell j
	  point_list_free(interior[j]);
	  point_list_free(boundary[j]);
	}else{ //If I_over_U_max==0.0 (or is -9999.0 which can happen if
	  //there aren't any upper cells), then we didn't find anything
	  //at all. In this case, make a new cell for the upper region
	  //and copy the lower onto it,
	  //but first do an extra quality cut since maybe the cell
	  //wasn't found in the upper region because it's a bad cell
	  if (n_points[i]<min_pixels_per_cell_hard){ //Remove the cell
	    point_list_free(interior[i]);
	    point_list_free(boundary[i]);
	    //Overwrite this point (note we haven't yet looked for
	    //a match at n_before_fret_copy-1)
	    if (i<(n_before_fret_copy-1)){
	      boundary[i]=boundary[n_before_fret_copy-1];
	      interior[i]=interior[n_before_fret_copy-1];
	      n_points[i]=n_points[n_before_fret_copy-1];
	      fret_copy_type[i]=fret_copy_type[n_before_fret_copy-1];
	    }
	    n_before_fret_copy--;
	    //Move all the upper region points down by 1 since
	    //they should all be offset from i by n_before_fret_copy, but
	    //we just decremented n_before_fret_copy;
	    for(j=n_before_fret_copy;j<(n_found-1);j++){
	      boundary[j]=boundary[j+1];
	      interior[j]=interior[j+1];
	      n_points[j]=n_points[j+1];
	      fret_copy_type[j]=fret_copy_type[j+1];
	    }
	    n_found--;
	    i--; //Must check the current value of i in loop
	    continue; //Don't do anything more for cell number i
	  }else{ //no match at all, but large enough to consider a cell
	    j=n_found; //Will be a new cell
	    k=n_found; //Best match to the cell for below
	    n_found++;
	  }
	}
	if (j>=n_before_fret_copy){
	  //The cell copy will be from i (which is lower) to j (which
	  //is upper)
	  boundary[j]=copy_cell_for_split_regions(boundary[i],0);
	  fret_copy_type[i]=0; //This cell was the original
	  fret_copy_type[j]=1; //This cell was a copy.
	}else{
	  //The cell copy will be from k (not i) (which is upper) to
	  //j (which is lower)
	  boundary[j]=copy_cell_for_split_regions(boundary[k],1);
	  fret_copy_type[k]=0; //This cell was the original
	  fret_copy_type[j]=1; //This cell was a copy.
	}
	interior[j]=find_interior_points(boundary[j]);
	//Number of interior points might be slightly different than
	//n_points[i or k].
	n_p=0;
	for(p=interior[j];p!=NULL;p=(p->next))n_p++;
	n_points[j]=n_p;
      }else{
	//We have a good match between lower and upper (i and k)
	//Mark both cells as non-copies. Below we'll re-arrange
	//the upper-cell location so it's offset by n_before_fret_copy
	//from the lower.
	fret_copy_type[i]=0;
	fret_copy_type[k]=0;
      }

      //End of check of if we don't have a match. If we're here then
      //either we have a match, or we didn't have one, but made one
      //by copying a cell and setting k appropriately.
      //Now make it so indices are such that
      //j=i+n_before_fret_copy (for both cases that we have a match
      //and that we don't have a match, but that we copied the one
      //regions's cell onto the others.
      j=i+n_before_fret_copy;
      ptmp=boundary[k];
      ptmp2=interior[k];
      n_p=n_points[k];
      ix=fret_copy_type[k];
      boundary[k]=boundary[j];
      interior[k]=interior[j];
      n_points[k]=n_points[j];
      fret_copy_type[k]=fret_copy_type[j];
      boundary[j]=ptmp;
      interior[j]=ptmp2;
      n_points[j]=n_p;
      fret_copy_type[j]=ix;
    }

    //We're now done the loop over matching the lower cells to the upper.
    //For each of the lower cells we started with, we either matched it
    //to an upper (by deleting upper and copying in some cases) or
    //added a new upper cell. In either case, we now have that the number
    //of upper cells is >= number of lower cells.
    //The excess (unmatched) upper cells all show up between
    //(n_before_fret_copy*2) and n_found since n_before_fret_copy is the
    //number of lower cells and each of those now has a match whose index
    //is i+n_before_fret_copy.
    for(i=(2*n_before_fret_copy);i<n_found;i++){
      //looping over unmatched upper cells
      //Do an extra quality cut since maybe the cell
      //wasn't found in the lower region because it's a bad cell
      if (n_points[i]<min_pixels_per_cell_hard){ //Remove the cell

	point_list_free(interior[i]);
	point_list_free(boundary[i]);

	//Overwrite this point (note we haven't yet looked for
	//a match at n_found-1) with last point
	boundary[i]=boundary[n_found-1];
	interior[i]=interior[n_found-1];
	n_points[i]=n_points[n_found-1];
	fret_copy_type[i]=fret_copy_type[n_found-1];
	i--;
	n_found--;
      }else{ //We're going to consider current a good cell and add a copy
	//to the lower regions
	//First move all the upper cells up one to make space for new
	//lower cell
	for (j=(n_found-1);j>=n_before_fret_copy;j--){
	  boundary[j+1]=boundary[j];
	  interior[j+1]=interior[j];
	  n_points[j+1]=n_points[j];
	  fret_copy_type[j+1]=fret_copy_type[j];
	}
	i++; //Must increment since we moved all cells up one
	//Put the new lower cell at the location n_before_fret_copy
	j=n_before_fret_copy;
	boundary[j]=copy_cell_for_split_regions(boundary[i],1);
	interior[j]=find_interior_points(boundary[j]);
	fret_copy_type[i]=0; //This cell is the original
	fret_copy_type[j]=1; //This cell is a copy
	//Number of interior points might be slightly different than
	//n_points[i].
	n_p=0;
	for(p=interior[j];p!=NULL;p=(p->next))n_p++;
	n_points[j]=n_p;
	n_before_fret_copy++;
	n_found++;
	//Note that we moved all the upper cells up one index, so
	//they should still all be correctly offset from the lower
	//cells since we just incremented n_before_fret_copy
      }
    }

  }
  printf("-----> After comparison: Total lower=%i and total upper=%i\n",
	 n_before_fret_copy,n_found-n_before_fret_copy);

  //Done with the top and bottom fret images.

  //Make interior lists for +- radii, etc. and calculate some
  //statistics.

  //Push the next part off to a subroutine so can change arrays more
  //easily.
  calculate_global_stats_from_interior_and_boundary();

  *boundary_out=boundary;
  *interior_out=interior;

  return n_found;
}


/*********************************************************/
void calculate_global_stats_from_interior_and_boundary(){

  int i,j;
  float r_save[n_points_r_vs_theta];

  struct point *p;
  float tmp,tmp1,tmp_x,tmp_y;
  int n_p;



  for(i=0;i<n_found;i++){
    fft_stat[i]=((float)(FFT_ratio(boundary[i])));
    r_vs_theta_from_boundary(boundary[i]);

    for(j=0;j<n_points_r_vs_theta;j++){
      r_save[j]=r_vs_theta[j];
    }
    //Pause here to get some statistics using r_vs_theta
    circumference[i]=0.0; //Flag so calculate everything
    statistics_from_r_vs_theta(1,
			       (circumference+i), //Address in array
			       (major_axis1+i),
			       (minor_axis1+i),
			       (major_axis2+i),
			       (minor_axis2+i),
			       (vol_rotation+i));


    //Plus 1
    for(j=0;j<n_points_r_vs_theta;j++){
      r_vs_theta[j]=(r_save[j]+1.0);
    }
    //Before making new boundary_p1,m1,m2,m3, etc, free the previous
    //points. Note these are static variables so they're still around
    //from the previous calls. Do NOT free boundary[] or interior[]
    //here since those are saved in the struct blob pointer lists for
    //each cell and if you free them here those pointers will end up
    //pointing to garbage.
    point_list_free(boundary_p1[i]);
    point_list_free(interior_p1[i]);
    point_list_free(boundary_m1[i]);
    point_list_free(interior_m1[i]);
    point_list_free(boundary_m2[i]);
    point_list_free(interior_m2[i]);
    point_list_free(boundary_m3[i]);
    point_list_free(interior_m3[i]);
    point_list_free(boundary_p5[i]);
    point_list_free(boundary_phalfminor[i]);

    boundary_p1[i]=boundary_from_r_vs_theta();
    interior_p1[i]=find_interior_points(boundary_p1[i]);

    //Minus 1
    for(j=0;j<n_points_r_vs_theta;j++){
      r_vs_theta[j]=(r_save[j]-1.0);
      if (r_vs_theta[j]<1.0)r_vs_theta[j]=1.0;
    }
    boundary_m1[i]=boundary_from_r_vs_theta();
    interior_m1[i]=find_interior_points(boundary_m1[i]);

    //Minus 2
    for(j=0;j<n_points_r_vs_theta;j++){
      r_vs_theta[j]=(r_save[j]-2.0);
      if (r_vs_theta[j]<1.0)r_vs_theta[j]=1.0;
    }
    boundary_m2[i]=boundary_from_r_vs_theta();
    interior_m2[i]=find_interior_points(boundary_m2[i]);

    //Minus 3
    for(j=0;j<n_points_r_vs_theta;j++){
      r_vs_theta[j]=(r_save[j]-3.0);
      if (r_vs_theta[j]<1.0)r_vs_theta[j]=1.0;
    }
    boundary_m3[i]=boundary_from_r_vs_theta();
    interior_m3[i]=find_interior_points(boundary_m3[i]);

    //Calculate volume and some other parameters
    calculate_volume(interior[i],boundary[i],
		     (vol_sphere+i),
		     (vol_cone+i),
		     (vol_eff_1+i),
		     (vol_eff_2+i),
		     (vol_eff_3+i),
		     (vol_eff_4+i),
		     (vol_eff_5+i),
		     (vol_eff_6+i)
		     );
    //Re-calculate volume at r+1 to get surface area
    calculate_volume(interior_p1[i],boundary_p1[i],
		     (surface_area+i),
		     &tmp,
		     &tmp,
		     &tmp,
		     &tmp,
		     &tmp,
		     &tmp,
		     &tmp);
    surface_area[i]-=vol_sphere[i];

    //And we do +5 and +minor/2 for local background stuff. Here
    //we only calculate boundary points though.
    //Plus 5
    for(j=0;j<n_points_r_vs_theta;j++){
      r_vs_theta[j]=(r_save[j]+5.0);
    }
    boundary_p5[i]=boundary_from_r_vs_theta();

    //Plus minor/2
    tmp1=(minor_axis1[i])/2.0;
    for(j=0;j<n_points_r_vs_theta;j++){
      r_vs_theta[j]=(r_save[j]+tmp1);
    }
    boundary_phalfminor[i]=boundary_from_r_vs_theta();
  }


  //Calculate mean_x and mean_y of interior points
  for(i=0;i<n_found;i++){
    n_p=0;
    tmp_x=0.0;
    tmp_y=0.0;
    for(p=interior[i];p!=NULL;p=p->next){
      tmp_x+=( (float)(p->i) );
      tmp_y+=( (float)(p->j) );
      n_p++;
    }
    tmp=(float)n_p;
    if(n_p>0){
      mean_x[i]=tmp_x/tmp;
      mean_y[i]=tmp_y/tmp;
    }else{
      mean_x[i]=-1.0;
      mean_y[i]=-1.0;
    }
  }

  //output_data_to_tif_file("test.tif",atest,xmax,ymax,NULL,0,16,0,"");

  printf("\nTotal this time=%i.\n",n_found);

  //Test-test-test-test-asg-test-asg-test
  //for(i=0;i<xmax;i++){
  //  for(j=0;j<ymax;j++){
  //    n_p=(j*xmax+i);
  //    c[n_p]=(float)work_array[n_p];
  //  }
  //}
  //output_data_to_tif_file("test.tif",c,xmax,ymax,NULL,0,8,0,"");
  //exit(0);

  return;
}

/*************************************************************/
int recombination_check(int i_t0,
			int n_cuts,
			int *cut_type,
			int *cut_flag,
			float *cuts){
  //This routine runs over all the cells found in cs[] in the latest
  //position. It uses the cuts to look for good cells that have no
  //nucleus and tries to recombine them with good cells that have a
  //nucleus. It recombines what it finds and deletes the original cell.
  //It does this backwards in b->i_time all the way back to (and including)
  //i_t0.

  //The cut_types are hard-coded:
  //cut_type[]=0 is cut on f_nuc/a_nuc-back on whether it's a nucleus or not.
  //cut_type[]=1 is cut on flraw/area-back on whether it's a good cell or not.

  struct blob *b0;
  struct blob *b1;
  struct blob *btmp;
  struct point *ptmp;
  struct blob *bloop;

  int try_again=0;
  int ix,iy;
  int u,u0;

  int icut;
  int j,k;
  int has_nucleus;
  int good_cell;
  int recombine;
  int icombine;
  int flag;
  float tmp,tmp_n;
  float max_overlap;
  int isection;
  int k_remove;
  float overlap_cut=0.0;
  int overlap_value_save;

  int *cell_marker;
  int recombination_min_pixels=75;

  if (n_known<=0) return 0; //It would drop through anyway, but
  //what will malloc(0) do?!?!?!?

  cell_marker=(int *)malloc(n_known*sizeof(int));

  //Loop over every cell in the cs[] and decide where it's a good cell
  //and has a nucleus. Examine all the cuts, etc.
  //cell_marker[]=0-->bad cell
  //cell_marker[]=1-->good cell, no nucleus
  //cell_marker[]=2-->good cell, has nucleus
  for(j=0;j<n_known;j++){
    cell_marker[j]=2; //default to good cell with nucleus.
    b0=cs[j];
    //We're going to consider all previous times for this cell back to
    //and including i_t0.
    //We're looking for a good cell that has no nucleus. Keep going backwards
    //(the flags can change as you go backwards) until we have determined
    //that we have a bad cell (defaults to good) or that we have a nucleus.
    btmp=b0;
    good_cell=1;//Default to good cell with nucleus
    has_nucleus=1;

    do{
      if ((btmp->i_time)<i_t0) break; //Keep going back until get to i_t0

      //Look for all the cuts that are set for this flag type
      flag=btmp->flag;
      for(icut=0;icut<n_cuts;icut++){
	if (cut_flag[icut]==flag){ //We have a cut for this guy
	  if (cut_type[icut]==1){ //Hardcoded cut_types
	    //Is this a good cell (using flraw/area-back)
	    if ((btmp->a)<=0.0){
	      good_cell=0;
	      break;
	    }else{
	      tmp=(btmp->fluor)/(btmp->a)-back_pixels[btmp->i_time];
	      if (tmp<cuts[icut]){ //bad cell
		good_cell=0;
		break; //Not interested further in this cell so no more
		//cut checks necessary
	      }
	    }
	  }else if (cut_type[icut]==0){
	    //Does this cell have a nucleus
	    //(Note there can potentially be multiple cuts here. If it
	    //fails any of them we declare no nucleus)
	    tmp=(btmp->fl_nucleus3)/(btmp->area_nucleus3)-
	      back_pixels[btmp->i_time];
	    if (tmp<cuts[icut]){ //no nucleus
	      has_nucleus=0;
	    }
	  }else if (cut_type[icut]==2){
	    //Does this cell have a nucleus-->In this case consider
	    //all cells produced after time 0 to have NO nucleus and
	    //all produced at i_t=0 to have a nucleus.
	    for(bloop=btmp;bloop->prev!=NULL;bloop=bloop->prev);
	    if ((bloop->i_time)==0){
	      has_nucleus=1;
	    }else{
	      has_nucleus=0;
	    }
	  }

	} //End if-check if cut applies to current flag
      } //End loop over all the cuts

      if (good_cell==0){
	break; //Done with this cell
      }
      btmp=btmp->prev;
    }while(btmp!=NULL);

    //If cell is very small then always consider it a bud to be
    //recombined
    if ((b0->n)<recombination_min_pixels){
      has_nucleus=0; //Force it to look like no nucleus
    }
    //We don't want to recombine a cell that's
    //already been separated. It's already been separated if the cells exists
    //for times earlier than i_t0.
    //We got out of the above loop when
    //btmp->i_time==i_t0 (the while(btmp!=NULL) is really just a check and
    //it should never get there), so just check that there is another earlier
    //time point.
    if (btmp!=NULL){
      if(btmp->prev!=NULL){ //Is a previous time, so don't combine this guy
	has_nucleus=1; //Force to look like it has a nucleus
      }
    }

    //We've gone all the way back, check results of cuts
    if (good_cell==0){
      cell_marker[j]=0; //bad cell
    }else{
      if (has_nucleus==0){
	cell_marker[j]=1;//good cell, no nucleus
      }else{
	cell_marker[j]=2;//good cell, has nucleus
      }
    }
  } //end loop over all known cells

  printf("RECOMBINATION INFO: {ID,(0=bad,1=good, no nuc, 2=good, nuc)}\n");
  for(j=0;j<n_known;j++){
    printf("{%i,%i} ",cs[j]->index,cell_marker[j]);
  }
  printf("\n");

  //Now compare all {good cell, no nucleus} with {good cell, has nucleus}
  //to try to match who should be combined up.
  for(j=0;j<n_known;j++){
    if (cell_marker[j]!=1) continue; //1==>good_cell, no nucleus
    b0=cs[j];
    //See if b0 should be attached to b1 (assigned below)
    //Try seeing if boundary+1 overlaps any point of boundary+1 of b1.
    //This is probably faster to do by hand rather than doing the
    //r_vs_theta stuff.
    //-->Actually do it with the interior points+1. Then below, we
    //compare interior-to-interior if a cell has a previous time point,
    //and will end up (effectively)
    //doing boundary-to-boundary is cell doesn't have a previous time point.
    update_overlap_value();
    fill_overlap_array_with_point_list(b0->interior);
    //Now also include the boundary+1
    for(ptmp=(b0->boundary);ptmp!=NULL;ptmp=(ptmp->next)){
      ix=ptmp->i;
      iy=ptmp->j;
      if ((ix>=0)&&(ix<xmax)&&(iy>=0)&&(iy<ymax)){
	u0=iy*xmax+ix;
	u=u0;
	over[u]=overlap_value;
	if ((ix+1)<xmax){
	  u=u0+1;
	  over[u]=overlap_value;
	}
	if ((ix-1)>=0){
	  u=u0-1;
	  over[u]=overlap_value;
	}
	if ((iy+1)<ymax){
	  u=u0+xmax;
	  over[u]=overlap_value;
	}
	if ((iy-1)>=0){
	  u=u0-xmax;
	  over[u]=overlap_value;
	}
      }
    }

  try_again=0;
 try_again:
    recombine=0;
    tmp_n=((float)(b0->n));
    max_overlap=0.0;
    k_remove=0;
    for(k=0;k<n_known;k++){
      if (cell_marker[k]!=2) continue; //2=good_cell, has nucleus
      b1=cs[k];
      //Start with cells which are at the current time point
      if ((b0->i_time)!=(b1->i_time)) continue;

      //Go back to previous time if we can for this cell (if not
      //continue on) ("previous time" is actually i_t0-1 not the time
      //int the b1-> list because the time there might just be the
      //previous fluorescence image and we want the previous brightfield.)
      while(b1->prev!=NULL){
	b1=b1->prev;
	if (b1->i_time==(i_t0-1))break;
      }
      if (b1->i_time!=(i_t0-1)){ //We broke out of loop above because
	//we reached the end. Reset it to current time
	b1=cs[k];
      }
      isection=overlap(b1->interior,0,0);
      tmp=( (float) isection )/tmp_n;
      if (tmp>max_overlap){
	max_overlap=tmp;
	k_remove=k;
      }
    }

    if (max_overlap==0.0){ //Increase footprint of cell we're trying to
      //match and try again

      if (try_again>=10){
	printf("Tried 10x to match %i, giving up.\n",b0->index);
      }else{
	try_again++;
	//Add an extra ring to the boundary
	overlap_value_save=overlap_value;
	update_overlap_value();
	for(ix=0;ix<xmax;ix++){
	  for(iy=0;iy<ymax;iy++){
	    u0=iy*xmax+ix;
	    if (over[u0]==overlap_value_save){
	      //Add 4 sides
	      over[u0]=overlap_value;
	      if ((ix+1)<xmax){
		u=u0+1;
		over[u]=overlap_value;
	      }
	      if ((ix-1)>=0){
		u=u0-1;
		over[u]=overlap_value;
	      }
	      if ((iy+1)<ymax){
		u=u0+xmax;
		over[u]=overlap_value;
	      }
	      if ((iy-1)>=0){
		u=u0-xmax;
		over[u]=overlap_value;
	      }
	    }
	  }
	}

	goto try_again;
      }
    }

    if (max_overlap>overlap_cut){
      printf("Looking to combine %i with %i (%e)\n",
	     b0->index,cs[k_remove]->index,max_overlap);
      icombine=combine_cells_in_cs_array(j,k_remove,i_t0);
      if (icombine==0){ //We totally deleted cell j
	cs[j]=cs[n_known-1];
	cell_marker[j]=cell_marker[n_known-1];
	j--;
	n_known--;
      }
    }

  } //end loop over b0

  free(cell_marker);
  return 0;
}

/*************************************************************/
int combine_cells_in_cs_array(i0,i1,i_t0){
  //Take cell cs[i0] and attach it to cs[i1]. Do it all the way back
  //to i_t0. After attaching for each time, delete that cell. If the
  //cell is totally deleted return 0 otherwise return 1;

  //This will screw-up the n_found arrays found in find_cells. Ie, after
  //calling this, you can't rely on the information from find_cell() and
  //have to run it again.

  struct blob *b0;
  struct blob *b1;


  struct point *p0;
  struct point *p1;
  struct point *p0tmp;
  struct point *p1tmp;
  struct point *p0max;
  struct point *p1max;
  struct point *p0min;
  struct point *p1min;
  struct point *p0_a;
  struct point *p0_b;
  struct point *p1_a;
  struct point *p1_b;
  struct point *ptmp;

  int i;

  unsigned int ix1,iy1,dx,dy,d,dmin,dmax;

  b0=cs[i0];
  b1=cs[i1];

  if ((b0==NULL)||(b1==NULL)){
    printf("**************************************\n");
    printf("*** Tried attaching NULL cells together!\n");
    printf("*** This shouldn't ever happen!\n");
    printf("**************************************\n");
    return 1;
  }
  printf("Attaching cell %i to cell %i and flagging cell %i.\n",
	 b0->index,b1->index,b0->index);

  if ((b0->i_time)!=(b1->i_time)){
    printf("Tried to combine two cells from different i_t!\n");
    return 1; //Do nothing
  }

  //First try to attach boundaries. If this fails then do nothing.
  //Remove all points within distance 2 of the other side.
  p0=b0->boundary; //May change
  p1=b1->boundary;
  //Start by circularizing both boundaries
  for (p0tmp=p0;p0tmp->next!=NULL;p0tmp=p0tmp->next);
  p0tmp->next=p0;
  p0->prev=p0tmp;
  for (p1tmp=p1;p1tmp->next!=NULL;p1tmp=p1tmp->next);
  p1tmp->next=p1;
  p1->prev=p1tmp;
  //Now find two points farthest apart to start
  dmax=0;
  for(p1tmp=p1;p1tmp->next!=p1;p1tmp=p1tmp->next){
    ix1=p1tmp->i;
    iy1=p1tmp->j;
    for(p0tmp=p0;p0tmp->next!=p0;p0tmp=p0tmp->next){
      dx=(p0tmp->i)-ix1;
      dy=(p0tmp->j)-iy1;
      d=(dx*dx+dy*dy);
      if (d>dmax){
	dmax=d;
	p0max=p0tmp;
	p1max=p1tmp;
      }
    }
  }
  //Now start at the maxima and find first points in either direction
  //whose distance is < 2 away from p0.
  p0_a=NULL;
  p1_a=NULL;
  for(p1tmp=p1max;p1tmp->next!=p1max;p1tmp=p1tmp->next){
    ix1=p1tmp->i;
    iy1=p1tmp->j;
    dmin=xmax*ymax+1; //default to something high
    for(p0tmp=p0;p0tmp->next!=p0;p0tmp=p0tmp->next){
      dx=(p0tmp->i)-ix1;
      dy=(p0tmp->j)-iy1;
      d=(dx*dx+dy*dy);
      if (d<dmin){
	dmin=d;
	p0min=p0tmp;
	p1min=p1tmp;
      }
    }
    //Find first point with distance<2 going in this direction.
    if (dmin<=5){
      p0_a=p0min; //First combine point
      p1_a=p1min;
      break; //Go other direction now
    }
  }
  if ((p0_a==NULL)||(p1_a==NULL)){
    printf("Cells too far apart to combine\n");
    //de-circularize before returning
    (p0->prev)->next=NULL;
    p0->prev=NULL;
    (p1->prev)->next=NULL;
    p1->prev=NULL;
    return 1;
  }

  //Now do the same going the other direction around p1
  p0_b=NULL;
  p1_b=NULL;
  dmin=xmax*ymax+1;
  for(p1tmp=p1max;p1tmp->prev!=p1max;p1tmp=p1tmp->prev){
    ix1=p1tmp->i;
    iy1=p1tmp->j;
    dmin=xmax*ymax+1;
    for(p0tmp=p0;p0tmp->next!=p0;p0tmp=p0tmp->next){
      dx=(p0tmp->i)-ix1;
      dy=(p0tmp->j)-iy1;
      d=(dx*dx+dy*dy);
      if (d<dmin){
	dmin=d;
	p0min=p0tmp;
	p1min=p1tmp;
      }
    }
    //Find first point with distance<2 going in this direction around p1
    if (dmin<=5){
      p0_b=p0min; //First combine point
      p1_b=p1min;
      break;
    }
  }
  if ((p0_b==NULL)||(p1_b==NULL)){
    printf("Cells too far apart to combine\n");
    //de-circularize before returning
    (p0->prev)->next=NULL;
    p0->prev=NULL;
    (p1->prev)->next=NULL;
    p1->prev=NULL;
    return 1;
  }

  //Check if the same point was found
  if (p0_a==p0_b)p0_b=p0_a->prev; //Put p0_b is upstream of p0_a,
  //so far p0_a,p0_b directionality has been irrelevant. Below however
  //we're assuming than when we go from a->b on p0 we cross p0max, which
  //implies that p0_b is after p0_a. The opposite is true for the p1 loop
  //where above we've made it so that p1_b should be downstream of p1_a.
  if (p1_a==p1_b)p1_b=p1_a->next; //but p1_b is downstream of p1_a

  //Delete all points between p1_a and p1_b (exclusive)
  p0tmp=p1_a->next;
  while(p0tmp!=p1_b) {
    p1tmp=p0tmp;
    p0tmp=p0tmp->next;
    point_free(p1tmp);
  }
  //We're going to connect p1_a to p0_a and p1_b to p0_b
  //Make sure the p0 loop is in the right direction. The right
  //direction is when we go from p0_a to p0_b we cross p0max.
  ix1=0;
  for(p0tmp=p0_a;p0tmp!=p0_b;p0tmp=p0tmp->next){
    if (p0tmp==p0max){
      ix1=1;
      break;
    }
  }
  if (ix1==0){ //Need to reverse direction of p0 loop.
    p0tmp=p0_a;
    do{
      ptmp=p0tmp->next;
      p0tmp->next=p0tmp->prev;
      p0tmp->prev=ptmp;
      p0tmp=ptmp;
    }while(p0tmp!=p0_a);
  }
  //Delete all points between p0_b and p0_a (exclusive)
  p0tmp=p0_b->next;
  while(p0tmp!=p0_a){
    p1tmp=p0tmp;
    p0tmp=p0tmp->next;
    point_free(p1tmp);
  }
  //connect up the lists
  p1_a->next=p0_a;
  p0_a->prev=p1_a;
  p0_b->next=p1_b;
  p1_b->prev=p0_b;

  //De-dircularize at arbitrary place
  (p1_a->prev)->next=NULL;
  (p1_a->prev)=NULL;
  p1=p1_a;
  p1tmp=find_interior_points(p1);

  //Change the n_found list partially (so is screwed up after this also)
  //because we might be outputting the cells later.
  ix1=-1;
  iy1=-1;
  for(i=0;i<n_found;i++){
    if (location_in_cs_array[i]==i0)ix1=i;
    if (location_in_cs_array[i]==i1)iy1=i;
  }
  if ((iy1<0)||(ix1<0)){
    printf("************** Didn't find cells in n_found arrays.\n");
  }

  //point_list_free(interior[ix1]);
  //point_list_free(boundary[ix1]);
  //point_list_free(interior[iy1]);
  //point_list_free(boundary[iy1]);
  boundary[iy1]=p1;
  interior[iy1]=p1tmp;
  boundary[ix1]=NULL;
  interior[iy1]=NULL;

  //Now loop over all the previous blobs in linked list
  do{
    b1->boundary=p1;
    b1->interior=p1tmp;
    b1->fluor+=b0->fluor;
    b1->a+=b0->a; //Won't be exactly the area because of combination
    //but is exactly the pixels used to calculate fluorescence (which I
    //added.
    b1->n+=b0->n;
    //A bunch of guys that approximately add
    b1->fluor_p1+=b0->fluor_p1;
    b1->fluor_m1+=b0->fluor_m1;
    b1->fluor_m2+=b0->fluor_m2;
    b1->fluor_m3+=b0->fluor_m3;
    b1->area_p1+=b0->area_p1;
    b1->area_m1+=b0->area_m1;
    b1->area_m2+=b0->area_m2;
    b1->area_m3+=b0->area_m3;
    b1->vacuole_fl+=b0->vacuole_fl;
    b1->vacuole_area+=b1->vacuole_area;
    b1->vol_rotation+=b0->vol_rotation;
    b1->vol_cone+=b0->vol_cone;
    b1->vol_sphere+=b0->vol_sphere;
    b1->surface_area+=b0->surface_area;
    b1->vol_eff_1+=b0->vol_eff_1;
    b1->vol_eff_2+=b0->vol_eff_2;
    b1->vol_eff_3+=b0->vol_eff_3;
    b1->vol_eff_4+=b0->vol_eff_4;
    b1->vol_eff_5+=b0->vol_eff_5;
    b1->vol_eff_6+=b0->vol_eff_6;
    b1->circumference+=b0->circumference;

    //Now delete b0. Effectively delete by making its (x0,y0) the negative
    //if the id we attached it to. We'll check for this as a flag below
    //in update_list_of_found_cells().
    b0->x=-1000.0-((float)(b1->index));
    b0->y=-1000.0-((float)(b1->index));
    b0=b0->prev;
    b1=b1->prev;
    if (b0==NULL) return 1;
    if (b1==NULL) return 1;
    if ((b1->i_time)<i_t0) break;

  } while(1);

  return 1;
}




/*************************************************************/
void internal_structure(int flag,int flag2){
  //Search through the found cells and within
  //each cell, make some statistics related to the internal structure
  //of the fluorescence in that cell.

  //First we look for a nucleus based on brightest circle of fixed
  //radius in fluorescence.
  //Then we calculate the membrane fluorescence and how many pixels
  //went into that sum.

  //flag2==0 means use geometric center, 1 means use image according to flag
  //flag==0 means use the FL image to find nucleus, or, for
  //        the fret case, use the lower part of fret image.
  //flag==1 means use the upper part of the FRET image


  int i;
  int k;
  int u;
  int itype;
  int ix,iy;

  float *array;
  int npoints;
  float sum1,sum2;

  int nloop_low,nloop_high;
  int *list_x;
  int *list_y;
  int n_list;

  float nucleus_radius;

  //int do_gauss=(nucleus_distribution_types-1);
  //  int do_contiguous=(nucleus_distribution_types-1);
  //int do_contiguous=-999; //Drop this method 7/29/03
  //int all_pixels=(nucleus_distribution_types-2);

  int vec_x,vec_y;
  int icopy;

  //float gauss_x,gauss_y,gauss_sig;
  //float back;

  float tmpx,tmpy;
  int n_p;

  //double dsum1,dsum2;
  //double f1,f2,wt,sum_wt;
  //double dtmp1,dtmp2,dr2,half_inv_sig2;
  //double mu_x,mu_y;
  struct point *p;

  static int first=1;

  if (first==1){
    first=0;
    printf("Initializing nuclear centers array.\n");
    for(i=0;i<max_cells;i++){
      p_nuclear_center[i]=NULL;
    }
  }

  if (flag>1){
    printf("Unknown flag value in internal_structure(%i)\n",flag);
    return;
  }

  //Loop over each of the found cells
  if (have_fret_image==1){
    if (flag==0){ //Using lower part of image
      nloop_low=0;
      nloop_high=n_before_fret_copy;
    }else if (flag==1){
      nloop_low=n_before_fret_copy;
      nloop_high=n_found;
    }
  }else{
    nloop_low=0;
    nloop_high=n_found;
  }

  if (flag2==0){//Use geometric center of the cell
    for (i=nloop_low;i<nloop_high;i++){
      //printf("----> %i %e %e\n",i,mean_x[i],mean_y[i]);fflush(stdout);
      if(p_nuclear_center[i]==NULL){
				p_nuclear_center[i]=(struct point *)malloc(sizeof(struct point));
      }
      //Recalculate geometric center in case re-aligning cells
      tmpx=0.0;
      tmpy=0.0;
      n_p=0;
      for(p=interior[i];p!=NULL;p=p->next){
				tmpx+=( (float)(p->i) );
				tmpy+=( (float)(p->j) );
				n_p++;
      }
      tmpx/=((float)n_p);
      tmpy/=((float)n_p);
      p_nuclear_center[i]->i=(int)(tmpx+0.5);
      p_nuclear_center[i]->j=(int)(tmpy+0.5);
    }
  }else{ // (flag2==1) use image acording to flag
    if (recalculate_internal==1){ //We haven't searched these
      //cells yet. Set centers to be NULL so they get set below
      for (i=nloop_low;i<nloop_high;i++){
				free(p_nuclear_center[i]);
				p_nuclear_center[i]=NULL;
      }
    }
  }
  recalculate_internal=0; //Don't do again unless reset in cell.c

  //If there is a nuclear_labelled image use that, otherwise
  //use the FL[] image itself.
  if ((third_image!=NULL)&&
      ((third_image_type==nuclear_label)||
       (third_image_type==vacuole_label))){
    array=third_image;
    printf("Internal structure using third image.\n");
  }else{
    array=fl;
    printf("Internal structure using fluorescence image.\n");
  }

  for(i=nloop_low;i<nloop_high;i++){ //Loop over all the cells

    for(itype=0;itype<nucleus_distribution_types;itype++){

			//V1.4.5 ussing the no Gaussian code only (see commented code below)
				//Now loop over each interior point and maximize a disk of
				//of fluorescence centered at that point.  Make sure the disk
				//is entirely within the cell.
	  		nucleus_radius=((float)nucleus_radii[itype]);

				maximum_pixels_within_fixed_radius(array,xmax,ymax,
						   interior[i],
					  	 nucleus_radius,
					   	 &list_x,&list_y,&n_list,
					   	 (p_nuclear_center+i));
				//list_x[],list_y[] are a set of n_list points that start at the
				//center of a circle of radius nucleus_radius that have the brightest
				//per_pixel values of array[].
				//*(p_nuclear_center) contains new value if started out NULL.

				if(itype==0){
					x_nucl[i]=p_nuclear_center[i]->i;
					y_nucl[i]=p_nuclear_center[i]->j;
				}

				area_nucleus[i][itype]=((float)n_list);
				sum1=0.0;
				sum2=0.0;
				for (k=0;k<n_list;k++){
	  			u=(list_y[k]*xmax)+(list_x[k]);
	  			sum1+=fl[u];
	  			sum2+=array[u];
				}
				fl_nucleus[i][itype]=sum1; //Value from fluorescence array
				fl_nucleus_from_search[i][itype]=sum2; //The value from search image

				if ((array==third_image)&&(third_image_type==nuclear_label)&&
	    		(itype==0)){
	  			//Mark nuclear boundaries in d[] array
	  			//printf("For cell %i: %i pixels\n",i,n_list);fflush(stdout);
	  			for (k=0;k<n_list;k++){
	    			u=(list_y[k]*xmax)+(list_x[k]);
	    			d[u]=cell_nucleus;
	    			//if (i==0){
	    			//  printf("Nucleus: (x,y)=(%i,%i)\n",list_x[k],list_y[k]);
	    			//}
	  			}
	  			//	}else if ((array==third_image)&&(itype==do_contiguous)&&
	 			  //		  (third_image_type==vacuole_label)){
	  			//	  //For the vacuole label, use the do_contiguous[] method
	  			//for (k=0;k<n_list;k++){
	  			//  u=(list_y[k]*xmax)+(list_x[k]);
	  			//  d[u]=cell_nucleus;
	  			//  //if (i==0){
	  			//  //  printf("Nucleus: (x,y)=(%i,%i)\n",list_x[k],list_y[k]);
	  			//  //}
	  			//	}
				}
				if (have_fret_image==1){
	  			//Use the region found above for the other part of the
	  			//split images (after offsetting to the other part).
	  			if (flag==0){
	    			icopy=i+n_before_fret_copy; //bottom to top
	  			}else{
	    			icopy=i-n_before_fret_copy; //top to bottom
	  			}
	  			vec_x=((int)(mean_x[icopy]-mean_x[i])); //To go from cell i
	  			vec_y=((int)(mean_y[icopy]-mean_y[i])); //to cell icopy
	  			//if ((fret_copy_type[i]==0)&&
	  			//    (fret_copy_type[icopy]==0)){
	  			//  printf("%e   %e   %e\n",mean_x[i],
	  			//		   (mean_x[icopy]-mean_x[i]),
	  			//	   (mean_y[icopy]-mean_y[i]));fflush(stdout);
	  			//}

	  			npoints=0;
	  			sum1=0.0;
	  			sum2=0.0;
	  			for (k=0;k<n_list;k++){
	    			ix=(list_x[k])+vec_x;
	    			iy=(list_y[k])+vec_y;
	    			//Make sure still in bounds (just in case)
	    			if ((ix>=0)&&(ix<xmax)&&(iy>=0)&&(iy<ymax)){
	      			npoints++;
	      			u=(iy*xmax)+ix;
	      			sum1+=fl[u];
	     				sum2+=array[u];
	      			if (itype==0)(d[u]=cell_nucleus);
	    			}
	  			}
	  			area_nucleus[icopy][itype]=((float)npoints);
	  			fl_nucleus[icopy][itype]=sum1;
	  			fl_nucleus_from_search[icopy][itype]=sum2;
				}





			/* V1.4.5 Commented out
			//find best center, next times use that center.
      if ((itype!=do_contiguous)&&(itype!=do_gauss)){
        //do_gauss=nucleus_distribution_type - 1
				//Now loop over each interior point and maximize a disk of
				//of fluorescence centered at that point.  Make sure the disk
				//is entirely within the cell.
				if (itype<all_pixels){ //all_pixels=nucleus_distribution_types - 2
	  			nucleus_radius=(2.0+((float)itype));
				}else{
	  			nucleus_radius=-1.0; //Will simply return all pixels in interior[]
				}
				maximum_pixels_within_fixed_radius(array,xmax,ymax,
						   interior[i],
					  	 nucleus_radius,
					   	 &list_x,&list_y,&n_list,
					   	 (p_nuclear_center+i));
				//list_x[],list_y[] are a set of n_list points that start at the
				//center of a circle of radius nucleus_radius that have the brightest
				//per_pixel values of array[].
				// *(p_nuclear_center) contains new value if started out NULL.

			}else if (itype==do_contiguous){ //Do contiguous cut method
				maximum_contiguous_pixels(array,xmax,ymax,interior[i],
						  &list_x,&list_y,&n_list);
      }	else if (itype==do_gauss){
				if (third_image_type==vacuole_label){
	  			maximum_contiguous_pixels(array,xmax,ymax,interior[i],
				    	&list_x,&list_y,&n_list);

				}else{
	  			//Assume that we've already done the maximum_pixels_.... stuff
	  			back=(area_nucleus[i][all_pixels]-area_nucleus[i][all_pixels-2]);
	  			if (back<=0.0){
	    			back=0.0;
	  			}else{
	    			back=(fl_nucleus_from_search[i][all_pixels]-
	    					fl_nucleus_from_search[i][all_pixels-2])/back;
	  			}
	  			//Back is now estimated fluorescence/pixels of background
	  			get_gauss_2d_parameters(array,xmax,ymax,back,interior[i],
						  &gauss_x,&gauss_y,&gauss_sig,
					  	&list_x,&list_y,&n_list);

	  			if (gauss_sig<=0.0)gauss_sig=1.0e5; //Essentially flat

	  			//printf("------->back=%e and (%i,%i)-->(%e,%e,%e)\n",
	  			//       back,
	  			//       p_nuclear_center[i]->i,
	  			//       p_nuclear_center[i]->j,
	  			//       gauss_x,gauss_y,gauss_sig);fflush(stdout);
				}

				if ((array==third_image)&&(third_image_type==vacuole_label)){
	  			//Mark gaussian sigma on images
	  			for (k=0;k<n_list;k++){
	    			u=(list_y[k]*xmax)+(list_x[k]);
	    			d[u]=cell_nucleus;
	    			//if (i==0){
	    			//  printf("Nucleus: (x,y)=(%i,%i)\n",list_x[k],list_y[k]);
	    			//}
	  			}
				}

      }


      //Sum-up the data for non-gaussian types (or for
      //vacuole_label types do for do_gauss also)
      if ((itype!=do_gauss)||((itype==do_gauss)&&
			     (third_image_type==vacuole_label))){
				area_nucleus[i][itype]=((float)n_list);
				sum1=0.0;
				sum2=0.0;
				for (k=0;k<n_list;k++){
	  			u=(list_y[k]*xmax)+(list_x[k]);
	  			sum1+=fl[u];
	  			sum2+=array[u];
				}
				fl_nucleus[i][itype]=sum1; //Value from fluorescence array
				fl_nucleus_from_search[i][itype]=sum2; //The value from search image

				if ((array==third_image)&&(third_image_type==nuclear_label)&&
	    		(itype==0)){
	  			//Mark nuclear boundaries in d[] array
	  			//printf("For cell %i: %i pixels\n",i,n_list);fflush(stdout);
	  			for (k=0;k<n_list;k++){
	    			u=(list_y[k]*xmax)+(list_x[k]);
	    			d[u]=cell_nucleus;
	    			//if (i==0){
	    			//  printf("Nucleus: (x,y)=(%i,%i)\n",list_x[k],list_y[k]);
	    			//}
	  			}
	  			//	}else if ((array==third_image)&&(itype==do_contiguous)&&
	 			  //		  (third_image_type==vacuole_label)){
	  			//	  //For the vacuole label, use the do_contiguous[] method
	  			//for (k=0;k<n_list;k++){
	  			//  u=(list_y[k]*xmax)+(list_x[k]);
	  			//  d[u]=cell_nucleus;
	  			//  //if (i==0){
	  			//  //  printf("Nucleus: (x,y)=(%i,%i)\n",list_x[k],list_y[k]);
	  			//  //}
	  			//	}
				}
				if (have_fret_image==1){
	  			//Use the region found above for the other part of the
	  			//split images (after offsetting to the other part).
	  			if (flag==0){
	    			icopy=i+n_before_fret_copy; //bottom to top
	  			}else{
	    			icopy=i-n_before_fret_copy; //top to bottom
	  			}
	  			vec_x=((int)(mean_x[icopy]-mean_x[i])); //To go from cell i
	  			vec_y=((int)(mean_y[icopy]-mean_y[i])); //to cell icopy
	  			//if ((fret_copy_type[i]==0)&&
	  			//    (fret_copy_type[icopy]==0)){
	  			//  printf("%e   %e   %e\n",mean_x[i],
	  			//		   (mean_x[icopy]-mean_x[i]),
	  			//	   (mean_y[icopy]-mean_y[i]));fflush(stdout);
	  			//}

	  			npoints=0;
	  			sum1=0.0;
	  			sum2=0.0;
	  			for (k=0;k<n_list;k++){
	    			ix=(list_x[k])+vec_x;
	    			iy=(list_y[k])+vec_y;
	    			//Make sure still in bounds (just in case)
	    			if ((ix>=0)&&(ix<xmax)&&(iy>=0)&&(iy<ymax)){
	      			npoints++;
	      			u=(iy*xmax)+ix;
	      			sum1+=fl[u];
	     				sum2+=array[u];
	      			if (itype==0)(d[u]=cell_nucleus);
	    			}
	  			}
	  			area_nucleus[icopy][itype]=((float)npoints);
	  			fl_nucleus[icopy][itype]=sum1;
	  			fl_nucleus_from_search[icopy][itype]=sum2;
				}

      }else{ //For summing the data with gaussian weights

				mu_x=((double)gauss_x);
				mu_y=((double)gauss_y);
				half_inv_sig2=0.5/((double)(gauss_sig*gauss_sig));
				sum_wt=0.0;
				dsum1=0.0;
				dsum2=0.0;
				for(p=interior[i];p!=NULL;p=(p->next)){
	  			ix=(p->i);
	  			iy=(p->j);
	  			if ((ix>=0)&&(ix<xmax)&&(iy>=0)&&(iy<ymax)){
	    			dtmp1=((double)ix)-mu_x;
	    			dtmp2=((double)iy)-mu_y;
	   	 			dr2=(dtmp1*dtmp1+dtmp2*dtmp2);
	    			wt=exp(-dr2*half_inv_sig2);
	    			sum_wt+=wt;

	    			u=(iy*xmax+ix);
	    			f1=((double)fl[u]);
	    			f2=((double)array[u]);

	    			dsum1+=(f1*wt);
	    			dsum2+=(f2*wt);

	  			}
				}
				area_nucleus[i][itype]=3.14159*gauss_sig*gauss_sig;
				//Average intensity from fluorescence array
				if ((npoints>0)&&(gauss_sig<1000.0)){
	  			fl_nucleus[i][itype]=((float)(dsum1/sum_wt));
	  			fl_nucleus_from_search[i][itype]=
	    			((float)(dsum2/sum_wt)); //The value from search image
				}else{
	  			fl_nucleus[i][itype]=0.0;
	  			fl_nucleus_from_search[i][itype]=0.0;
				}
				if (have_fret_image==1){
	  			//Use the region found above for the other part of the
	  			//split images (after offsetting to the other part).
	 				if (flag==0){
	    			icopy=i+n_before_fret_copy; //bottom to top
	  			}else{
	    			icopy=i-n_before_fret_copy; //top to bottom
	  			}
	  			mu_x+=((double)(mean_x[icopy]-mean_x[i]));
	  			mu_y+=((double)(mean_y[icopy]-mean_y[i]));

	  			//Repeat sum for icopy image
	  			sum_wt=0.0;
	  			dsum1=0.0;
	  			dsum2=0.0;
	  			npoints=0;
	  			for(p=interior[icopy];p!=NULL;p=(p->next)){
	    			ix=(p->i);
	    			iy=(p->j);
	    			if ((ix>=0)&&(ix<xmax)&&(iy>=0)&&(iy<ymax)){
	      			npoints++;
	      			dtmp1=((double)ix)-mu_x;
	     	 			dtmp2=((double)iy)-mu_y;
	      			dr2=(dtmp1*dtmp1+dtmp2*dtmp2);
	      			wt=exp(-dr2*half_inv_sig2);
	      			sum_wt+=wt;

	      			u=(iy*xmax+ix);
	      			f1=((double)fl[u]);
	      			f2=((double)array[u]);

	      			dsum1+=(f1*wt);
	      			dsum2+=(f2*wt);

	    			}
	  			}
	  			area_nucleus[icopy][itype]=3.14159*gauss_sig*gauss_sig;
	  			if ((npoints>0)&&(gauss_sig<1000.0)){
	    			fl_nucleus[icopy][itype]=
	      				((float)(dsum1/sum_wt)); //Value from fluorescence array
	    			fl_nucleus_from_search[icopy][itype]=
	      				((float)(dsum2/sum_wt)); //The value from search image
	  			}else{
	    			fl_nucleus[icopy][itype]=0.0;
	    			fl_nucleus_from_search[icopy][itype]=0.0;
	  			}


				}
      }//End if of whether doing normal sum or gaussian weighted sums
			*/ //End V1.4.5 commented out

    }//End loop over different types of searches


    // printf("----->Cell %i: %e-->%e\n",
    //	   i,
    //	   fl_nucleus[i][3]/area_nucleus[i][3],
    //	   fl_nucleus[i][do_gauss]/area_nucleus[i][do_gauss]);fflush(stdout);


  }//End loop over all the found cells
  return;
}



/*************************************************************/
void next_prev_fl_comparison(){
  //Collect information by comparing
  //pixel-by-pixel the found cells and the same
  //pixel locations in the previous image.
  int i,j;
  int k,l;
  int n;
  int u1;
  int u2;
  int offx,offy;
  float sum0,sum1;

  double mu_x,mu_y;
  double pos_mean_x;
  double pos_mean_y;
  double pos_mean_r;
  double pos_rms_r;
  double neg_mean_x;
  double neg_mean_y;
  double neg_mean_r;
  double neg_rms_r;
  double tot_pos,tot_neg;
  double curprev2;
  double dtmp1,dtmp2;

  struct point *p;

  if (flprev==NULL){
    flprev=(float *)malloc(xmax_ymax*sizeof(float));
    cur_prev=(float *)malloc(xmax_ymax*sizeof(float));
  }
  first_cur_prev_comparison++;
  if(first_cur_prev_comparison==1){
    for(i=0;i<n_found;i++){
      //Set all statistics we would have calculated here to -1.0
      pos_sig_mean_x[i]=0.0;
      pos_sig_mean_y[i]=0.0;
      pos_sig_mean_r[i]=0.0;
      pos_sig_rms_r[i]=0.0;
      neg_sig_mean_x[i]=0.0;
      neg_sig_mean_y[i]=0.0;
      neg_sig_mean_r[i]=0.0;
      neg_sig_rms_r[i]=0.0;
    }
    goto save_current_as_previous;
  }

  //Calculate an offset from this image to the previous
  i=0;
  j=0;

  do{
    offx=i;
    offy=j;
    sum0=diff_fl_prev(i,j);
    while((sum1=diff_fl_prev(++i,j))<sum0){ //++i increments before call
      sum0=sum1;
    }
    i--; //Went one too far.
    if (i==offx){ //Didn't move, try other direction
      while((sum1=diff_fl_prev(--i,j))<sum0){ //++i increments before call
	sum0=sum1;
      }
      i++; //Went one too far.
    }

    while((sum1=diff_fl_prev(i,++j))<sum0){
      sum0=sum1;
    }
    j--; //Went one too far.
    if (j==offy){ //Didn't move, try other direction
      while((sum1=diff_fl_prev(i,--j))<sum0){
	sum0=sum1;
      }
      j++; //Went one too far.
    }

  }while((i!=offx)||(j!=offy));

  printf("offset=(%i,%i)\n",offx,offy);

  //Now calculate a "significance" of the difference. Ie, calculate
  //delta/sigma pixel by pixel. Assume sigma=sqrt(N1+N2).
  for(i=0;i<xmax;i++){
    for(j=0;j<ymax;j++){
      k=i+offx;
      l=j+offy;
      if((k<0)||(k>=xmax))k=i;
      if((l<0)||(l>=ymax))l=j;
      if(j<=2){  //For j==1 (and ~ 2) the images seem to be set to zero
        k=i;
        l=j;
      }
      //4/22/03--Changed this so doesn't divide by sum, but is
      //just difference.
      u1=(j*xmax+i);
      u2=(l*xmax+k);
      cur_prev[u1]=(fl[u2]-flprev[u1]);
      //	/((float)sqrt((double)((fl[k][l]+flprev[i][j]))));
    }
  }

  /*
  //asg--test
  for(i=0;i<xmax;i++){
    for(j=0;j<ymax;j++){
      c[i][j]=cur_prev[i][j]+1000.0;
    }
  }
  if(first_cur_prev_comparison==2){
    output_data_to_tif_file("tmp1.tif",0,16,"");
  }else
  if(first_cur_prev_comparison==3){
    output_data_to_tif_file("tmp2.tif",0,16,"");
  }else
  if(first_cur_prev_comparison==4){
    output_data_to_tif_file("tmp3.tif",0,16,"");
  }else
  if(first_cur_prev_comparison==5){
    output_data_to_tif_file("tmp4.tif",0,16,"");
  }
  */


  //Use the new "significance" plot to calculate some statistics for
  //each cell we found.
  for(n=0;n<n_found;n++){ //Loop over all the cells

    mu_x=((double)mean_x[n]);
    mu_y=((double)mean_y[n]);

    pos_mean_x=0.0;
    pos_mean_y=0.0;
    pos_mean_r=0.0;
    pos_rms_r=0.0;
    neg_mean_x=0.0;
    neg_mean_y=0.0;
    neg_mean_r=0.0;
    neg_rms_r=0.0;
    tot_pos=0.0;
    tot_neg=0.0;
    for(p=interior[n];p!=NULL;p=p->next){
      i=(p->i);
      j=(p->j);
      if((i>=0)&&(j>=0)&&(i<xmax)&&(j<ymax)){
	u1=(j*xmax+i);
	curprev2=((double)cur_prev[u1]);
	//curprev2*=curprev2;
	dtmp1=(((double)i)-mu_x);
	dtmp2=(((double)j)-mu_y);

	//if (cur_prev[i][j]>((float)sqrt((double)(fl[i][j])))){
	if (cur_prev[u1]>0.0){
	  pos_mean_x+=(dtmp1*curprev2);
	  pos_mean_y+=(dtmp2*curprev2);
	  pos_mean_r+=(sqrt(dtmp1*dtmp1+dtmp2*dtmp2)*curprev2);
	  pos_rms_r+=((dtmp1*dtmp1+dtmp2*dtmp2)*curprev2);
	  tot_pos+=curprev2;
	}else{// if (cur_prev[i][j]<=0.0){
	  neg_mean_x+=(dtmp1*curprev2);
	  neg_mean_y+=(dtmp2*curprev2);
	  neg_mean_r+=(sqrt(dtmp1*dtmp1+dtmp2*dtmp2)*curprev2);
	  neg_rms_r+=((dtmp1*dtmp1+dtmp2*dtmp2)*curprev2);
	  tot_neg+=curprev2;
	}
      }
    }

    if (tot_pos>0.0){
      pos_mean_x/=tot_pos;
      pos_mean_y/=tot_pos;
      pos_mean_r/=tot_pos;
      pos_rms_r=(pos_rms_r/tot_pos-(pos_mean_r*pos_mean_r));
      if (pos_rms_r>=0.0){
	pos_rms_r=sqrt(pos_rms_r);
      }else{
	printf("Negative positive RMS! (%e)\n",pos_rms_r);
      }
    }else{
      pos_mean_x=-10.0;
      pos_mean_y=-10.0;
      pos_mean_r=-10.0;
      pos_rms_r=-10.0;
    }
    if (tot_neg>0.0){
      neg_mean_x/=tot_neg;
      neg_mean_y/=tot_neg;
      neg_mean_r/=tot_neg;
      neg_rms_r=(neg_rms_r/tot_neg-(neg_mean_r*neg_mean_r));
      if (neg_rms_r>=0.0){
	neg_rms_r=sqrt(neg_rms_r);
      }else{
	printf("Negative negative RMS! (%e)\n",neg_rms_r);
      }
    }else{
      neg_mean_x=-10.0;
      neg_mean_y=-10.0;
      neg_mean_r=-10.0;
      neg_rms_r=-10.0;
    }

    pos_sig_mean_x[n]=((float)pos_mean_x);
    pos_sig_mean_y[n]=((float)pos_mean_y);
    pos_sig_mean_r[n]=((float)pos_mean_r);
    pos_sig_rms_r[n]=((float)pos_rms_r);
    neg_sig_mean_x[n]=((float)neg_mean_x);
    neg_sig_mean_y[n]=((float)neg_mean_y);
    neg_sig_mean_r[n]=((float)neg_mean_r);
    neg_sig_rms_r[n]=((float)neg_rms_r);

  }

 save_current_as_previous:
  for(i=0;i<xmax_ymax;i++){
    flprev[i]=fl[i];
  }

  return;
}

/*************************************************************/
float diff_fl_prev(int offx, int offy){
  //Calculate a statistic to look to see if should move images
  //by overall offset. (Comparing two fluorescent images for this
  //purpose....)
  int i,j,k,l;
  int u1,u2;
  float sum;
  float tmp;
  float total;

  sum=0.0;
  total=0.0;
  for(i=2;i<xmax-2;i++){
    for(j=2;j<ymax-2;j++){
      k=i+offx;
      l=j+offy;
      if((k>=0)&&(k<xmax)&&(l>=0)&&(l<ymax)){
	u1=(j*xmax+i);
	u2=(l*xmax+k);
        tmp=(fl[u2]-flprev[u1]);
        sum+=(tmp*tmp);
        total+=1.0;
      }
    }
  }

  sum/=total;
  return sum;
}



/*************************************************************/
struct point *make_boundary_list(int x, int y, int v){

  int xcur,ycur;
  struct point *p;
  struct point *pstart;

  int xstart,ystart;
  int xtmp,ytmp;
  int n,nmax,imax,jmax;
  int xloop[]={-1, -1, -1, 0, 0, 1, 1, 1};
  int yloop[]={-1,  0,  1,-1, 1,-1, 0, 1};
  int i;
  int u;

  //printf("Boundary list for (%i,%i): ",x,y);
  //fflush(stdout);

  u=(y*xmax+x);
  //Make sure we're in the region we're interested in
  if(d[u]!=v) return NULL;

  p=point_malloc();
  p->prev=NULL; //To mark start of list
  pstart=p;

  //Go one pixel past the edge in any direction
  xcur=x;
  ycur=y;
  u=ycur*xmax;
  while(d[(u+xcur)]==v){
    xcur--;
    if(xcur<0){
      xcur++;
      break;
    }
  }
  //Now find maximum around this point for starting point
  nmax=-1;
  imax=-1;
  jmax=-1;
  for(i=0;i<8;i++){
    xtmp=xcur+xloop[i];
    ytmp=ycur+yloop[i];
    if((xtmp>=0)&&(xtmp<xmax)&&(ytmp>=0)&&(ytmp<ymax)){
      u=(ytmp*xmax+xtmp);
      if((d[u]!=v)||(xtmp==0)||(ytmp==0)||(xtmp==xmax-1)||(ytmp==ymax-1)){
	      //Search points that aren't in the region, or the boundary
	      n=neighboring_points(xtmp,ytmp,v);
	      if(n>nmax){
	        nmax=n;
	        imax=xtmp;
	        jmax=ytmp;
	      }
      }
    }
  }

  if(nmax>0){
    p->i=imax;
    p->j=jmax;
    p->next=point_malloc();
    (p->next)->prev=p;
    p=p->next; //Ready for next point
    u=(jmax*xmax+imax);
    work_array[u]=v;
    xcur=imax;
    ycur=jmax;
  }
  xstart=xcur;
  ystart=ycur;

  //Continue until we make a full loop back to the start, don't allow
  //previously used point to be used again.
  for(;;){

    nmax=-1;
    imax=-1;
    jmax=-1;
    //printf("(%i,%i) ",xcur,ycur);
    //fflush(stdout);
    for(i=0;i<8;i++){
      xtmp=xcur+xloop[i];
      ytmp=ycur+yloop[i];
      if((xtmp>=0)&&(xtmp<xmax)&&(ytmp>=0)&&(ytmp<ymax)){
	      u=(ytmp*xmax+xtmp);
	      if(work_array[u]!=v){
	        if((d[u]!=v)||(xtmp==0)||(ytmp==0)||(xtmp==xmax-1)||(ytmp==ymax-1)){
	          //Search points that aren't in the region, or the boundary
	          n=neighboring_points(xtmp,ytmp,v);
	          if(n>nmax){
	            nmax=n;
	            imax=xtmp;
	            jmax=ytmp;
	          }
	        }
	      }
      }
    }

    if(nmax>0){
      if((imax==(pstart->i))&&(jmax==(pstart->j))) break; //done loop
      p->i=imax;
      p->j=jmax;
      p->next=point_malloc();
      (p->next)->prev=p;
      p=p->next; //Ready for next point
      xcur=imax;
      ycur=jmax;
      u=(ycur*xmax+xcur);
      work_array[u]=v;
    }else{

      //Check if we've gone full circle
      if(((xcur-xstart)<=1)&&((xstart-xcur)<=1)&&((ycur-ystart)<=1)&&
	                                                     ((ystart-ycur)<=1)){
	      //printf("Full circle ");
	      //fflush(stdout);
	      break;
      }
      //We've gone down a tunnel.  We need to back up until we
      //have a point again.  We won't repeat this because we now
      //have set work_array[xcur][ycur]=v;
      if((p->prev)==NULL){
	      point_free(p);
	      //printf("Done boundary, no points\n");
	      //fflush(stdout);
	      return NULL; //No points
      }
      p=(p->prev);
      point_free(p->next);
      //p is now the next free location as it should be.
      //Set (xcur,ycur) to the one just before this one
      if((p->prev)==NULL){
	      point_free(p);
	      //printf("Done boundary !!!!, no points\n");
	      //fflush(stdout);
	      return NULL; //No points
      }
      xcur=(p->prev)->i;
      ycur=(p->prev)->j;
      //printf("Backup ");
      //fflush(stdout);

    }
  }

  if((p->prev)==NULL){
    point_free(p);
    //printf("Done boundary, SURPRISINGLY no points\n");
    //fflush(stdout);
    return NULL;
  }

  //Fix up dangling ends
  (p->prev)->next=NULL;
  point_free(p);

  //printf("Done boundary, return normal.\n");
  //fflush(stdout);
  return pstart;
}

/*************************************************************/
int neighboring_points(int x, int y, int v){
  //Count how many of the 8 neighboring points are equal to v.  If we leave
  //the image area, then count those non-pixels are not equal to v.

  //Just write out all the possibilities
  int n;

  int u1,u2,u3;

  n=0;

  u1=y*xmax+x;
  u2=u1-xmax;
  u3=u1+xmax;

  if((x+1)<xmax){
    if(d[u1+1]==v)n++;
  }
  if((x-1)>=0){
    if(d[u1-1]==v)n++;
  }

  if((y-1)>=0){
    if(d[u2]==v)n++;
    if((x+1)<xmax){
      if(d[u2+1]==v)n++;
    }
    if((x-1)>=0){
      if(d[u2-1]==v)n++;
    }
  }

  if((y+1)<ymax){
    if(d[u3]==v)n++;
    if((x+1)<xmax){
      if(d[u3]==v)n++;
    }
    if((x-1)>=0){
      if(d[u3]==v)n++;
    }
  }

  return n;
}

/*****************************************************************/
struct point *fix_spirals(struct point *p_begin){
  //Remove cells that kind of spiral inwards a few points but removing
  //points at the end of the list.  If the penultimate point is closer to
  //the first than the last, then we make that the last point.  Similarly
  //for the third to last, if the penultimate is closer than the last, etc.
  //Note we could just have easily have compared the last point to the first
  //instead of the first to the last, so there's some arbitrariness here.

  //I'm not sure this is such a good idea since it can remove an entire list
  //of cells.


  struct point *p;
  struct point *ptmp;
  float d1,d2;
  int a,b;

  a=(p_begin->i);
  b=(p_begin->j);

  p=p_begin;
  while((p->next)!=NULL){
    p=p->next;
  }
  d1=(float)( ((p->i)-a)*((p->i)-a) + ((p->j)-b)*((p->j)-b) );
  ptmp=p->prev;
  while(ptmp!=NULL){
    d2=(float)( ((ptmp->i)-a)*((ptmp->i)-a) + ((ptmp->j)-b)*((ptmp->j)-b) );
    if(d2<=d1){ //penultimate is closer than last point
      d1=d2;  //d1 is now distance of last point to start point
      ptmp->next=NULL; //Get rid of last point
      ptmp=ptmp->prev; //Try again
    }else{
      break;
    }
  }


  return p_begin; //Didn't adjust the start
}





/******************************************************/
struct point *clean_up_tails(struct point *p_begin){

  struct point *start;
  struct point *ptmp;
  struct point *p;
  int isect;

  float x,y; //do_segments_intersect() wants to return intersection
  //location

  start=p_begin;
  p=start;
  while((p->next)!=NULL){
    ptmp=(p->next)->next;  //Don't want adjacent segment so skip a point
    while(ptmp!=NULL){
      if((ptmp->next)!=NULL){
	isect=do_segments_intersect(
				    (float)p->i, (float)p->j,
				    (float)p->next->i, (float)p->next->j,
				    (float)ptmp->i, (float)ptmp->j,
				    (float)ptmp->next->i, (float)ptmp->next->j,
				    &x,&y);

	if(isect==1){ //Have intersection, get rid of all points before and
	  //including p, and all points after ptmp, but keep ptmp.

	  ptmp->next=NULL; //Memory leak! (but our weird allocation
	  //scheme makes it hard to clean this up, probably won't matter....)
	  start=p->next;
	  start->prev=NULL;

	  return start;
	}
      }
      ptmp=ptmp->next;
    }
    p=p->next;
  }
  //Note if some full-cells finish with last point being equal to first,
  //then the above loop will have removed the last point, which is a nice
  //feature.

  return start;

}

/*****************************************************/
struct point *find_interior_points(struct point *begin){
  //begin points to a linked list of points.  Here we find all the
  //points interior to the polygon that's formed by connecting all the
  //points pointed to by begin.  Create a linked list of interior points.

  struct point *p0;
  struct point *p1;

  int jmin,jmax;
  int j;
  int i,k;

  struct point *interior_start;
  struct point *r;

  float tmp;

  int isect_found;
  int isect;
  float x[isect_max];
  float y[isect_max];

  float xmin_use;  //=-5000.0;
  float xmax_use; //=5000.0;
  //Make bigger than xmax and make negative smaller than any point.
  //(Note that points can go off the screen because of fret offsets)
  xmax_use=xmax+10;
  xmin_use=-10;

  if(begin==NULL){
    printf("A NULL call sent to look for interior points!\n");
    return NULL;
  }
  if((begin->next)==NULL){
    printf("A 1-point call sent to look for interior points!\n");
    return NULL;
  }

  interior_start=point_malloc();
  r=interior_start;//Will be list of points (r points to next free location)
  r->prev=NULL;
  r->next=NULL; //In case want to break out from error in middle

  //Find min and max value of j in the point list
  jmin=ymax;
  jmax=0;
  for(p0=begin;p0!=NULL;p0=p0->next){
    if((p0->j)<jmin) jmin=p0->j;
    if((p0->j)>jmax) jmax=p0->j;
  }

  //We're going to loop from jmin to jmax and consider segments going
  //horizontally across entire picture.  We then look for intersections
  //of our line segments with this big segment.  We list the locations
  //of all found intersections.  If we write these points sorted in x as
  //i1,i2,i3,i4,i5,..., then the interior points are i1 to i2 and then
  //i3 to i4, etc.

  for(j=jmin;j<=jmax;j++){

    isect_found=0;
    for(p0=begin;p0!=NULL;p0=p0->next){
      p1=p0->next;
      if(p1==NULL)p1=begin; //Keep loop going around begin point

      isect=do_segments_intersect(
				  (float)p0->i,(float)p0->j,
				  (float)p1->i,(float)p1->j,
				  xmin_use,((float) j) + 0.1,
				  xmax_use,((float) j) + 0.1,
				  (x+isect_found),
				  (y+isect_found) );
      //Add .1 to j so don't have any questions about tangential
      //intersections.

      if(isect==1){
	      isect_found++;
	      if(isect_found>isect_max){
	        printf("Too many intersections found.\n");
	        isect_found--;
	      }
      }
    }//End loop over segment list for this value of j
    /*
    for(i=0;i<isect_found;i++){
      printf("       i-point: (%e,%e)\n",x[i],y[i]);
      fflush(stdout);
    }
    */

    //Make sure number found is even
    if((isect_found%2)!=0){
      printf("Odd number of intersection points!!!!!\n");
      fflush(stdout);
      point_list_free(interior_start); //Note last point points to NULL
      return NULL;
    }

    //If we've found at least two points, then sort them in x[] so can
    //calculate interior points.
    for(i=0;i<(isect_found-1);i++){
      for(k=i+1;k<isect_found;k++){
	if(x[i]>x[k]){
	  tmp=x[i];
	  x[i]=x[k];
	  x[k]=tmp;
	  tmp=y[i];
	  y[i]=y[k];
	  y[k]=tmp;
	}
      }
    }

    //Now add all the points between points 0 and 1,  2 and 3, etc.
    for(i=0;i<(isect_found-1);i+=2){
      for(k=(int)(x[i]+0.5);k<=((int)(x[i+1]+0.5));k++){
	r->i=k;
	r->j=j;
	r->next=point_malloc();
	(r->next)->prev=r;
	(r->next)->next=NULL;
	r=r->next;

      }
    }

  } //End of loop over j from jmin to jmax of segment points

  //r is an empty point, close up list and free r.
  if((r->prev)!=NULL){ //Make sure we found at least 1
    (r->prev)->next=NULL;
  }else{
    //If r->prev==NULL then there were no points, return NULL
    interior_start=NULL;
  }
  point_free(r);

  return interior_start; //Start of list of points

}

/*********************************************/
int do_segments_intersect(
			  float a1x, float a1y,
			  float b1x, float b1y,
			  float a2x, float a2y,
			  float b2x, float b2y,
			  float *x,  float *y){
  //Check if the segment connecting point a1 to b1 intersects the
  //segment connecting point a2 to b2.  Return 0 if don't intersect
  //and 1 if they do.
  //If they do intersect, then (*x,*y) will be the point of
  //intersection.  If the segments are parallel and also intersect (ie,
  //if they intersect in more than one point, then (x,y) will be
  //returned as a1.

  //We parameterize each segment as t*M+A where the capital letters are
  //2-vectors.  Using s for the a1-b1 line and t for the a2-b2 line, we want
  //to find s or t such that
  // s*M1+A1 = t*M2+A2 where M1=B1-A1 and M2=B2-A2
  //Multiplying by a vector perpendicular to M2 (call it M2') gives
  // s*(M2'*M1)=M2'*(A2-A1)
  //Which then gives s.  Similarly for t.
  //For our line parameterization, we set s=0 (t=0) to give one of the
  //endpoints, and s=1 (t=1) the other end point.  Then the segments
  //intersect iff 0<=s<1 and 0<=t<1.  (Note I only include the one
  //endpoint in the intersection.

  float m1x,m1y,m2x,m2y,a21x,a21y;
  float s,t;
  float m2_perp_m1;
  float m2_perp_a21;
  float m1_perp_a21;

  *x=-1.0;
  *y=-1.0;

  //Difference of intercept vectors (B2-B1) (We use a1 and a2 as start
  //points for parameterization)
  a21x=a2x-a1x;
  a21y=a2y-a1y;

  //Slope vectors
  m1x=b1x-a1x;
  m1y=b1y-a1y;
  m2x=b2x-a2x;
  m2y=b2y-a2y;
  //So t*M+B goves from point a1 to point b1 as t goes from 0 to 1
  //and similarly for points a2 and b2.

  m2_perp_m1=(m2y*m1x-m2x*m1y);
  m2_perp_a21=(m2y*a21x-m2x*a21y);
  m1_perp_a21=(m1y*a21x-m1x*a21y); //Note m1_perp_m2 = -m2_perp_m1

  //Calculate s and t
  if(m2_perp_m1==0.0){ //Lines are parallel
    if(m2_perp_a21==0.0){ //lines are the same
      //Check if A2 or B2 is located on the 1 segment somewhere for this
      //special case.
      if(m1x!=0.0){
	s=a21x/m1x; //Location of A2 on 1 line
      }else if(m1y!=0.0){
	s=a21y/m1y; //Location of A2 on 1 line
      }else{
	return 0;  //Slope vector is 0!
      }
      if(m1x!=0.0){
	t=(b2x-a1x)/m1x; //Location of b2 on 1 line
      }else if(m1y!=0.0){
	t=(b2y-a1y)/m1y; //Location of B2 on 1 line
      }else{
	return 0;  //Slope vector is 0!
      }
      //If segments are on same line then we have an intersection if
      //either A2 or B2 is in the A1-B1 segment--ie, if either s or t as
      //defined above is between 0 and 1.
      if(((0.0<=s)&&(s<1.0))||((0.0<=t)&&(t<1.0))){
	*x=a1x;
	*y=a1y;
	return 1;
      }else{
	return 0;
      }
    }else{
      return 0; //lines are parallel and are different lines so
      //no intersection
    }

  }else{ //Now do the most common, non-parallel case
    s=m2_perp_a21/m2_perp_m1;
    t=m1_perp_a21/m2_perp_m1;
  }

  if((0.0<=s)&&(s<1.0)&&(0.0<=t)&&(t<1.0)){
    *x=a1x+m1x*s;
    *y=a1y+m1y*s;
    return 1;
  }else{
    return 0;
  }

}

/***********************************************************/
void remove_overlaps(void){
  //Compare the interior points of the different cells, and if two cells
  //overlap too much, remove the smaller.

  int a;
  int i,j;
  int remove_point;

  struct point *istart;
  int isection;

  float tmp;
  float overlap_cut=0.3;

  // float circ,area;
  //  float c_over_a1,c_over_a2;

  for(i=n_before_fret_copy;i<(n_found-1);i++){

    istart=interior[i]; //Use to clear the array below since i index

    /*
    p=istart;
    while(p!=NULL){
      printf("(%i,%i)(p=%i,n=%i,c=%i)\n",
	     p->i,p->j,p->prev,p->next,p);
      p=p->next;
    }
    */

    /*
    //Calculate circumference**2/area for comparison below
    area=area_of_cell(boundary[i]);
    circ=circumference_of_cell(boundary[i]);
    if(area!=0.0){
      c_over_a1=circ*circ/area;
    }else{
      c_over_a1=1.0e30;
    }
    */

    //might change when we replace cells.
    update_overlap_value();
    fill_overlap_array_with_point_list(istart);
    for(j=i+1;j<n_found;j++){
      isection=overlap(interior[j],0,0);//Compare jth with ith point

      if(isection>0){
       	if(n_points[i]<n_points[j]){
	  remove_point=i;
	}else{
	  remove_point=j;
	}
	tmp=( (float) isection ) / ( (float) n_points[remove_point] );

	if(tmp>overlap_cut){
	  //Remove cell with fewer interior points.
	  //Do this by putting last cell in this location and decrementing
	  //n_found.

	  /*
	  //Actually sometimes a real fucked up cell overlaps a lot of
	  //other ones.  So let's remove the more circular cell, as
	  //measured by the smaller value of circumference^2/area
	  area=area_of_cell(boundary[j]);
	  circ=circumference_of_cell(boundary[j]);
	  if(area!=0.0){
	    c_over_a2=circ*circ/area;
	  }else{
	    c_over_a2=-1.0;
	  }
	  if(c_over_a2>c_over_a1){
	    remove_point=j;
	  }else{
	    remove_point=i;
	  }
	  printf(
		 "O: (%i,%i)(N=%i,c2/a=%e) (%i,%i)(N=%i,c2/a=%e),I=%i,f=%e\n",
		 (int)(mean_x[i]+.5),(int)(mean_y[i]+.5),n_points[i],
		 c_over_a1,
		 (int)(mean_x[j]+.5),(int)(mean_y[j]+.5),n_points[j],
		 c_over_a2,
		 isection,tmp);
	  fflush(stdout);
	  */

	  a=n_found-1;
	  n_found--;
	  point_list_free(interior[remove_point]);
	  point_list_free(boundary[remove_point]);
	  interior[remove_point]=interior[a];
	  boundary[remove_point]=boundary[a];
	  n_points[remove_point]=n_points[a];

	  //mean_x[remove_point]=mean_x[a];
	  //mean_y[remove_point]=mean_y[a];
	  if(remove_point==i){
	    i--;
	    break; //Break out of j loop to start again with new i
	  }else{
	    j--; //Do this j again since we just replaced it
	  }
	} //end if past overlap_cut
      }//end if are any overlap points (isection>0)

    }//End loop over inner loop of bubble sort (j)

  }//End loop over outer loop over bubble sort (i)

  return;
}

/***********************************************************/
void update_overlap_value(void){

  int u;

  //Make it so overlap_value is never 0 and initialize array to be 1.
  //This way 0 can be used to consistently elsewhere to mark bad points
  //or some such thing.
  if(over==NULL){
    over=(unsigned char *)malloc(over_array_max*sizeof(unsigned char));
    overlap_value=1;
    for(u=0;u<over_array_max;u++) over[u]=overlap_value;
  }

  overlap_value++;
  if (overlap_value>overlap_value_max){
    overlap_value=2;
    for(u=0;u<over_array_max;u++) over[u]=1; //Reset array

  }
  return;
}

/***********************************************************/
void fill_overlap_array_with_point_list(struct point *begin){

  int u;
  struct point *r;
  int ix,iy;

  r=begin;
  while(r!=NULL){
    ix=(r->i);
    iy=(r->j);
    if ((ix>=0)&&(ix<xmax)&&(iy>=0)&&(iy<ymax)){
      u=(iy*xmax+ix);
      over[u]=overlap_value;
    }else{
      //printf("u too high!!!! (%i,%i).\n",r->i,r->j);
    }
    r=r->next;
  }

  return;
}

/***********************************************************/
void fill_overlap_array_with_point_list_offset(struct point *begin,
					       int offx,
  					       int offy){

  int u;
  struct point *r;
  int ix,iy;

  r=begin;
  while(r!=NULL){
    ix=(r->i)+offx;
    iy=(r->j)+offy;
    if ((ix>=0)&&(ix<xmax)&&(iy>=0)&&(iy<ymax)){
      u=(iy*xmax+ix);
      over[u]=overlap_value;
    }else{
      //printf("u too high!!!! (%i,%i).\n",r->i,r->j);
    }
    r=r->next;
  }

  return;
}

/***********************************************************/
int overlap(struct point *begin, int offset_i, int offset_j){
  //Calculate the intersection of two linked lists of points.
  //The first point should already have been written into the over[] array
  //using fill_overlap_array_with_point_list().
  //this since this routine will usually be used in the middle of a
  //bubble sort, so don't want to waste time refilling and clearing memory.

  int u;

  int n;
  int ix,iy;
  struct point *r;

  r=begin;
  n=0;
  while(r!=NULL){
    ix=(r->i)+offset_i;
    iy=(r->j)+offset_j;
    if ((ix>=0)&&(ix<xmax)&&(iy>=0)&&(iy<ymax)){
      u=iy*xmax+ix;
      if(over[u]==overlap_value) n++;
    }else{
      //Can happen without an error if the offsets aren't 0
      //printf("u too high in comparison!!!! (%i,%i,%i).\n",r->i,r->j,u);
    }
    r=r->next;
  }

  return n;
}
/***********************************************************/
int output_individual_cells_to_file(int i_t,
				    char *basefile,
				    float *input_data,
				    int xmax_data,
				    int ymax_data,
				    int type,
				    int bit_size,
				    int invert){


  int xmax_out=50;
  int ymax_out=50;
  float *output_data=NULL;
  int *output_labels=NULL;

  char *file=NULL;
  int filelen=1000;

  int i;
  struct blob *b;
  int i0,j0;
  char cellnum[100];

  int ix,iy,u;
  int iuse,juse,uuse;
  struct point *p;

  file=(char *)malloc(filelen*sizeof(char));
  output_data=(float *)malloc(xmax_out*ymax_out*sizeof(float));
  output_labels=(int *)malloc(xmax_out*ymax_out*sizeof(int));

  for(i=0;i<n_known;i++){
    b=cs[i]; //cs[] points to most recently added cell
    if(b->i_time==i_t){
      //Write number i at position (b->x,b->y)
      //      printf("Adding number %i at (%e,%e).\n",b->index,b->x,b->y);
      //fflush(stdout);

      if (((b->x)>=0.0)&&((b->y)>=0.0)){
	i0=((int)(b->x))-xmax_out/2;
	j0=((int)(b->y))-ymax_out/2;

	//Zero out this box first (ie, mark all pixels as "deleted")
	for(ix=0;ix<xmax_out;ix++){
	  for(iy=0;iy<ymax_out;iy++){
	    u=(iy*xmax_out+ix);
	    output_labels[u]=delete_pixel;
	    output_data[u]=0.0;
	  }
	}

	//Add in our cell and its boundary
	for(p=b->interior;p!=NULL;p=p->next){
	  ix=(p->i)-i0;
	  iy=(p->j)-j0;
	  if ((ix>=0)&&(ix<xmax_out)&&(iy>=0)&&(iy<ymax_out)){
	    u=(iy*xmax_out+ix);
	    output_labels[u]=0; //transparent
	  }
	}
	for(p=b->boundary;p!=NULL;p=p->next){
	  ix=(p->i)-i0;
	  iy=(p->j)-j0;
	  if ((ix>=0)&&(ix<xmax_out)&&(iy>=0)&&(iy<ymax_out)){
	    u=(iy*xmax_out+ix);
	    output_labels[u]=found_border;
	  }
	}

	//Now copy over to our box
	//Zero out this box first (ie, mark all pixels as "deleted")
	for(ix=0;ix<xmax_out;ix++){
	  for(iy=0;iy<ymax_out;iy++){
	    //Now zero out d[] array in this region (see we can
	    //use add_boundary_pixels()....
	    u=(iy*xmax_out+ix);
	    iuse=ix+i0;
	    juse=iy+j0;
	    output_data[u]=0.0;
	    if ((iuse>=0)&&(iuse<xmax_data)&&(juse>=0)&&(juse<ymax_data)){
	      uuse=(juse*xmax_data+iuse);
	      output_data[u]=input_data[uuse];
	    }
	  }
	}

	//Put number in middle of bottom
	add_numbers_to_data(b->index,
			    (xmax_out/2),(ymax_out-10),
			    output_labels,
			    xmax_out,ymax_out);

	strcpy(file,basefile);
	strcat(file,"_");
	digits_to_string(cellnum,bit_size,10);
	strcat(file,cellnum);
	strcat(file,"b_id_");
	digits_to_string(cellnum,b->index,9999);
	strcat(file,cellnum);
	strcat(file,".tif");
	if(output_data_to_tif_file(file,
                             output_data,
                             xmax_out,
                             ymax_out,
                             output_labels,
                             type,
                             bit_size,
                             invert,
                             0)==0){
	  printf("Couldn't output cell %i to %s.\n",b->index,file);
	}

      }
    }

  }

  free(file);
  free(output_data);
  free(output_labels);


  return 1;
}






/***********************************************************/
void background_level(int i_time){
  //Calculate the mode of all the pixels that aren't found in the
  //interior of any cells.  Also calculate the total fluorescence in
  //all cells that are found.

  float h[nbins_for_cut_calculation];
  float w;

  float min_pixel_value;
  float max_pixel_value;

  float tmp;
  int i,j;
  int k;
  int u;

  float mu;
  float mu1,mu2,mu3;

  int nbins;
  int nbins_max;
  float bin_width;
  float inv_bin_width;

  int i1,i2,i3;
  float hmax1,hmax2,hmax3;
  float max_bin_width;

  int max_calc_done;

  int fret_region;

  nbins_max=nbins_for_cut_calculation;

  //  double back_calc;

  //Try a new method of averaging the largest contiguous region that's
  //not a border.  This is calculated in find_cells() above and is
  //in the variable back_average.
  /*
  if(n_pixels_for_background>0){
    back_calc=0.0;
    for(k=0;k<n_pixels_for_background;k++){
      i=(int)(pixels_for_background_x[k]);
      j=(int)(pixels_for_background_y[k]);
      back_calc+=(double)(fl[i][j]);
      printf("%i %i %e %e\n",i,j,fl[i][j],c[i][j]);fflush(stdout);
    }
    back_calc/=((double)n_pixels_for_background);
    back_pixels[i_time]=(float)(back_calc);;
    exit(0);
    return;
  }
  */

  //If no pixles found from above, try the old method of finding the
  //mode of the unused cells.

  back_pixels[i_time]=0.0;
  back_pixels_1[i_time]=0.0;
  back_pixels_2[i_time]=0.0;

  //First add all known cells to overlap array so we can remove them
  //from the mode calculation
  update_overlap_value();
  for(i=0;i<n_found;i++){//was "<n_known" before 03/12/04-which is
    //wrong. interior[] is only defined for n_found.
    fill_overlap_array_with_point_list(interior[i]);
  }

  //Get max and min pixel value for histogramming
  //Only use pixels that aren't in any other cells
  if (have_fret_image==1){
    if (fret_labels==NULL){
      printf("No split-labels for split image in background");
      printf(" calculation.\n");
      perror(0); //exit(0); //http://r-pkgs.had.co.nz/src.html
    }
    fret_region=1; //The upper part of the image (lower y) is labelled 2
  }

 start_fret_goto:

  min_pixel_value=1.0e30;
  max_pixel_value=-1.0;
  if (have_fret_image==1){
    for(u=0;u<xmax_ymax;u++){
      if (fret_labels[u]==fret_region){
	if(over[u]!=overlap_value){
	  tmp=fl[u];
	  if (tmp<=0.1) continue;
	  if(tmp<min_pixel_value) min_pixel_value=tmp;
	  if(tmp>max_pixel_value) max_pixel_value=tmp;
	}
      }
    }
  }else{
    for(u=0;u<xmax_ymax;u++){
      if(over[u]!=overlap_value){
	tmp=fl[u];
	if (tmp<=0.1) continue;
	if(tmp<min_pixel_value) min_pixel_value=tmp;
	if(tmp>max_pixel_value) max_pixel_value=tmp;
      }
    }
  }

  if(min_pixel_value==max_pixel_value){
    printf("min=%e=%e=max, returning cut value of zero. \n",
	   min_pixel_value,max_pixel_value);
    return;
  }

  max_bin_width=(max_pixel_value-min_pixel_value)/10.0;

  //width of bins has to be an integer since fluorescence values
  //are actually integers from ccd camera.
  for(i=1;;i++){
    w=(max_pixel_value-min_pixel_value)/((float)i)+1.0; //~number of bins
    if (w<((float)nbins_max)) break; //is ok
  }
  bin_width=((float)i);

  back_pixels[i_time]=-1.0e6; //Initialization just in case can't find
  //adequate bin size.
  max_calc_done=0;
  while(max_calc_done==0){
    //Will increase bin_width in loop if think too small

    inv_bin_width=1.0/bin_width;
    w=min_pixel_value;
    for(i=0;w<=max_pixel_value;i++){
      h[i]=0.0;
      w+=bin_width;
    }

    if (have_fret_image==1){
      for(u=0;u<xmax_ymax;u++){
	if(fret_labels[u]==fret_region){
	  if(over[u]!=overlap_value){
	    k=(int)(  (fl[u] - min_pixel_value)*inv_bin_width );
	    if((k>=0)&&(k<nbins_for_cut_calculation)){
	      h[k]+=1.0;
	    }
	  }
	}
      }
    }else{
      for(u=0;u<xmax_ymax;u++){
	if(over[u]!=overlap_value){
	  k=(int)(  (fl[u] - min_pixel_value)*inv_bin_width );
	  if((k>=0)&&(k<nbins_for_cut_calculation)){
	    h[k]+=1.0;
	  }
	}
      }
    }

    //Find highest three bins
    hmax1=-1.0;
    hmax2=-1.0;
    hmax3=-1.0;
    i1=0;
    i2=5;
    i3=10;//Initialize far apart just in case.
    w=min_pixel_value;
    nbins=0;
    for(i=0;w<=max_pixel_value;i++){
      //printf("%i %e\n",i,h[i]);fflush(stdout);
      tmp=h[i];
      if(tmp>hmax1){ //Find top three bins
	hmax3=hmax2;
	hmax2=hmax1;
	hmax1=tmp;
	i3=i2;
	i2=i1;
	i1=i;
      }else if (tmp>hmax2){
	hmax3=hmax2;
	hmax2=tmp;
	i3=i2;
	i2=i;
      }else if (tmp>hmax3){
	hmax3=tmp;
	i3=i;
      }
      w+=bin_width;
      nbins++;
    }

    //Now check that the top three bins are all near each other. If
    //not then bin width might have been too low
    i=i1-i2;
    if (i<0) i=(-i);
    j=i2-i3;
    if (j<0) j=(-j);
    k=i1-i3;
    if (k<0) k=(-k);
    if ((i>1)||((j>1)&&(k>1))){ //Make bin_width bigger
      bin_width+=1.0;
      if(bin_width>=max_bin_width){
	max_calc_done=2; //Break out of loop but flag error
      }
    }else{
      max_calc_done=1; //Done calculation, calculation ok
    }
  } //End while(max_calc_done==0)

  if (max_calc_done==2){
    printf("***** Couldn't get good mode calculation for background");
    printf("pixels. *****\n");
    back_pixels[i_time]=-1.0e6;
    if (have_fret_image==1){
      printf("***** (for split-image region=%i\n",fret_region);
      if (fret_region==1){
	back_pixels_1[i_time]=-1.0e6;
      }else{
	back_pixels_2[i_time]=-1.0e6;
      }
    }
  }else{ //mode calculation went ok
    //Average the three top locations
    mu1=min_pixel_value+( (((float)i1)+0.5)*bin_width );
    mu2=min_pixel_value+( (((float)i2)+0.5)*bin_width );
    mu3=min_pixel_value+( (((float)i3)+0.5)*bin_width );

    mu=(mu1*hmax1+mu2*hmax2+mu3*hmax3)/(hmax1+hmax2+hmax3);

    back_pixels[i_time]=mu;
    back_cur=mu;
    if (have_fret_image==1){
      if (fret_region==1){
	back_pixels_1[i_time]=mu;
      }else{
	back_pixels_2[i_time]=mu;
      }
    }
  }
  if (have_fret_image==1){
    printf("Background in non-cells at time t=%i for region %i is %e.\n",
	   i_time,fret_region,back_pixels[i_time]);
    if (fret_region==1){
      fret_region=2;
      goto start_fret_goto; //To do next fret region
    }
  }else{
    printf("Background in non-cells at time t=%i is %e (bin_width=%e).\n",
	   i_time,back_pixels[i_time],bin_width);
  }

  return;
}
/***********************************************************/
void calculate_cut(){
//Calculates cut for the intensity level when finding cell borders.


  //float h[nbins_for_cut_calculation];
  float w;

  float min_pixel_value=1.0e30;
  float max_pixel_value=0.0;

  float tmp;
  int i;
  //int k;

  float mu,sig,scale;

  double mean_other,rms_other;
  int n_other;
  float cut_other;

  //int nbins;

  cut=0.0; //In case don't find anything
  cut_low=0.0;
  cut_high=0.0;

  if((have_fret_image==1)&&(fret_labels==NULL)){
    printf("Haven't calculated upper and lower regions from");
    printf("fluorescence image.\n");
    fflush(stdout);
    perror(0); //exit(0); //http://r-pkgs.had.co.nz/src.html
  }

  //Get max and min pixel value for histogramming
  mean_other=0.0;
  rms_other=0.0;
  n_other=0;
  cut_other=-1.0;
  if ((image_type==fret_bf_bottom_only)||
      (image_type==fret_bf_top_only)){
    //Calculate other region first for a cut value (do this since if
    //it's only the bottom part, for exapmle, the upper can be very low
    //compared to the bottom part and that can screw up the histogram
    //below, and vice versa.
    for(i=0;i<xmax_ymax;i++){
      tmp=c[i];
      if (fret_labels[i]!=fret_region_use){ //other region (upper if
	     //looking at lower, lower if looking at upper)
	     mean_other+=tmp;
	     rms_other+=(tmp*tmp);
	     n_other++;
      }
    }
    mean_other/=((double)n_other);
    rms_other=(rms_other/((double)n_other))-(mean_other*mean_other);
    if (rms_other<=0.0){
      printf("Negative rms of opposite region: %e\n",rms_other);
      cut_other=-1.0;
    }else{
      rms_other=sqrt(rms_other);
      //Will only consider pixels in lower region above 3*upper region
      //(Since the exact division of lower and upper might have some
      //stray pixels)
      cut_other=((float)(mean_other+3.0*rms_other));
      printf("Cut from opposite region in cut calculation: %e\n",cut_other);
    }
  }

  //Add cut_other as a requirement on the lower_fret_region
  //for
  mu=0.0;
  sig=0.0;
  scale=0.0;

  //Loop over all the pixels of the bright field image array
  for(i=0;i<xmax_ymax;i++){
    tmp=c[i];
    if(have_fret_image==1){
      if (fret_labels[i]!=fret_region_use){
	      continue;
      }
      if (tmp<cut_other) continue; //for FRET case only
    }
    if (tmp<=0.1) continue;
    mu+=tmp;
    sig+=(tmp*tmp);
    scale+=1.0;
    if(tmp<min_pixel_value) min_pixel_value=tmp;
    if(tmp>max_pixel_value) max_pixel_value=tmp;
  }

  mu/=scale; //mu=sum(pix)/nPix
  sig/=scale;//sig=sum(pix^2)/nPix
  sig=(float)sqrt((double)(sig-mu*mu));
  //sig=( sum(pix^2)/nPix - (sum(pix)/nPix)^2 )^1/2
  printf("mean and sig: %e,%e\n",mu,sig);
  //  for (i=0;i<xmax_ymax;i++){
  // printf("%e\n",c[i]);
  //}
  //exit(0);
  //Just use mean and sig and skip gaussian fit.--4/23/04
  mu-=min_pixel_value;
  w=1.0;

 /*
  goto skip;

  //Most points are background, so make a histogram of all points, and
  //the data should show up as a high end tail.  Find it by looking for
  //first crossing of zero of the derivative from it's most negative point.
  //(Ie, the data is peaked around the background value, the falling part
  //of that is the most negative point for the derivative.  We're going to
  //put the cut just where the derivative first crosses zero).


  printf("Min and max=%e, %e.\n",min_pixel_value,max_pixel_value);
  if (min_pixel_value<0.5){ //Don't allow a 0.0 min value in histogram.
    //This is because if we're fitting a metamorph-"deconvolution"
    //image then there's a big spike at 0 which we want to be our
    //zero.
    if((image_type!=bright_field) && (image_type!=confocal_transmission) ){
      printf("!!!!! Found minimum pixel value of 0, but not doing\n");
      printf("!!!!! metamorph deconvolution case.\n");
    }
  }


  if(min_pixel_value==max_pixel_value){
    printf("min=%e=%e=max, returning cut value of zero. \n",
	   min_pixel_value,max_pixel_value);
    return;
  }

  nbins=nbins_for_cut_calculation;

 gauss_fit_start:

  w=((float) nbins)/(max_pixel_value-min_pixel_value);

  for(i=0;i<nbins;i++){
    h[i]=0.0;
  }

  for(i=0;i<xmax_ymax;i++){
    if(have_fret_image==1){
      if (fret_labels[i]!=fret_region_use) continue;
    }
    if (c[i]>=min_pixel_value){
      k=(int)(  (c[i] - min_pixel_value)*w );
      if((k>=0)&&(k<nbins)){
	h[k]+=1.0;
      }
    }
  }

  //Find max bin and height (put in mu and scale)
  mu=0.0;
  tmp=0.0;
  sig=0.0;
  scale=0.0;
  for(i=0;i<nbins;i++){
    if(h[i]>scale){
      scale=h[i];
      mu=(float) i;
    }
  }
  //
  //if (scale<100.0){
  //  nbins=nbins/10;
  //  if (nbins>10){
  //    goto gauss_fit_start;
  //  }
  //}


  //Fit the peak to a gaussian.  Use mu=peak location and scale=size of
  //peak for first guesses.  Guess sig to be something reasonable.
  if((image_type==bright_field)||(image_type==confocal_transmission)
                               ||(have_fret_image==1)){
    sig=10.0;
  }else if(image_type==metamorph_deconvolution){
    sig=1.0;
  }

  printf("Starting values for gauss fit: mu=%e, sig=%e, scale=%e.\n",
	 mu,sig,scale);fflush(stdout);
  k=gauss_fit(h, nbins, &mu, &sig, &scale);
  if(k==0){
    printf("!!!!!!! A problem doing gaussian fit to get cut value.\n");
    printf("!!!!!!! mu=%e, sig=%e, scale=%e\n",mu,sig,scale);
    nbins=nbins/2;
    if (nbins<10){
      printf("Couldnt get good gaussian fit.\n");
    }else{
      goto gauss_fit_start;
    }
  }
  if ((mu<0.0)||(sig<=0.0)){
    printf("---------->(%e,%e,%e)\n",mu,sig,max_pixel_value);fflush(stdout);
    max_pixel_value*=.8; //Maybe have an outlier so divide and do again
    if ((max_pixel_value<min_pixel_value)||(max_pixel_value<1.0)){
      printf("Couldnt get good gaussian fit in calculate_cut.\n");
      mu=min_pixel_value*w;
      sig=5.0*w;
    }else{
      goto gauss_fit_start;
    }
  }
 skip:
  */

  printf("In Pixel units: mu=%e, sig=%e.\n",min_pixel_value+mu/w,sig/w);

  //Set cut to be 3*sigma above mean.  mu and sig are in units of
  //histogram bin numbers.  Change them to be pixel size units
  //cut=min_pixel_value+( (mu+2.*sig)/w );
  cut=min_pixel_value+( (mu+3.0*sig)/w );

  //cut-values;
  cut_low=min_pixel_value+( (mu-background_reject_factor*sig)/w );
  cut_high=min_pixel_value+( (mu+background_reject_factor*sig)/w );

  if (cut_low<min_pixel_value){
    cut_low=min_pixel_value;
    cut_high=min_pixel_value+( 3.0*sig/w );
  }
  mid=(cut_high+cut_low)/2.0;

  return;
}


/**********************************************************/
/*
void check_mem(void){

  int i;
  for(i=0;i<100;i++){
    if((pmem+i)->i!=7777){
      printf("EARLY OVERWRITE!!!!! %i %i.\n",i,(pmem+i)->i);
    }
  }
  for(i=mem_size;i<mem_size+100;i++){
    if((pmem+i)->i!=8888){
      printf("LATE OVERWRITE!!!!! %i %i.\n",i,(pmem+i)->i);
    }
  }

  i=max_cells-1;

  if(boundary[i]!=NULL) printf("Screwup in boundary[](%i,%i)!!!!!\n",
			       i,(int) boundary[i]);
  if(interior[i]!=NULL) printf("Screwup in interior[](%i,%i)!!!!!\n",
			       i,(int) interior[i]);
  if(mean_x[i]!=333.0) printf("Screwup in x[](%i,%e)!!!!!\n",i,mean_x[i]);
  if(mean_y[i]!=222.0) printf("Screwup in y[](%i,%e)!!!!!\n",i,mean_y[i]);
  if(n_points[i]!=111) printf("Screwup in n_points(%i,%i)!!!!!\n",
			      i,n_points[i]);

  return;
}
*/

/**********************************************************/
struct point *point_malloc(void){
  //My own allocation system to keep things fast since I think I call
  //malloc a lot otherwise.

  //We need one more sizeof(struct point) bit of memory.  Next free
  //location in the point_mem array is pointed at by pmem_cur.
  struct point *p;
  struct point_mem *ptmp;
  int i;
  int j;

  //This system assumed that you would de-allocate points in reverse order
  //of their allocation, but I don't think that's really true anymore.
  //So just do malloc() and free().
  p=(struct point *)malloc(sizeof(struct point));
  return p;


  if(pmem_start==NULL){ //first time through
    pmem_start=(struct point_mem *)malloc(sizeof(struct point_mem));
    pmem_start->next=NULL;
    pmem_cur=pmem_start;
    i=(int)pmem_start;
    j=(int)(&(pmem_start->p)); //Make sure these are the same, otherwise
    //need to calculate an offset for point_free().
    if(i!=j){
      printf("Need to introduce an offset into point_free() to since\n");
      printf("the structure and the memory-unit holding the structure\n");
      printf("have different addresses (%i != %i) in point_malloc().\n",i,j);
      perror(0); //exit(0); //http://r-pkgs.had.co.nz/src.html
    }
  }

  if(pmem_cur->next==NULL){  //Allocate new block of memory
    ptmp=(struct point_mem *)malloc((mem_size)*sizeof(struct point_mem));

    //Set up the list
    for (i=0;i<mem_size-1;i++){
      ptmp[i].next=(ptmp+i+1);
    }
    ptmp[mem_size-1].next=NULL; //Last in list

    pmem_cur->next=ptmp;

    if(next_block>=n_mem_blocks){
      printf("Giant memory leak!  Giving up keeping track of large");
      printf("blocks of memory!");
      next_block=0; //start over
    }

    mem_blocks[next_block]=pmem_cur; //Save all the allocated blocks so can
    //delete completely later.
    next_block++;  //Location to put next large block of memory
  }

  p=&(pmem_cur->p); //Address of "struct point" in the structure
  pmem_cur=pmem_cur->next;

  return p;
}

/**********************************************************/
void point_list_free(struct point *p_in){

  struct point *p;

  if(p_in==NULL)return;
  for(p=p_in;(p->next)!=NULL;p=(p->next)){
    point_free(p->prev);
  }
  point_free(p->prev);
  point_free(p);
  return;
}
/**********************************************************/
void point_free(struct point *p){

  struct point_mem *ptmp;


  //This system assumed that you would de-allocate points in reverse order
  //of their allocation, but I don't think that's really true anymore.
  //So just do malloc() and free().
  free(p);
  return;

  //Here's a hack--I'm going to assume that the address of p is also
  //the address (or at least a constant offset). (Maybe could add an
  //offset if this doesn't work.)
  if(p==NULL)return;

  //printf("f(%i)--",(int)p);
  ptmp=pmem_cur;
  pmem_cur=(struct point_mem *)p; //Make the new first location

  pmem_cur->next=ptmp;  //point to previous first location

  return;
}

/**********************************************************/
void update_list_of_found_cells(int i_t, int secs, int flag){
  //Go through all the known cells and compare to the currently found
  //cells.  If a current cell isn't near any of them, then start a
  //new one in the list of known cells, otherwise add the current one to
  //the known list using time index i_t.
  //secs is the time in seconds since the first image.

  struct blob *b;
  struct blob *bnew;
  struct blob *bsave;
  int i,j;
  int iloc;

  int loop_total;

  float I_over_U; //Intersection over union
  float I_over_U_max;
  //float c_over_a1;
  //float c_over_a2;
  //float area,circ;
  int j_max;
  int id_fret;
  int n_p;
  int isection;
  int uon;

  int current_total;
  int offset_i,offset_j;
  int tmp_i,tmp_j;

  int add_new_cell;
  float area;

  printf("Number found: %i\n",n_found);

  total_time++; //Keep track of how many pictures we've checked
  loop_total=n_known;
  //Use loop_total because n_known will change as we add cells, and we
  //only want to compare the cells we just found with cells we found in
  //previous times.  We don't want to compare cells from the same times.

  offset_i=0;
  offset_j=0;

  if(new_phase==1){ //Have a new image, do overlap search
    //Compare each current cell too all the known cells

    printf("Comparing new cells to known cells.\n");

    //Calculate an offset for the current image with the previous one

    //Prepare for calculation over intersection with all known points by
    //filling the over[] array with all the interior points from the
    //previous image
    i=0;
    if (total_time>1){ //Was a previous image at all
      //Note that we just incremented total_time above
      update_overlap_value();
      for(j=0;j<loop_total;j++){
	b=cs[j];  //cs[j] is the previous time point found for cell j
	if (b->x<=0) continue; //As a flag
	if((b->i_time)==(i_t-1)){ //We only want the previous image
	  i++;
	  fill_overlap_array_with_point_list(b->interior);
	}
      }
      //The overlap array now contains all the points found from the
      //previous image.  We now vary an offset to maximize the total overlap
      //with the current image.
    }
    if(i>0){ //If we had a previous image and there were cells
      current_total=0;

      if(image_type!=hexagonal_grid){
	//Loop over all possible overlaps, start with big steps and focus in.
	//Don't do this for hexagonal grid since it's regularity could lead
	//to false, large offsets.
	for(tmp_i=-250;tmp_i<250;tmp_i+=10){
	  //printf("%i: (%i,%i)=%i.\n",tmp_i,offset_i,offset_j,current_total);
	  for(tmp_j=-250;tmp_j<250;tmp_j+=10){
	    if((j=total_overlap_of_all_cells(tmp_i,tmp_j))>current_total){
	      current_total=j;
	      offset_i=tmp_i;
	      offset_j=tmp_j;
	    }
	  }
	}
      }

      tmp_i=offset_i;
      tmp_j=offset_j;

      do{
	offset_i=tmp_i;
	offset_j=tmp_j;
	while(
	      (isection=total_overlap_of_all_cells(tmp_i+1,tmp_j))>
	      current_total ){
	  current_total=isection;
	  tmp_i++;
	}
	if(tmp_i==offset_i){ //Didn't gain anything going right, check left
	  while(
		(isection=total_overlap_of_all_cells(tmp_i-1,tmp_j))>
		current_total ){
	    current_total=isection;
	    tmp_i--;
	  }
	}
	while(
	      (isection=total_overlap_of_all_cells(tmp_i,tmp_j+1))>
	      current_total ){
	  current_total=isection;
	  tmp_j++;
	}
	if(tmp_j==offset_j){ //Didn't gain anything going right, check left
	  while(
		(isection=total_overlap_of_all_cells(tmp_i,tmp_j-1))>
		current_total ){
	    current_total=isection;
	    tmp_j--;
	  }
	}
      } while((tmp_i!=offset_i)||(tmp_j!=offset_j)) ;

    } //If there was a previous image to calculate offset against

    //Below we put the current image in the over[] array,
    //so meaning of offsets gets reversed.
    offset_i=-offset_i;
    offset_j=-offset_j;

    printf("Offset from previous image: (%i,%i).\n",offset_i,offset_j);
    //We'll only use the offsets to calculate overlaps to match to cells
    //from previous image.

  }else{
    I_over_U_max=999.0; //Always will pass cut below for finding a match
  }

  total_area[i_t]=0.0;
  total_fluorescence[i_t]=0.0;
  for(i=0;i<n_found;i++){
    total_fluorescence[i_t]+=fluorescence[i];
    area=cell_area[i];
    total_area[i_t]+=area;
    //Create a new "blob" pointer to put this cell data in.  Below we'll
    //decide where in the cs[] array of linked lists to put it.
    bnew=(struct blob *)malloc(sizeof(struct blob));
    bnew->next=NULL; //Will be last in the linked list one way or other
    bnew->x=mean_x[i];
    bnew->y=mean_y[i];
    bnew->circumference=circumference[i];
    bnew->vol_rotation=vol_rotation[i];
    bnew->vol_cone=vol_cone[i];
    bnew->vol_sphere=vol_sphere[i];
    bnew->surface_area=surface_area[i];
    bnew->vol_eff_1=vol_eff_1[i];
    bnew->vol_eff_2=vol_eff_2[i];
    bnew->vol_eff_3=vol_eff_3[i];
    bnew->vol_eff_4=vol_eff_4[i];
    bnew->vol_eff_5=vol_eff_5[i];
    bnew->vol_eff_6=vol_eff_6[i];
    bnew->major_axis_length1=major_axis1[i];
    bnew->minor_axis_length1=minor_axis1[i];
    bnew->major_axis_length2=major_axis2[i];
    bnew->minor_axis_length2=minor_axis2[i];
    bnew->fluor=fluorescence[i];
    bnew->i_time=i_t;
    bnew->flag=flag;
    bnew->secs=secs;
    //Make area negative if we have a fret image and the cell we found
    //is a copy of another.
    bnew->a=area;
    bnew->n=n_points[i]; //is area in pixle units (without any requirements
    //that data be in the image or anything else). bnew->a is the actual
    //number of pixels used in the fluorescence sum calculation.
    if (have_fret_image==1){
      if(fret_copy_type[i]==1){
	bnew->a=(-(bnew->a));
	bnew->n=(-(bnew->n));
      }else if(fret_copy_type[i]<-1){
	printf("The %ith cell in the current image isn't marked",i);
	printf("whether it's an original or a copy.\n");
	fflush(stdout);
	perror(0); //exit(0); //http://r-pkgs.had.co.nz/src.html
      }
    }
    bnew->interior=interior[i];
    bnew->boundary=boundary[i];

		bnew->x_nucleus=x_nucl[i]; //V1.4.5
		bnew->y_nucleus=y_nucl[i];

    bnew->fft_stat=fft_stat[i];
    bnew->vacuole_area=vacuole_area[i];
    bnew->vacuole_fl=vacuole_fl[i];
    bnew->area_nucleus1=area_nucleus[i][0];
    bnew->area_nucleus2=area_nucleus[i][1];
    bnew->area_nucleus3=area_nucleus[i][2];
    bnew->area_nucleus4=area_nucleus[i][3];
    bnew->area_nucleus5=area_nucleus[i][4];
    bnew->area_nucleus6=area_nucleus[i][5];
//    bnew->area_nucleus7=area_nucleus[i][6];
//    bnew->area_nucleus8=area_nucleus[i][7];

		//V1.4 TODO
    bnew->fl_nucleus1=fl_nucleus[i][0];
    bnew->fl_nucleus2=fl_nucleus[i][1];
    bnew->fl_nucleus3=fl_nucleus[i][2];
    bnew->fl_nucleus4=fl_nucleus[i][3];
    bnew->fl_nucleus5=fl_nucleus[i][4];
    bnew->fl_nucleus6=fl_nucleus[i][5];
//    bnew->fl_nucleus7=fl_nucleus[i][6];
//    bnew->fl_nucleus8=fl_nucleus[i][7];
    bnew->fl_nucleus_from_search1=fl_nucleus_from_search[i][0];
    bnew->fl_nucleus_from_search2=fl_nucleus_from_search[i][1];
    bnew->fl_nucleus_from_search3=fl_nucleus_from_search[i][2];
    bnew->fl_nucleus_from_search4=fl_nucleus_from_search[i][3];
    bnew->fl_nucleus_from_search5=fl_nucleus_from_search[i][4];
    bnew->fl_nucleus_from_search6=fl_nucleus_from_search[i][5];
//    bnew->fl_nucleus_from_search7=fl_nucleus_from_search[i][6];
//    bnew->fl_nucleus_from_search8=fl_nucleus_from_search[i][7];

    bnew->pos_sig_mean_x=pos_sig_mean_x[i];
    bnew->pos_sig_mean_y=pos_sig_mean_y[i];
    bnew->pos_sig_mean_r=pos_sig_mean_r[i];
    bnew->pos_sig_rms_r=pos_sig_rms_r[i];
    bnew->neg_sig_mean_x=neg_sig_mean_x[i];
    bnew->neg_sig_mean_y=neg_sig_mean_y[i];
    bnew->neg_sig_mean_r=neg_sig_mean_r[i];
    bnew->neg_sig_rms_r=neg_sig_rms_r[i];
    //Stuff related to changing radius by different pixels
    bnew->area_p1=cell_area_p1[i];
    bnew->area_m1=cell_area_m1[i];
    bnew->area_m2=cell_area_m2[i];
    bnew->area_m3=cell_area_m3[i];
    bnew->fluor_p1=fluorescence_p1[i];
    bnew->fluor_m1=fluorescence_m1[i];
    bnew->fluor_m2=fluorescence_m2[i];
    bnew->fluor_m3=fluorescence_m3[i];

    //Some local background stuff
    bnew->back_p5=back_p5[i];
    bnew->pixels_back_p5=pixels_back_p5[i];
    bnew->pixels_total_p5=pixels_total_p5[i];
    bnew->back_phalfminor=back_phalfminor[i];
    bnew->pixels_back_phalfminor=pixels_back_phalfminor[i];
    bnew->pixels_total_phalfminor=pixels_total_phalfminor[i];

    //Find cell in known list that is closest.
    id_fret=9999;
    if(new_phase==1){ //This is a new list of cells

      iloc=i;

      if ((have_fret_image==1)&&(i>=n_before_fret_copy)){
	//For the FRET image type, only do the comparison for the
	//id<fret_offset cells. Since the n_found cells that i is looping
	//over all have the id+fret_offset cells in the upper half of the
	//list, we can assume that if i>=n_before_fret_copy, then we can use
	//the overlap from the cell=i-n_before_fret_copy cell. This cell is
	//labelled by i-n_before_fret_copy, and we can assume it's already
	//been calculated.
	//We have the id number that the cell (i-n_before_fret_copy) was
	//identified with through location_in_cs_array[]. We'll just use
	//that below, so do nothing here. (Note that this assumes that
	//we run over the non-fret_offset cells _first_ in the loop
	//i=1,n_found above.)
	j=location_in_cs_array[i-n_before_fret_copy];
	if (j<0){
	  free(bnew);
	  continue; //This cell was skipped
	}
	id_fret=cs[j]->index;
	id_fret+=fret_offset; //For the >fret_offset cells
	//Now find the cs[] element that has this id number.
	j_max=-1;
	I_over_U_max=-1.0; //To signal below for new cell
	for(j=(n_known-1);j>=0;j--){
	  if((cs[j]->index)==id_fret){
	    j_max=j; //Where we should put this cell
	    I_over_U_max=999.0; //Will always pass cut below for match
	    break;
	  }
	}

      }else{

	I_over_U_max=-1.0;
	//Prepare for calculation over intersection with all known points below
	update_overlap_value(); //Do each cell separately
	fill_overlap_array_with_point_list(interior[i]);
	n_p=n_points[i];

	/*
	  area=area_of_cell(boundary[i]);
	  circ=circumference_of_cell(boundary[i]);
	  if(area!=0.0){
	  c_over_a1=circ*circ/area;
	  }else{
	  c_over_a1=1.0e30;
	  }
	*/

	j_max=-1;
	tmp_j=(fret_offset+overall_id_offset);
	for(j=0;j<loop_total;j++){
	  //Look at previous time point.  Last time point is cs[j], however
	  //it's possible that we just added a point, so go to previous
	  //time.
	  if ((cs[j]->index)>=tmp_j){
	    continue; //No need to compare to id's above fret_offset
	  }
	  if (cs[j]->x<0) continue; //A flag for removing cells
	  for(b=cs[j];(b!=NULL)&&(b->i_time==i_t);b=b->prev);
	  if(b==NULL){ //Should never happen
	    printf("Chose a list that had no previous elements!!!\n");
	    printf("!!!!!!!!!! (%i, %i) !!!!!!!!!\n",i_t,j);
	    break;
	  }
	  isection=overlap(b->interior,offset_i,offset_j);
	  uon=((b->n)+n_p) - isection; //union of points
	  I_over_U=((float)isection)/((float)uon);
	  if(I_over_U>I_over_U_max){
	    I_over_U_max=I_over_U;
	    j_max=j;
	  }
	}
	if ((j_max<0)&&(loop_total>0)){
	  printf("Failed to find any points to compare to!\n");
	  I_over_U_max=-1.0;
	}
      } //if we have a fret-image and we're in the id+fret_offset range

    }else{ //if we don't have a new list of cells
      j_max=location_in_cs_array[i];
      if (j_max<0){
	free(bnew);
	continue; //Continue with next cell since we
	//decided to skip this cell when we looked at it before.
      }
    }

    bnew->i_over_u=I_over_U_max;
    //Check if we're close enough to the one we think matches
    if(I_over_U_max>I_over_U_for_match){
      add_new_cell=0;
    }else{
      add_new_cell=1;
    }
    if (add_new_cell==0){ //Not a new cell since matched to previous time
      b=cs[j_max]; //Add new cell just after this one in linked list
      //Check if we've already added a cell at this time point
      if(b->i_time==i_t){ //Already matched a cell here
	add_new_cell=1; //For check below
	//Should we use the current cell or the one we added before?
	if(I_over_U_max>(b->i_over_u)){ //Match current, make previous
	  //cell a new one (otherwise make current a new cell, which will
	  //happen automatically since add_new_cell is now set to 1)
	  bsave=b; //Cell we're replacing
	  (b->prev)->next=bnew;
	  bnew->prev=(b->prev); //b now removed from list
	  bnew->index=(b->prev)->index;
	  cs[j_max]=bnew; //cs[] always points to latest time point
	  //Find index that we're replacing
	  for(iloc=0;iloc<i;iloc++){
	    if (location_in_cs_array[iloc]==j_max) break;
	  }
	  //iloc is location_in_cs_array[] for cell we're moving.
	  location_in_cs_array[i]=j_max; //Save this location for cell i
	  bnew=bsave; //add_new_cell=1 so will add this one to the list
	}
      }
    }
    if (add_new_cell==0){ //Might have been changed above
      (b->next)=bnew;
      bnew->prev=b;
      bnew->index=b->index;
      cs[j_max]=bnew; //cs[] always points to latest time point
      location_in_cs_array[i]=j_max; //Save this location
    }else{

      if (n_known<max_cells){

	  //Add a new cell to the list
	location_in_cs_array[iloc]=n_known; //Save this location
	cs[n_known]=bnew;
	bnew->prev=NULL;
	//bnew->index=n_known+overall_id_offset;
	if ((have_fret_image==1)&&(i>=n_before_fret_copy)){
	  bnew->index=id_fret;
	}else{
	  if ((next_index<fret_offset)||(have_fret_image!=1)||
	      (i>=n_before_fret_copy)){
	    bnew->index=next_index+overall_id_offset;
	    next_index++; //Will be different from n_known for fret guys
	  }else{ //Skip this cell
	    printf("************** Too many known cells*****: %i.\n",n_known);
	    location_in_cs_array[iloc]=-999;
	    free(bnew);
	    continue;
	  }
	}
	n_known++;
      }else{ //Skip this cell
	printf("************** Too many known cells*****: %i.\n",n_known);
	location_in_cs_array[iloc]=-999;
	free(bnew);
	continue;
      }

    }

  } //End loop over found cells

  printf("Currently found total of %i separate cells.\n",n_known);

  return;
}
/**********************************************************/
void free_pixels_from_earlier_time_points(){
  //Go through all the known cells and free interior[] and boundary[]
  //points from the cs[]... array. Keep the latest points but free
  //all the previous. (Note that they might point to the same list if
  //there were multiple FL images per bright field.)
  int j;
  struct point *p;
  struct blob *b;

  for(j=0;j<n_known;j++){
    p=cs[j]->interior;
    for(b=cs[j];(b!=NULL)&&((b->interior)==p);b=b->prev);
    if (b!=NULL){
      point_list_free(b->interior);
      point_list_free(b->boundary);
      b->interior=NULL;
      b->boundary=NULL;
    }
  }
  return;
}
/**********************************************************/
int total_overlap_of_all_cells(offset_i,offset_j){

  int total;
  int i;

  total=0;
  for(i=0;i<n_found;i++){
    total += overlap(interior[i],offset_i,offset_j);
  }

  //printf("Current overlap: (%i,%i) gives %i.\n",offset_i,offset_j,total);

  return total;
}

/**********************************************************/
void add_cell_number_to_the_data(int i_t){
  //Loop over the list of known cells and add to the data a number label
  //based on the cells location in the known list.  Only do those cells
  //that were found during iteration i_t.

  int i;
  struct blob *b;

  int alternate;
  int count;

  if(image_type==hexagonal_grid){
    alternate=10;
  }else{
    alternate=1;
  }

  count=0;
  for(i=0;i<n_known;i++){
    b=cs[i]; //cs[] points to most recently added cell
    if(b->i_time==i_t){
      //Write number i at position (b->x,b->y)
      //     printf("Adding number %i at (%e,%e).\n",b->index,b->x,b->y);
      //fflush(stdout);

      if((count%alternate)==0){
	       if (((b->x)>=0.0)&&((b->y)>=0.0)){
	         add_numbers_to_data(b->index,
			       ((int)((b->x))),
			       ((int)((b->y))),
			       d,xmax,ymax);
	       }
      }
      count++;

    }
  }

  return;
}


/**********************************************************/
//V1.2a All output goes to a single file (no parts)
int output_cells_single_file(char *basename, char *append, int *time_index){
  //basename is where file goes+beginning of output name
  //append is "" for overwrite or "a" for append

  int i,j;
  struct blob *b;
  FILE *fp=NULL;
  char file[500];

  //char sys_cmd[500];
  int time1,time2;
  float back;

  char wa[3];

  if ((append[0]=='a')||(append[0]=='A')){
    strcpy(wa,"a"); //write and append
  }else{
    strcpy(wa,"w"); //write and overwrite
  }


  //First sort the cs[] list based on how many timepoints are in each
  //of the linked lists (ie, how many nodes in each linked list).
  for(i=0;i<(n_known-1);i++){
    b=cs[i];
    time1=0;
    while((b->prev)!=NULL){
      b=b->prev;
      time1++;
    }
    for(j=i+1;j<n_known;j++){
      b=cs[j];
      time2=0;
      while((b->prev)!=NULL){
	b=b->prev;
	time2++;
      }
      if(time2>time1){
	b=cs[i];
	cs[i]=cs[j];
	cs[j]=b;
	time1=time2;
      }
    }
  }

  //Put all cells into one file --4/24/03, so no more writing out a
  //separate file for each and every cell.
  strcpy(file,basename);
  strcat(file,"_all");

  if((fp=fopen(file,wa))==NULL){
    printf("Couldn't open single output file %s\n",file);
    fflush(stdout);
    return 0;
  }

  //Printing header to file
  fprintf(fp,"cellID \tt.frame\t time \t xpos \t ypos \t f.tot  \t a.tot  \t");
  fprintf(fp,"num.pix\tfft.stat \t perim  \t maj.axis \t min.axis \t f.nucl \t ");
  fprintf(fp,"a.nucl \t flag \t rot.vol \t con.vol \t a.vacuole \t f.vacuole \t f.bg   \t ");

  fprintf(fp,"f.tot.p1 \t a.tot.p1 \t f.tot.m1 \t a.tot.m1 \t f.tot.m2 \t");
  fprintf(fp,"a.tot.m2 \t f.tot.m3 \t a.tot.m3 \t");

	fprintf(fp,"xpos.nucl\typos.nucl\t ");
  fprintf(fp,"f.nucl1 \t f.nucl.tag1 \t a.nucl1 \t f.nucl2 \t f.nucl.tag2 \t ");
  fprintf(fp,"a.nucl2 \t f.nucl3 \t f.nucl.tag3 \t a.nucl3 \t f.nucl4 \t ");
  fprintf(fp,"f.nucl.tag4 \t a.nucl4 \t f.nucl5 \t f.nucl.tag5 \t a.nucl5 \t ");
  fprintf(fp,"f.nucl6 \t f.nucl.tag6 \t a.nucl6 \t ");

  fprintf(fp,"f.local.bg \t a.local.bg \t a.local \t f.local2.bg \t");
  fprintf(fp,"a.local2.bg \t a.local2 \t ");

  fprintf(fp,"a.surf \t con.vol \t sphere.vol ");

  fprintf(fp,"\n");


  for(i=0;i<n_known;i++){
    b=cs[i];
    time1=1;

    while(b->prev!=NULL){
      b=b->prev; //Go to first time point and count how many points
      time1++;
    }

    while(b!=NULL){
      fprintf(fp,"%4i\t%4i\t%4i\t%4i\t%4i\t%10.6e\t%10.6e\t%4i\t",
        b->index,
	      time_index[b->i_time],
	      b->secs,
	      ((int)((b->x)+0.5)),
	      ((int)((b->y)+0.5)),
	      b->fluor,
	      b->a,
	      b->n);
      fprintf(fp,"%10.6e\t",
	      b->fft_stat);
      fprintf(fp,"%10.6e\t",b->circumference);
      fprintf(fp,"%10.6e\t",b->major_axis_length1);
      fprintf(fp,"%10.6e\t",b->minor_axis_length1);

      fprintf(fp,"%10.6e\t",b->fl_nucleus3);
      fprintf(fp,"%10.6e\t",b->area_nucleus3);

      fprintf(fp,"%4i\t",b->flag);
      fprintf(fp,"%10.6e\t",b->vol_rotation);
      fprintf(fp,"%10.6e\t",b->vol_cone);

      fprintf(fp,"%10.6e\t",b->vacuole_area);
      fprintf(fp,"%10.6e\t",b->vacuole_fl);
      if (have_fret_image!=1){ //Don't have fret image
        back=back_pixels[b->i_time];
      }else{
        if ((b->index)<fret_offset){ //high y-cells (lower in image)
          back=back_pixels_1[b->i_time];
        }else{
          back=back_pixels_2[b->i_time];
        }
      }
      fprintf(fp,"%10.6e\t",back);
      //fprintf(fp,"\n");

      //And stuff related to changing the radius of the cells
      //(Put in new file so paw can read it all in)
      fprintf(fp,"%10.6e\t",b->fluor_p1);
      fprintf(fp,"%10.6e\t",b->area_p1);
      fprintf(fp,"%10.6e\t",b->fluor_m1);
      fprintf(fp,"%10.6e\t",b->area_m1);
      fprintf(fp,"%10.6e\t",b->fluor_m2);
      fprintf(fp,"%10.6e\t",b->area_m2);
      fprintf(fp,"%10.6e\t",b->fluor_m3);
      fprintf(fp,"%10.6e\t",b->area_m3);
      //fprintf(fp,"\n");


			//V1.4.5
			fprintf(fp,"%9i\t%9i\t",
	      b->x_nucleus,b->y_nucleus);
      //no se si tiene sentido usar estos valores
			fprintf(fp,"%10.6e\t%10.6e\t%10.6e\t",
	      b->fl_nucleus1,b->fl_nucleus_from_search1,b->area_nucleus1);
      fprintf(fp,"%10.6e\t%10.6e\t%10.6e\t",
	      b->fl_nucleus2,b->fl_nucleus_from_search2,b->area_nucleus2);
      fprintf(fp,"%10.6e\t%10.6e\t%10.6e\t",
	      b->fl_nucleus3,b->fl_nucleus_from_search3,b->area_nucleus3);
      fprintf(fp,"%10.6e\t%10.6e\t%10.6e\t",
	      b->fl_nucleus4,b->fl_nucleus_from_search4,b->area_nucleus4);
      fprintf(fp,"%10.6e\t%10.6e\t%10.6e\t",
	      b->fl_nucleus5,b->fl_nucleus_from_search5,b->area_nucleus5);
      fprintf(fp,"%10.6e\t%10.6e\t%10.6e\t",
	      b->fl_nucleus6,b->fl_nucleus_from_search6,b->area_nucleus6);
      //fprintf(fp,"\n");

      //Local background stuff
      //Repeat id, time, and flag in file
      //fprintf(fp5,"% 4i % 4i % 4i",
	    //  b->index,
	    //  time_index[b->i_time],
	    //  b->flag);
      fprintf(fp,"%10.6e\t%10.6e\t%10.6e\t",
	      b->back_p5,
	      b->pixels_back_p5,
	      b->pixels_total_p5);
      fprintf(fp,"%10.6e\t%10.6e\t%10.6e\t",
	      b->back_phalfminor,
	      b->pixels_back_phalfminor,
	      b->pixels_total_phalfminor);
      //fprintf(fp5,"\n");

      //Volume calculations and v_effective calculations
      //Repeat id, time, and flag in file
      //fprintf(fp6,"% 4i % 4i % 4i ",
	    //  b->index,
	    //  time_index[b->i_time],
	    //  b->flag);
      fprintf(fp,"%10.6e\t%10.6e\t%10.6e",
	      b->surface_area,
	      b->vol_cone, //V1.3 se puede sacar, esta repetido
	      b->vol_sphere);
      //fprintf(fp,"%10.6e\t%10.6e\t%10.6e\t%10.6e\t%10.6e\t%10.6e\t",
	    //  b->vol_eff_1,
	    //  b->vol_eff_2,
	    //  b->vol_eff_3,
	    //  b->vol_eff_4,
	    //  b->vol_eff_5,
	    //  b->vol_eff_6);

      fprintf(fp,"\n");


      b=b->next;
    }
  }
  if (fp!=NULL)fclose(fp);

  printf("Done writing output, file closed.");
  return 1;
}

/**********************************************************/
int output_cells(char *basename, char *append, int *time_index){

  //basename is where file goes+beginning of output name
  //append is "" for overwrite or "a" for append

  int i,j;
  struct blob *b;
  FILE *fp=NULL;
  FILE *fp2=NULL;
  FILE *fp3=NULL;
  FILE *fp4=NULL;
  FILE *fp5=NULL;
  FILE *fp6=NULL;
  char file[500];
  char file2[500];
  char file3[500];
  char file4[500];
  char file5[500];
  char file6[500];
  //char sys_cmd[500];
  int time1,time2;
  float back;

  char wa[3];

  if ((append[0]=='a')||(append[0]=='A')){
    strcpy(wa,"a"); //write and append
  }else{
    strcpy(wa,"w"); //write and overwrite
  }


  //First sort the cs[] list based on how many timepoints are in each
  //of the linked lists (ie, how many nodes in each linked list).
  for(i=0;i<(n_known-1);i++){
    b=cs[i];
    time1=0;
    while((b->prev)!=NULL){
      b=b->prev;
      time1++;
    }
    for(j=i+1;j<n_known;j++){
      b=cs[j];
      time2=0;
      while((b->prev)!=NULL){
	b=b->prev;
	time2++;
      }
      if(time2>time1){
	b=cs[i];
	cs[i]=cs[j];
	cs[j]=b;
	time1=time2;
      }
    }
  }

  //Put all cells into one file --4/24/03, so no more writing out a
  //separate file for each and every cell.
  strcpy(file,basename);
  strcat(file,"_all");

  strcpy(file2,file);
  strcat(file2,"_part2");

  strcpy(file3,file);
  strcat(file3,"_part3");

  strcpy(file4,file);
  strcat(file4,"_part4");

  strcpy(file5,file);
  strcat(file5,"_part5");

  strcpy(file6,file);
  strcat(file6,"_part6");

  if((fp=fopen(file,wa))==NULL){
    printf("Couldn't open fp file %s\n",file);
    fflush(stdout);
    return 0;
  }
  if((fp2=fopen(file2,wa))==NULL){
    printf("Couldn't open fp2 file %s\n",file2);
    fflush(stdout);
    return 0;
  }
  if((fp3=fopen(file3,wa))==NULL){
    printf("Couldn't open fp3 file %s\n",file3);
    fflush(stdout);
    return 0;
  }
  if((fp4=fopen(file4,wa))==NULL){
    printf("Couldn't open fp4 file %s\n",file5);
    fflush(stdout);
    return 0;
  }
  if((fp5=fopen(file5,wa))==NULL){
    printf("Couldn't open fp5 file %s\n",file5);
    fflush(stdout);
    return 0;
  }
  if((fp6=fopen(file6,wa))==NULL){
    printf("Couldn't open fp6 file %s\n",file5);
    fflush(stdout);
    return 0;
  }

  for(i=0;i<n_known;i++){
    b=cs[i];
    time1=1;

    while(b->prev!=NULL){
      b=b->prev; //Go to first time point and count how many points
      time1++;
    }

    /*
    //The names have two numbers.  First is number of nodes in the list
    //and the second is the number we've labelled this cell list with
    strcpy(file,basename);
    strcat(file,"_");
    digits_to_string(file,time1,total_time);
    strcat(file,"_");
    digits_to_string(file,b->index,n_known-1);

    strcpy(file2,file);
    strcat(file2,"_part2");

    if((fp=fopen(file,"w"))==NULL){
      printf("Couldn't open file %s\n",file);
      fflush(stdout);
      return 0;
    }
    if((fp2=fopen(file2,"w"))==NULL){
      printf("Couldn't open file %s\n",file2);
      fflush(stdout);
      return 0;
    }
    */

    while(b!=NULL){
      fprintf(fp,"% 4i % 4i % 4i % 4i  % 4i  % 10.6e  % 10.6e  % 4i ",
	      b->index,
	      time_index[b->i_time],
	      b->secs,
	      ((int)((b->x)+0.5)),
	      ((int)((b->y)+0.5)),
	      b->fluor,
	      b->a,
	      b->n);
      fprintf(fp,"% 10.6e ",
	      b->fft_stat);
      //asg--7/31/03--Don't really use these statistics any more.
      //Start replacing them with circumference, etc.
      //fprintf(fp,"% 10.6e ",b->pos_sig_mean_x);
      //fprintf(fp,"% 10.6e ",b->pos_sig_mean_y);
      //fprintf(fp,"% 10.6e ",b->pos_sig_mean_r);
      //fprintf(fp,"% 10.6e ",b->pos_sig_rms_r);
      //fprintf(fp,"% 10.6e ",b->neg_sig_mean_x);
      fprintf(fp,"% 10.6e ",b->circumference);
      fprintf(fp,"% 10.6e ",b->major_axis_length1);
      fprintf(fp,"% 10.6e ",b->minor_axis_length1);

      //fprintf(fp,"% 10.6e ",b->major_axis_length2);
      //fprintf(fp,"% 10.6e ",b->minor_axis_length2);
      //asg--9/14/04-->Fill this position with the nuclear stuff.
      //Changed in output_...9.kumac's also.
      fprintf(fp,"% 10.6e ",b->fl_nucleus3);
      fprintf(fp,"% 10.6e ",b->area_nucleus3);

      //      fprintf(fp,"% 10.6e ",b->neg_sig_mean_y);
      fprintf(fp,"% 4i ",b->flag);
      fprintf(fp,"% 10.6e ",b->vol_rotation);
      //fprintf(fp,"% 10.6e ",b->neg_sig_mean_r);
      //fprintf(fp,"% 10.6e ",b->neg_sig_rms_r);
      fprintf(fp,"% 10.6e ",b->vol_cone);

      //fprintf(fp,"% 10.6e ",total_area[b->i_time]);
      //fprintf(fp,"% 10.6e ",total_fluorescence[b->i_time]);
      fprintf(fp,"% 10.6e ",b->vacuole_area);
      fprintf(fp,"% 10.6e ",b->vacuole_fl);
      if (have_fret_image!=1){ //Don't have fret image
				back=back_pixels[b->i_time];
      }else{
				if ((b->index)<fret_offset){ //high y-cells (lower in image)
	  			back=back_pixels_1[b->i_time];
				}else{
	  			back=back_pixels_2[b->i_time];
				}
      }
      fprintf(fp,"% 10.6e ",back);
      fprintf(fp,"\n");
      //And stuff related to changing the radius of the cells
      //(Put in new file so paw can read it all in)
      fprintf(fp2,"% 10.6e ",b->fluor_p1);
      fprintf(fp2,"% 10.6e ",b->area_p1);
      fprintf(fp2,"% 10.6e ",b->fluor_m1);
      fprintf(fp2,"% 10.6e ",b->area_m1);
      fprintf(fp2,"% 10.6e ",b->fluor_m2);
      fprintf(fp2,"% 10.6e ",b->area_m2);
      fprintf(fp2,"% 10.6e ",b->fluor_m3);
      fprintf(fp2,"% 10.6e ",b->area_m3);
      fprintf(fp2,"\n");

      fprintf(fp3,"% 10.6e % 10.6e % 10.6e ",
	      b->fl_nucleus1,b->fl_nucleus_from_search1,b->area_nucleus1);
      fprintf(fp3,"% 10.6e % 10.6e % 10.6e ",
	      b->fl_nucleus2,b->fl_nucleus_from_search2,b->area_nucleus2);
      fprintf(fp3,"% 10.6e % 10.6e % 10.6e ",
	      b->fl_nucleus3,b->fl_nucleus_from_search3,b->area_nucleus3);
      fprintf(fp3,"% 10.6e % 10.6e % 10.6e ",
	      b->fl_nucleus4,b->fl_nucleus_from_search4,b->area_nucleus4);
      fprintf(fp3,"% 10.6e % 10.6e % 10.6e ",
	      b->fl_nucleus5,b->fl_nucleus_from_search5,b->area_nucleus5);
      fprintf(fp3,"% 10.6e % 10.6e % 10.6e ",
	      b->fl_nucleus6,b->fl_nucleus_from_search6,b->area_nucleus6);
      fprintf(fp3,"% 10.6e % 10.6e % 10.6e ",
	      0.0,0.0,0.0);
      fprintf(fp3,"% 10.6e % 10.6e % 10.6e ",
	      0.0,0.0,0.0);
      fprintf(fp3,"\n");

			//V1.4.5
      fprintf(fp4,"% 4i % 4i ",
	      b->x_nucleus,b->y_nucleus);
      fprintf(fp4,"\n");


      //Local background stuff
      //Repeat id, time, and flag in file
      fprintf(fp5,"% 4i % 4i % 4i",
	      b->index,
	      time_index[b->i_time],
	      b->flag);
      fprintf(fp5,"% 10.6e % 10.6e % 10.6e ",
	      b->back_p5,
	      b->pixels_back_p5,
	      b->pixels_total_p5);
      fprintf(fp5,"% 10.6e % 10.6e % 10.6e ",
	      b->back_phalfminor,
	      b->pixels_back_phalfminor,
	      b->pixels_total_phalfminor);
      fprintf(fp5,"\n");

      //Volume calculations and v_effective calculations
      //Repeat id, time, and flag in file
      fprintf(fp6,"% 4i % 4i % 4i ",
	      b->index,
	      time_index[b->i_time],
	      b->flag);
      fprintf(fp6,"% 10.6e % 10.6e % 10.6e ",
	      b->surface_area,
	      b->vol_cone,
	      b->vol_sphere);
      fprintf(fp6,"% 10.6e % 10.6e % 10.6e % 10.6e % 10.6e % 10.6e ",
	      b->vol_eff_1,
	      b->vol_eff_2,
	      b->vol_eff_3,
	      b->vol_eff_4,
	      b->vol_eff_5,
	      b->vol_eff_6);
      fprintf(fp6,"\n");


      b=b->next;
    }
  }
  if (fp!=NULL)fclose(fp);
  if (fp2!=NULL)fclose(fp2);
  if (fp3!=NULL)fclose(fp3);
  if (fp4!=NULL)fclose(fp4);
  if (fp5!=NULL)fclose(fp5);
  if (fp6!=NULL)fclose(fp6);

  return 1;
}

/**********************************************************/
void all_free(){
  //Free all the memory that we've allocated

  pmem_cur=pmem_start;
  pmem_cur->next=NULL;

  while((--next_block)>=0) free(mem_blocks[next_block]);

  next_block=0; //Just to be sure
  return;
}

/**********************************************************/
int calculate_fluorescence_with_r_info(void){

  int i,j;
  struct point *p;
  double total;
  float count;
  int u;
  int ix,iy;

  //Stuff for local background calculation
  struct point *p1;
  struct point *p2;
  float dx,dy,x0,y0,x,y;
  int ix2,iy2;
  float step=.01;
  float s;
  float *flback[2];
  float *pixelsback[2];
  float *pixelstotal[2];
  struct point **bstart[2];
  struct point *b;

  int fret_region;

  int *clist_x;
  int *clist_y;
  int n_vacuole;

  //double dtmp,dtmp1,dtmp2;
  //double mu1,mu2;
  //int j;
  static int asg_time_tmp=-1;
  //static FILE *fp_tmp=NULL;

  //double mu,rms,val,tmp;
  //int j,k,iy2,u2;

  asg_time_tmp++;

  //Loop over all the interior points (assume they're already calculated)
  //and sum up fl[] array (assume it's already filled with fluorescence
  //data.
  //All calculate the distance of each pixel from the boundary.

  //If we have a third image that is a vacuole-label, then we
  //want to exclude the vacuole part of the image. To do this, we
  //have set the d[] array in the internal_structure. To speed things
  //up I do the check up front, and then just write the code twice (so
  //in case we don't want to do the check, we don't have to figure that
  //out at every point).
  //  if (third_image_type==vacuole_label){
  if (0==1){ //Asg--never do this part

    for(i=0;i<n_found;i++){
      total=0.0;
      count=0.0;
      for(p=interior[i];p!=NULL;p=p->next){
	ix=(p->i);
	iy=(p->j);
	if ((ix>=0)&&(ix<xmax)&&(iy>=0)&&(iy<ymax)){
	  u=iy*xmax+ix;
	  if (d[u]!=cell_nucleus){
	    total+=(double)(fl[u]);
	    count+=1.0;
	  }
	}
      }
      fluorescence[i]=((float)total);
      cell_area[i]=count; //See comment above in calculate_fluorecence()

      //radius increased by 1
      total=0.0;
      count=0.0;
      for(p=interior_p1[i];p!=NULL;p=p->next){
	ix=(p->i);
	iy=(p->j);
	if ((ix>=0)&&(ix<xmax)&&(iy>=0)&&(iy<ymax)){
	  u=iy*xmax+ix;
	  if (d[u]!=cell_nucleus){
	    total+=(double)(fl[u]);
	    count+=1.0;
	  }
	}
      }
      fluorescence_p1[i]=((float)total);
      cell_area_p1[i]=count;

      //radius decreased by 1
      total=0.0;
      count=0.0;
      for(p=interior_m1[i];p!=NULL;p=p->next){
	ix=(p->i);
	iy=(p->j);
	if ((ix>=0)&&(ix<xmax)&&(iy>=0)&&(iy<ymax)){
	  u=iy*xmax+ix;
	  if (d[u]!=cell_nucleus){
	    total+=(double)(fl[u]);
	    count+=1.0;
	  }
	}
      }
      fluorescence_m1[i]=((float)total);
      cell_area_m1[i]=count;

      //radius decreased by 2
      total=0.0;
      count=0.0;
      for(p=interior_m2[i];p!=NULL;p=p->next){
	ix=(p->i);
	iy=(p->j);
	if ((ix>=0)&&(ix<xmax)&&(iy>=0)&&(iy<ymax)){
	  u=iy*xmax+ix;
	  if (d[u]!=cell_nucleus){
	    total+=(double)(fl[u]);
	    count+=1.0;
	  }
	}
      }
      fluorescence_m2[i]=((float)total);
      cell_area_m2[i]=count;

      //radius decreased by 3
      total=0.0;
      count=0.0;
      for(p=interior_m3[i];p!=NULL;p=p->next){
	ix=(p->i);
	iy=(p->j);
	if ((ix>=0)&&(ix<xmax)&&(iy>=0)&&(iy<ymax)){
	  u=iy*xmax+ix;
	  if (d[u]!=cell_nucleus){
	    total+=(double)(fl[u]);
	    count+=1.0;
	  }
	}
      }
      fluorescence_m3[i]=((float)total);
      cell_area_m3[i]=count;

    }

  }else{

    for(i=0;i<n_found;i++){

      total=0.0;
      count=0.0;
      for(p=interior[i];p!=NULL;p=p->next){
	ix=(p->i);
	iy=(p->j);
	if ((ix>=0)&&(ix<xmax)&&(iy>=0)&&(iy<ymax)){
	  u=iy*xmax+ix;
	  total+=(double)(fl[u]);
	  count+=1.0;
	}
      }
      fluorescence[i]=((float)total);
      cell_area[i]=count; //See comment above in calculate_fluorecence()

      //radius increased by 1
      total=0.0;
      count=0.0;
      for(p=interior_p1[i];p!=NULL;p=p->next){
	ix=(p->i);
	iy=(p->j);
	if ((ix>=0)&&(ix<xmax)&&(iy>=0)&&(iy<ymax)){
	  u=iy*xmax+ix;
	  total+=(double)(fl[u]);
	  count+=1.0;
	}
      }
      fluorescence_p1[i]=((float)total);
      cell_area_p1[i]=count;

      //radius decreased by 1
      total=0.0;
      count=0.0;
      for(p=interior_m1[i];p!=NULL;p=p->next){
	ix=(p->i);
	iy=(p->j);
	if ((ix>=0)&&(ix<xmax)&&(iy>=0)&&(iy<ymax)){
	  u=iy*xmax+ix;
	  total+=(double)(fl[u]);
	  count+=1.0;
	}
      }
      fluorescence_m1[i]=((float)total);
      cell_area_m1[i]=count;

      //radius decreased by 2
      total=0.0;
      count=0.0;
      for(p=interior_m2[i];p!=NULL;p=p->next){
	ix=(p->i);
	iy=(p->j);
	if ((ix>=0)&&(ix<xmax)&&(iy>=0)&&(iy<ymax)){
	  u=iy*xmax+ix;
	  total+=(double)(fl[u]);
	  count+=1.0;
	}
      }
      fluorescence_m2[i]=((float)total);
      cell_area_m2[i]=count;

      //radius decreased by 3
      total=0.0;
      count=0.0;
      for(p=interior_m3[i];p!=NULL;p=p->next){
	ix=(p->i);
	iy=(p->j);
	if ((ix>=0)&&(ix<xmax)&&(iy>=0)&&(iy<ymax)){
	  u=iy*xmax+ix;
	  total+=(double)(fl[u]);
	  count+=1.0;
	}
      }
      fluorescence_m3[i]=((float)total);
      cell_area_m3[i]=count;

      find_vacuole(interior_m3[i],fl,xmax,ymax,
		   &clist_x,
		   &clist_y,
		   &n_vacuole);
      //Sum of fluorescence associated with my vacuole location
      total=0.0;
      count=0.0;
      for (j=0;j<n_vacuole;j++){
	ix=clist_x[j];
	iy=clist_y[j];
	if ((ix>=0)&&(ix<xmax)&&(iy>=0)&&(iy<ymax)){
	  u=iy*xmax+ix;
	  total+=(double)(fl[u]);
	  count+=1.0;
	}
      }

      vacuole_area[i]=count;
      vacuole_fl[i]=((float)total);

    }

  }

  //Now do local background calculations using
  //boundary_p5 and boundary_phalfminor. To do this, mark the
  //pixels that correspond to found cells. Use the +1 radii
  update_overlap_value();
  for(i=0;i<n_found;i++){
    fill_overlap_array_with_point_list(interior_p1[i]);
  }
  //Copy stuff from "add_boundary_points_to_data" code. The boundary points
  //aren't necessarily contiguous....We could do interior_p5 and then
  //take the difference from interior_p4, but this should do the same.
  //Put the pointers into an array so can use the same code.

  flback[0]=back_p5;
  pixelsback[0]=pixels_back_p5;
  pixelstotal[0]=pixels_total_p5;
  bstart[0]=boundary_p5;

  flback[1]=back_phalfminor;
  pixelsback[1]=pixels_back_phalfminor;
  pixelstotal[1]=pixels_total_phalfminor;
  bstart[1]=boundary_phalfminor;

  for (i=0;i<n_found;i++){
    //For fret image, don't cross out of fret-region. Use the
    //the fret region for this cell.
    if (have_fret_image==1){
      p=interior[i];
      ix=(p->i);
      iy=(p->j);
      if ((ix>=0)&&(ix<xmax)&&(iy>=0)&&(iy<ymax)){
	u=iy*xmax+ix;
	fret_region=fret_labels[u];
      }
    }
    for(j=0;j<2;j++){
      (flback[j])[i]=0.0;
      (pixelsback[j])[i]=0.0;
      (pixelstotal[j])[i]=0.0;
      b=(bstart[j])[i];
      for(p1=b;p1!=NULL;p1=p1->next){
	p2=p1->next;
	if(p2==NULL) p2=b; //So make full circuit

	dx=(float)( (p2->i)-(p1->i) );
	dy=(float)( (p2->j)-(p1->j) );
	ix=p1->i;
	iy=p1->j;
	x0=(float) ix;
	y0=(float) iy;
	for(s=step;s<=1.0;s+=step){
	  x=x0+dx*s+0.5; //0.5 to round off (place integers in middle of bins)
	  y=y0+dy*s+0.5;
	  ix2=(int) x;
	  iy2=(int) y;
	  if(((ix!=ix2)||(iy!=iy2))&&
	     (ix2>=0)&&(ix2<xmax)&&(iy2>=0)&&(iy2<ymax)){
	    //New point
	    ix=ix2;
	    iy=iy2;
	    if ((ix>=0)&&(ix<xmax)&&(iy>=0)&&(iy<ymax)){
	      u=iy*xmax+ix;
	      (pixelstotal[j])[i]+=1.0;
	      if (over[u]==overlap_value) continue;
	      if (have_fret_image==1){//If we have a fret image
		if (fret_labels[u]!=fret_region) continue;
	      }
	      //If we're here, then we're not in another cell, or
	      //we haven't crossed a fret-boundary for fret images.
	      (pixelsback[j])[i]+=1.0;
	      (flback[j])[i]+=(fl[u]);
	    }
	  }
	}
      }
      //Make into an average
      if (((pixelsback[j])[i])>0.0){
	((flback[j])[i])/=((pixelsback[j])[i]);
      }else{
	(flback[j][i])=0.0;
      }
    }

  }

return 1;
}

/**********************************************************/
int get_data_from_text_file(char *file){
  //f is location to write the data

  FILE *fp_in;
  int i,j;
  int k;
  float tmp;

  if( (fp_in=fopen(file,"r"))==NULL ){
    printf("Couldnt open file %s (get_data_from_text_file).\n",file);
    return 0;
  }

  //Zero out the data to start since file might not contain every point.
  for(i=0;i<xmax_ymax;i++){
      c[i]=0.0;
  }

  //Read the data in form (i_index, j_index, data).
  //  printf("Reading data...\n");
  while(fscanf(fp_in, "%i %i %i", &i, &j, &k )==3){
    //The "l" means its a double rather than a float
    //3=number of items read
    if((i<0)||(i>=xmax)||(j<0)||(j>=ymax)){
      printf("Point out of bounds: %i %i %i.\n",i,j,k);
      return 0;
    }
    tmp=(float) k;
    c[(j*xmax+i)]=tmp;
    //    d[i][j]=tmp; //Use d[][] to find starting points.
  }

  fclose(fp_in);
  //  printf("done reading.\n");

return 1;
}

/**********************************************************/
void debug_test(int flag){

  printf("Here_debug: %i\n",flag);fflush(stdout);
  printf("--->%i",(int)interior_m1[0]);
  if (interior_m1[0]!=NULL){
    printf("----------------: (%i,%i)\n",interior_m1[0]->i,
	   interior_m1[0]->j);
  }else{
    printf("\n");
  }
}

/**********************************************************/
void align_found_cells_to_fl(int flag){

  int i;
  int itmpx,itmpy;
  int offx[max_cells];
  int offy[max_cells];
  int overlap_value_save;
  int change;
  int max_trials=10;
  int i_trials;

  printf("In align\n");fflush(stdout);
  if (flag==1){ //Just re-set the originals and return
    for(i=0;i<n_found;i++){
      itmpx=(-save_realign_offset_x[i]); //Re-set
      itmpy=(-save_realign_offset_y[i]);
      point_list_adjust(interior[i],itmpx,itmpy);
      point_list_adjust(interior_m1[i],itmpx,itmpy);
      point_list_adjust(interior_m2[i],itmpx,itmpy);
      point_list_adjust(interior_m3[i],itmpx,itmpy);
      point_list_adjust(interior_p1[i],itmpx,itmpy);
      point_list_adjust(boundary[i],itmpx,itmpy);
      point_list_adjust(boundary_m1[i],itmpx,itmpy);
      point_list_adjust(boundary_m2[i],itmpx,itmpy);
      point_list_adjust(boundary_m3[i],itmpx,itmpy);
      point_list_adjust(boundary_p1[i],itmpx,itmpy);
    }

    return;
  }

  //For each cell, calculate new set of interior points

  //Initialize offsets to 0
  for(i=0;i<n_found;i++){
    offx[i]=0;
    offy[i]=0;
    save_realign_offset_x[i]=0;
    save_realign_offset_y[i]=0;
  }

  i_trials=0;
  do{
    i_trials++;
    if (i_trials>max_trials){
      printf("Too many trials for individual-cell re-alignment:%i>%i\n",
	     i_trials,max_trials);
      break;
    }

    change=0;//Change will be set if any offset has changed

    //Set overlap array with current positions
    update_overlap_value();
    overlap_value_save=overlap_value;
    for(i=0;i<n_found;i++){
      fill_overlap_array_with_point_list_offset(interior[i],
						offx[i],
						offy[i]);
    }
    for(i=0;i<n_found;i++){
      //Remove this one from the overlap list
      overlap_value=0;
      fill_overlap_array_with_point_list_offset(interior[i],
						offx[i],
						offy[i]);
      overlap_value=overlap_value_save;
      itmpx=offx[i];
      itmpy=offy[i];
      if (flag==0) {
	       re_align_cell(interior[i],(offx+i),(offy+i));
      }else{
	       re_align_cell(boundary[i],(offx+i),(offy+i));
      }
      if ((offx[i]!=itmpx)||(offy[i]!=itmpy)) change=1;
      //Add it back
      fill_overlap_array_with_point_list_offset(interior[i],
						offx[i],
						offy[i]);
    }
  } while(change>0);

  //Done with re-alignment
  //Save the offsets for returning to original, and
  //adjust all lists with the offset
  for(i=0;i<n_found;i++){
    save_realign_offset_x[i]=offx[i];
    save_realign_offset_y[i]=offy[i];
    itmpx=offx[i];
    itmpy=offy[i];
    point_list_adjust(interior[i],itmpx,itmpy);
    point_list_adjust(interior_m1[i],itmpx,itmpy);
    point_list_adjust(interior_m2[i],itmpx,itmpy);
    point_list_adjust(interior_m3[i],itmpx,itmpy);
    point_list_adjust(interior_p1[i],itmpx,itmpy);
    point_list_adjust(boundary[i],itmpx,itmpy);
    point_list_adjust(boundary_m1[i],itmpx,itmpy);
    point_list_adjust(boundary_m2[i],itmpx,itmpy);
    point_list_adjust(boundary_m3[i],itmpx,itmpy);
    point_list_adjust(boundary_p1[i],itmpx,itmpy);


  }

  return;
}

/****************************************************/
void point_list_adjust(struct point *p_in,int offx,int offy){

  struct point *ptmp;

  if (p_in==NULL) return;

  for (ptmp=p_in;ptmp!=NULL;ptmp=ptmp->next){
    ptmp->i=(ptmp->i)+offx;
    ptmp->j=(ptmp->j)+offy;
  }
  return;
}

/****************************************************/
void re_align_cell(struct point *p_in,int *offx_out, int *offy_out){

  //Just look left and right, minimum one pixels steps
  int offx,offy;

  struct point *p;

  int u,ix,iy;
  int offxmax,offymax;
  int offxstart,offystart;
  double totmax;
  double tot;

  //(*offx_out)=0;
  //(*offy_out)=0;

  if (p_in==NULL){
    return;
  }

  totmax=-1.0e30;
  offxmax=10;
  offymax=10;

  offxstart=(*offx_out);
  offystart=(*offy_out);

  offy=offystart;
  offx=offxstart-1;

 top:

  //First to the right
  for(offx=(offxstart+1);offx<=offxmax;offx++){
    tot=0.0;
    for(p=p_in;p!=NULL;p=p->next){
      ix=(p->i)+offx;
      iy=(p->j)+offy;
      if ((ix>=0)&&(ix<xmax)&&(iy>=0)&&(iy<ymax)){
	       u=(iy*xmax+ix);
	       if (over[u]!=overlap_value){
	         tot+=(fl[u]);
	       }
      }
    }
    if (tot>totmax){
      totmax=tot;
    }else{
      break;
    }
  }
  offx--;
  if (offx==offxstart){ //Check to the left
    for(offx=(offxstart-1);offx>=(-offxmax);offx--){
      tot=0.0;
      for(p=p_in;p!=NULL;p=p->next){
	     ix=(p->i)+offx;
	     iy=(p->j)+offy;
	       if ((ix>=0)&&(ix<xmax)&&(iy>=0)&&(iy<ymax)){
	         u=(iy*xmax+ix);
	         if (over[u]!=overlap_value){
	           tot+=(fl[u]);
	         }
	       }
      }
      if (tot>totmax){
	     totmax=tot;
      }else{
	     break;
      }
    }
    offx++;
  }

  //Now check up
  for(offy=(offystart+1);offy<=offymax;offy++){
    tot=0.0;
    for(p=p_in;p!=NULL;p=p->next){
      ix=(p->i)+offx;
      iy=(p->j)+offy;
      if ((ix>=0)&&(ix<xmax)&&(iy>=0)&&(iy<ymax)){
	u=(iy*xmax+ix);
	if (over[u]!=overlap_value){
	  tot+=(fl[u]);
	}
      }
    }
    if (tot>totmax){
      totmax=tot;
    }else{
      break;
    }
  }
  offy--;
  if (offy==offystart){ //Check to the left
    for(offy=(offystart-1);offy>=(-offymax);offy--){
      tot=0.0;
      for(p=p_in;p!=NULL;p=p->next){
	ix=(p->i)+offx;
	iy=(p->j)+offy;
	if ((ix>=0)&&(ix<xmax)&&(iy>=0)&&(iy<ymax)){
	  u=(iy*xmax+ix);
	  if (over[u]!=overlap_value){
	    tot+=(fl[u]);
	  }
	}
      }
      if (tot>totmax){
	totmax=tot;
      }else{
	break;
      }
    }
    offy++;
  }

  //printf("(%i,%i)-->(%i,%i)\n",offxstart,offystart,offx,offy);fflush(stdout);


  if ((offx!=offxstart)||(offy!=offystart)){
    offxstart=offx;
    offystart=offy;
    goto top;
  }

  (*offx_out)=offxstart;
  (*offy_out)=offystart;

  return;
}
/****************************************************/
void add_boundary_points_to_data(struct point *p_in, int blank_out_bg){

  struct point *p1;
  struct point *p2;
  int i;
  float step=.01;
  float mx,my,s;
  float x0,y0,x,y;
  int a,b,a2,b2;

  int border;
  border=found_border;  // tif_routines.h says: #define found_border 5

  struct point *p_start;

  //Add boundary points for border list p_in.
  //if p_in==NULL then do all n_found borders.
  //add found_border to d[] array in appropriate place.

  p_start=p_in;
  //for(i=0;i<xmax_ymax;i++)d[i]=0;
  for(i=0;i<n_found;i++){
    //calculate_volume_cone_old(interior[i]);
    if (p_in==NULL) p_start=boundary[i];  //if p_in==NULL then do all n_found borders.
    
    for(p1=p_start;p1!=NULL;p1=p1->next){
      p2=p1->next;
      if(p2==NULL) p2=p_start; //So make full circuit
      mx=(float)( (p2->i)-(p1->i) );
      my=(float)( (p2->j)-(p1->j) );
      a=p1->i;
      b=p1->j;
      x0=(float) a;
      y0=(float) b;
      for(s=step;s<=1.0;s+=step){
        x=x0+mx*s+0.5; //0.5 to round off (place integers in middle of bins)
        y=y0+my*s+0.5;
	      a2=(int) x;
        b2=(int) y;
        if(((a!=a2)||(b!=b2))&&(a2>=0)&&(a2<xmax)&&(b2>=0)&&(b2<ymax)){
          //New point
          a=a2;
          b=b2;
          // The following data in "d[(b*xmax+a)]=border;" ends up in labels[u] at tif.c
          // Parameter border=found_border; defined above.
          //   And tif_routines.h defines found_border as: #define found_border 5
          // If the blank_out_bg flag is set to 1, label the cells with their "number" in the for loop. This may not be the CellID.
          if(blank_out_bg==1){
            d[(b*xmax+a)]=i;
          } else {
            d[(b*xmax+a)]=border;
          }
	      }
      }
    }
    //A hack--break out of loop if only doing the passed in list.
    if(p_in!=NULL) break;
  }

  return;
}

/**********************************************************/
struct point *circularize_points(struct point *begin, float cut){
  //We start with the list pointed at by begin.  We calculate the value
  //of r=circumference^2/area for removing the first point and the last
  //point.  If r of either is below the value of cut, then we go ahead
  //and remove one of the points and repeat the process.  This should
  //remove weird cells that have big spikes at the end from connecting
  //the first and last points into a closed curve.
  //Note that (I believe) r is minimized for a circular shape (is 4*pi in
  //that case).  So r should always be above 4*pi, except for variations
  //from round-off errors and things being integerized, etc.


  struct point *pnew;
  struct point *p;
  struct point *psave;
  float r;
  float area;
  float circ;

  float rfirst,rlast;

  if(begin==NULL){
    return NULL;
  }

  //First try removing first point
  psave=begin;
  p=psave->next;
  if(p!=NULL){
    p->prev=NULL;
    area=area_of_cell(p);
    circ=circumference_of_cell(p);
    if(area!=0.0){
      r=circ*circ/area;
    }else{
      r=1.0e30;
    }
  }else{
    r=1.0e30;
  }
  rfirst=r;

  //Restore first point
  p->prev=psave;

  //Now try removing last point
  psave=begin;
  while(psave->next!=NULL)psave=psave->next;
  p=psave->prev; //p is now penultimate point, psave is last
  if(p!=NULL){
    p->next=NULL;
    area=area_of_cell(begin);
    circ=circumference_of_cell(begin);
    if(area!=0.0){
      r=circ*circ/area;
    }else{
      r=1.0e30;
    }
  }else{
    r=1.0e30;
  }
  rlast=r;
  //Restore point (Note p hasn't been changed from above)
  p->next=psave;

  //  printf("   Current=%e, first=%e, last=%e.\n",cut,rfirst,rlast);
  if((rfirst<cut)||(rlast<cut)){
    //Go ahead and remove first or last point
    if(rlast<rfirst){
      p->next=NULL; //Remove last (p is already pointing at penultimate point)
      pnew=begin; //pnew is now start of the new list
      r=rlast;
    }else{
      pnew=(begin->next);
      pnew->prev=NULL; //Remove first point
      r=rfirst;
    }
    return circularize_points(pnew,r); //Recursive iteration
  }else{ //Nothing to be gained by doing anything more
    return begin;
  }

}









/**********************************************************/
float circumference_of_cell(struct point *begin){
  //Calculate circumference for linked list of points pointed at by
  //begin.

  int j,k,u;
  float circ;
  struct point *p0;
  struct point *p1;

  //the circumference
  circ=0.0;
  for(p0=begin;p0!=NULL;p0=p0->next){
    p1=p0->next;
    if(p1==NULL)p1=begin; //So make full circuit
    j=(p1->i)-(p0->i);
    k=(p1->j)-(p0->j);
    u=j*j+k*k;
    circ+=( (float) sqrt( (double) u ) );
  }

  return circ;
}


/*****************************************************/
float max_d_over_s(struct point *begin,
		   struct point **max0,
		   struct point **max1){

  //For each point loop over all the other points and find the maximum
  //distanc over circumference travelled to get to that point

  struct point *p0;
  struct point *p1;
  struct point *ptmp;

  float ctot,cc,dd,tmp;
  float r,r_max;
  float r_global_max;
  int i,j,di,dj;
  int icur,jcur;
  int iprev,jprev;

  (*max0)=NULL;
  (*max1)=NULL;

  //Calculate total circumference
  ctot=0.0;
  p0=begin;
  while((p0->next)!=NULL)p0=p0->next; //p0 is now last in list
  iprev=p0->i;
  jprev=p0->j;
  for(p0=begin;p0!=NULL;p0=p0->next){
    di=(p0->i)-iprev;
    dj=(p0->j)-jprev;
    ctot+=((float)sqrt((double)(di*di+dj*dj)));
    iprev=p0->i;
    jprev=p0->j;
  }

  //printf("circ: %e out of %e.\n",ctot,smallest_circumference);
  if(ctot<smallest_circumference) return -1.0;



  r_global_max=-1.0;
  for(p0=begin;(p0->next)!=NULL;p0=p0->next){
    cc=0.0;
    r=0.0;
    r_max=-1.0;
    i=p0->i;
    j=p0->j;
    iprev=i;
    jprev=j;
    for(p1=p0->next;p1!=NULL;p1=p1->next){
      icur=(p1->i);
      jcur=(p1->j);

      di=(icur-i);
      dj=(jcur-j);
      dd=(float)(di*di+dj*dj);

      di=(icur-iprev);
      dj=(jcur-jprev);
      cc+=((float)sqrt((double)(di*di+dj*dj)));

      iprev=icur;
      jprev=jcur;

      //It might have been shorter to go the other way around
      if(cc>(ctot-cc)){
	      tmp=ctot-cc;
      }else{
	      tmp=cc;
      }
      if(tmp>5.0){
	      r=tmp*tmp/dd;  //Squared of (circumference/distance)

	      if(r>r_max){
	        r_max=r;
	        ptmp=p1; //Save location
	      }
      }
    }

    if(r_max>r_global_max){
      r_global_max=r_max;
      (*max0)=p0; //Save location to return to calling function.  max0
      (*max1)=ptmp; //is the earlier in the list and max1 the later.
    }
  } //End loop over p1

  if(r_global_max>0.0){
    r_global_max=(float)sqrt((double)r_global_max);
  }else{
    r_global_max=-1.0;
  }
  return r_global_max;
}



/*****************************************************/
float area_of_cell(struct point *begin){
  //begin points to a linked list of points.  Here we find all the
  //points interior to the polygon that's formed by connecting all the
  //points pointed to by begin.  This routine is essentially identical
  //to find_interior_points() except that it only calculates the number
  //of points (also allowing fractional points).
  //(I could have somehow hacked this into find_interior_points(), but they'll
  //both probably go faster like this, and it's not as ugly.)

  struct point *p0;
  struct point *p1;

  int jmin,jmax;
  int i,k;

  float tmp;

  int isect_found;
  int isect;
  float x[isect_max];
  float y[isect_max];
  float yj,dely;
  double doub_yj;

  float area;

  if(begin==NULL){
    printf("A NULL call sent to find area!\n");
    return 0.0;
  }
  if((begin->next)==NULL){
    printf("A 1-point call sent to find area!\n");
    return 0.0;
  }

  //Find min and max value of j in the point list
  jmin=ymax;
  jmax=0;
  for(p0=begin;p0!=NULL;p0=p0->next){
    if((p0->j)<jmin) jmin=p0->j;
    if((p0->j)>jmax) jmax=p0->j;
  }

  //We're going to loop from jmin to jmax and consider segments going
  //horizontally across entire picture.  We then look for intersections
  //of our line segments with this big segment.  We list the locations
  //of all found intersections.  If we write these points sorted in x as
  //i1,i2,i3,i4,i5,..., then the interior points are i1 to i2 and then
  //i3 to i4, etc.
  area=0.0;
  yj=((float)jmin)-0.1;
  while(yj<=(float)jmax){
    dely=0.1;
    yj+=dely;
    doub_yj=(double)yj;
    if(  (fabs(floor(doub_yj)-doub_yj))<0.001 ){
      //Don't want integer values of yj
      //since we'll have to deal with tangential overlaps at segments edges
      //(floor(x)=largest integer not greater than x (as a double))
      yj+=dely; //Step another dely
      dely=2.0*dely; //Will use for integral sum below
    }

    isect_found=0;
    for(p0=begin;p0!=NULL;p0=p0->next){
      p1=p0->next;
      if(p1==NULL)p1=begin; //Keep loop going around begin point

      isect=do_segments_intersect(
				  (float)p0->i,(float)p0->j,
				  (float)p1->i,(float)p1->j,
				  -1.0,yj,
				  (float) (xmax+1),yj,
				  (x+isect_found),
				  (y+isect_found) );

      if(isect==1){
	isect_found++;
	if(isect_found>isect_max){
	  printf("Too many intersections found.\n");
	  isect_found--;
	}
      }
    }//End loop over segment list for this value of yj

    //Make sure number found is even
    if((isect_found%2)!=0){
      printf("Odd number of intersection points in area search!!!!!\n");
      return 0.0;
    }

    //If we've found at least two points, then sort them in x[] so can
    //calculate interior points.
    for(i=0;i<(isect_found-1);i++){
      for(k=i+1;k<isect_found;k++){
	if(x[i]>x[k]){
	  tmp=x[i];
	  x[i]=x[k];
	  x[k]=tmp;
	  tmp=y[i];
	  y[i]=y[k];
	  y[k]=tmp;
	}
      }
    }

    //Now add all the points between points 0 and 1,  2 and 3, etc.
    for(i=0;i<(isect_found-1);i+=2){
      area+=(x[i+1]-x[i])*dely;
    }

  }//End of loop over j from jmin to jmax of segment points

  return area;

}

/************************************************************/
void fill_cos_sin_arrays(void){

  int i;
  double theta;
  free(sin_theta);
  free(cos_theta);

  sin_theta=(double *)malloc(n_points_r_vs_theta*sizeof(double));
  cos_theta=(double *)malloc(n_points_r_vs_theta*sizeof(double));

  dtheta=twopi/((double)n_points_r_vs_theta);
  for(i=0;i<n_points_r_vs_theta;i++){
    theta=dtheta*((double)i);
    cos_theta[i]=cos(theta);
    sin_theta[i]=sin(theta);
  }

}

/***********************************************************/
void r_vs_theta_from_boundary(struct point *p){
  //p is a linked list of boundary points. Use them to calculate
  //radius vs theta (starting at centroid).

  //Take a list of (x,y) points as a polygon.  Using the centroid, calculate
  //r vs theta and then calculate FFT.

  double length,total_length;
  double tmp1,tmp2,sumx,sumy;
  double x0,y0;
  double x1,y1;
  struct point *start;
  struct point *p0;
  struct point *p1;
  struct point *psave;

  int isect;

  float tmpx,tmpy;

  int i;

  if ((sin_theta==NULL)||(cos_theta==NULL)){
    fill_cos_sin_arrays();
  }

  //Find start of list
  while((p->prev)!=NULL) p=p->prev;
  start=p;

  //Calculate centroid.  The x-component of the centroid is the sum of
  //the average x positions weighted by the length of each segment, and
  //similarly for the y-component.
  p0=start;
  p1=start->next;
  if(p1==NULL){  //Only two points, not a polygon
    for(i=0;i<n_points_r_vs_theta;i++){
      r_vs_theta[i]=0.0;
    }
    return;
  }
  total_length=0.0;
  sumx=0.0;
  sumy=0.0;
  while(p0!=NULL){ //Loop all the way around
    x0=(double)( p0->i );
    x1=(double)( p1->i );
    y0=(double)( p0->j );
    y1=(double)( p1->j );
    tmp1=x1-x0;
    tmp2=y1-y0;
    length=sqrt( tmp1*tmp1 + tmp2*tmp2 );
    total_length += length;
    sumx += (length*((x0+x1)/2.0));
    sumy += (length*((y0+y1)/2.0));
    p0=p0->next;
    p1=p1->next;
    if(p1==NULL) p1=start;
  }

  x0=sumx/total_length;
  y0=sumy/total_length;
  //  printf("Centroid: (%e,%e)\n",x0,y0);
  centroid_x=((float)x0); //Save centroid position in case want to
  centroid_y=((float)y0);//alter r-vs-theta and reconstruct a boundary

  p0=start;
  p1=start->next;
  //Now transform the polygon into npoints of (r,theta) for theta from 0 to
  //twopi
  for(i=0;i<n_points_r_vs_theta;i++){
    x1=x0+cos_theta[i]*1024.0; //1024 to get out of image region
    y1=y0+sin_theta[i]*1024.0;

    //Find polygon segment that this intersects
    psave=p0; //To make sure don't go around more than once
    while((isect=do_segments_intersect(
				       (float)(p0->i), (float)(p0->j),
                                       (float)(p1->i), (float)(p1->j),
                                       (float)(x0), (float)(y0),
                                       (float)(x1), (float)(y1),
                                       &tmpx,&tmpy)
           )==0){ //0 means they don't intersect
      p0=p1;
      p1=p1->next;
      if(p1==NULL)p1=start;
      if(p0==psave){ //Didn't find any intersection for this theta.  Our
        //centroid is probably not inside the polygon.
        //printf("No intersection, returning NULL\n");
	tmpx=((float)x0);
	tmpy=((float)y0); //Return radius of zero
	break;

      }
    }
    //printf("%i %e, %e\n",i,tmpx,tmpy);
    x1=((double)tmpx)-x0;
    y1=((double)tmpy)-y0;
    r_vs_theta[i]=(float)(sqrt( x1*x1 + y1*y1));
  }

  return;
}


/***********************************************************/
struct point *boundary_from_r_vs_theta(void){

  struct point *p;
  struct point *pstart;

  int i;

  if ((sin_theta==NULL)||(cos_theta==NULL)){
    fill_cos_sin_arrays();
  }

  p=point_malloc();
  p->next=NULL;
  p->prev=NULL;
  pstart=p;

  for(i=0;i<n_points_r_vs_theta;i++){
    (p->i)=((int)(r_vs_theta[i]*((float)cos_theta[i])+centroid_x+0.5));
    (p->j)=((int)(r_vs_theta[i]*((float)sin_theta[i])+centroid_y+0.5));
    //Make next boundary point
    p->next=point_malloc();
    ((p->next)->prev)=p;
    p=p->next;
  }
  //Close up list (we've made one too many points)
  p=p->prev;
  point_free(p->next);
  p->next=NULL;

  //Return start of list
  return pstart;

}

/***********************************************************/
float calculate_volume_cone_old(struct point *p_interior){
  //Calculate volume by stripping off pixels one at a time at
  //the boundary. Assume that every new set of interior pixels that
  //gets shrunk down to a smaller boundary is included in every previous
  //set (at +-z distances).
  //I could do this with r_vs_theta, but I think it might be sensitive
  //to where I put the centroid.
  //We only need the boundary[] points because I want to include them
  //as part of the cell.
  //First fill over[] with all boundary and interior points.
  //next remove the boundary points. Next remove the interior
  //points that are on the edge. Continue doing that until there
  //are not points left. Add up 1 voxel for every point looked at
  //every time.
  //This method is actually the volume of a cone. The volume of a
  //cone is (1/3)*(Area_base*height). For the case of a sphere, the
  //volume of half the sphere is (2/3*Area_base*height). So I scale
  //my result by 2.0 at the end to match a sphere (and also an
  //ellipsoid I think).

  int ix,iy,u;
  int u1;
  struct point *p;

  int voxels_first;
  int voxels;
  int n_found;
  float f_voxels;
  unsigned char interior_value;
  unsigned char current_outer_value;

  int test_label;

  test_label=0;
  //The interior seems to include at least some of the boundary
  //points, so I'm only going to use the interior points here.
  //We look at neighboring points, so clear out the over[] array
  update_overlap_value();
  memset(over,overlap_value,(xmax_ymax)*sizeof(unsigned char));

  update_overlap_value();
  fill_overlap_array_with_point_list(p_interior);
  interior_value=overlap_value;
  current_outer_value=interior_value; //will increment below

  voxels=0;
  voxels_first=-999;
  //Loop over interior and add good guys to voxels
  //and remove boundary points


  do {
    n_found=0;
    test_label++;
    if (test_label==5) test_label=1;
    //update_overlap_value();
    //current_outer_value=overlap_value;
    //can't use update_overlap_value since it resets over[] when it
    //overlap_value reaches its maximum.
    do{
      current_outer_value++;
      if(current_outer_value>=overlap_value_max){
	current_outer_value=2;
      }
    }while(current_outer_value==interior_value);

    for(p=p_interior;p!=NULL;p=p->next){
      ix=(p->i);
      iy=(p->j);
      //check if this guy already been removed
      if((ix>=0)&&(iy>=0)&&(ix<xmax)&&(iy<ymax)){
	u=(iy*xmax+ix);
	if (over[u]!=interior_value)continue;
      }else{
	continue; //Just ignore points outside boundary.
      }
      n_found++; //Sum up all the guys which are still around

      //if(test_label==1){
      //	d[u]=found_border;
      //}else{
      //	d[u]=0;
      //}

      //Check 4 neighboring squares and 4 diagonal squares
      if ((ix+1)<xmax){
	u1=u+1;             //(x+1,y+0)
	if ((over[u1]!=interior_value)
	    &&(over[u1]!=current_outer_value)){
	  over[u]=current_outer_value;
	  continue;
	}
	/* -->Ignore diagonal squares
	if ((iy+1)<ymax){  //(x+1,y+1)
	  u1=u+1+xmax;
	  if ((over[u1]!=interior_value)
	      &&(over[u1]!=current_outer_value)){
	    over[u]=current_outer_value;
	    continue;
	  }
	}
	if ((iy-1)>=0){    //(x+1,y-1)
	  u1=u+1-xmax;
	  if ((over[u1]!=interior_value)
	      &&(over[u1]!=current_outer_value)){
	    over[u]=current_outer_value;
	    continue;
	  }
	}
	*/

      }

      if ((ix-1)>=0){
	u1=u-1;        //(x-1,y+0)
	if ((over[u1]!=interior_value)
	    &&(over[u1]!=current_outer_value)){
	  over[u]=current_outer_value;
	  continue;
	}
	/* -->Ignore diagonal squares
	if ((iy+1)<ymax){  //(x-1,y+1)
	  u1=u-1+xmax;
	  if ((over[u1]!=interior_value)
	      &&(over[u1]!=current_outer_value)){
	    over[u]=current_outer_value;
	    continue;
	  }
	}
	if ((iy-1)>=0){    //(x-1,y-1)
	  u1=u-1-xmax;
	  if ((over[u1]!=interior_value)
	      &&(over[u1]!=current_outer_value)){
	    over[u]=current_outer_value;
	    continue;
	  }
	}
	*/

      }

      if ((iy+1)<ymax){  //(x+0,y+1)
	u1=u+xmax;
	if ((over[u1]!=interior_value)
	    &&(over[u1]!=current_outer_value)){
	  over[u]=current_outer_value;
	  continue;
	}
      }
      if ((iy-1)>=0){    //(x+0,y-1)
	u1=u-xmax;
	if ((over[u1]!=interior_value)
	    &&(over[u1]!=current_outer_value)){
	  over[u]=current_outer_value;
	  continue;
	}
      }

    }

    if (voxels_first<0){
      voxels_first=n_found;
    }
    voxels+=n_found;
  }while(n_found>0);

  //Scale by 2 to correct cone to a sphere and by 2 again since we want
  //the top and bottom direction. But first subtract out the first
  //voxels calculated since those are for the original image and I don't
  //want to double them.
  voxels=(voxels-voxels_first)*4+voxels_first;

  f_voxels=((float)voxels);
  return f_voxels;
}

/***********************************************************/
int calculate_volume(struct point *p_interior,
		     struct point *p_boundary,
		     float *v_sphere_out,
		     float *v_cone_out,
		     float *v_eff_1_out,
		     float *v_eff_2_out,
		     float *v_eff_3_out,
		     float *v_eff_4_out,
		     float *v_eff_5_out,
		     float *v_eff_6_out
		     ){
  //For every pixel calculate its minimum distance to the boundary.
  //This method is actually the volume of a cone. The volume of a
  //cone is (1/3)*(Area_base*height). For the case of a sphere, the
  //volume of half the sphere is (2/3*Area_base*height). So I scale
  //my result by 2.0 at the end to match a sphere (and also an
  //ellipsoid I think).

  int ix,iy,u,uh;
  struct point *p1;
  struct point *p2;

  int ix_min,iy_min,ix_max,iy_max,xmax_cell,ymax_cell;
  float *height=NULL;
  float *height0=NULL;

  float wt;
  float ht,ht_max;

  int dx,dy;
  int dist_min2,dist_cur2;
  int dist_min2_start;

  float v_cone;
  float v_sphere;
  float v_eff1,v_eff2,v_eff3,v_eff4,v_eff5,v_eff6;

  //For calculating effective volume observed by cross section
  float del_z=0.227; //Pixel size in microns
  float sigma_a=0.2; //sigma in microns of beam width
  float sigma_b=0.6;
  float sigma_c=1.0;
  float sigma_d=1.4;
  float sigma_e=1.8;

  float gauss_weight;

  sigma_a/=del_z; //Width in pixel units
  sigma_b/=del_z;
  sigma_c/=del_z;
  sigma_d/=del_z;
  sigma_e/=del_z;
  gauss_weight=((float)(sqrt(2.0*3.1415926)));

  //During loop over interior, keep track of min and max of {ix,iy} so
  //can make another array to keep track of heights.
  ix_min=2*xmax;
  iy_min=2*ymax;
  ix_max=-xmax;
  iy_max=-ymax;

  v_cone=0.0;
  dist_min2_start=xmax_ymax;
  for(p1=p_interior;p1!=NULL;p1=p1->next){
    ix=(p1->i);
    iy=(p1->j);
    if ((ix<0)||(iy<0)||(ix>=xmax)||(iy>=ymax)) continue;
    if (ix<ix_min)ix_min=ix;
    if (iy<iy_min)iy_min=iy;
    if (ix>ix_max)ix_max=ix;
    if (iy>iy_max)iy_max=iy;

    u=(iy*xmax+ix);
    //Loop over boundary points to find minimum distance to this point
    dist_min2=dist_min2_start;
    for(p2=p_boundary;p2!=NULL;p2=p2->next){
      dx=(p2->i)-ix;
      dy=(p2->j)-iy;
      dist_cur2=dx*dx+dy*dy;
      if (dist_cur2<dist_min2)dist_min2=dist_cur2;
    }

    work_array[u]=dist_min2;

  }

  //Now combine a bunch of spheres where sphere_i is centered on point
  //i and has a radius sqrt(dist_min) for that point. The combination is
  //to set the height of each point to the maximum of all the radii.

  //Create an array to keep track of heights
  xmax_cell=(ix_max-ix_min)+1;
  ymax_cell=(iy_max-iy_min)+1;
  height=(float *)malloc(xmax_cell*ymax_cell*sizeof(float));
  height0=(float *)malloc(xmax_cell*ymax_cell*sizeof(float));
  for(uh=0;uh<xmax_cell*ymax_cell;uh++){
    height[uh]=0.0;
  }
  //add 1.0 to heights and initialize height0[] and calculate v_cone
  for(p1=p_interior;p1!=NULL;p1=p1->next){
    ix=(p1->i);
    iy=(p1->j);
    if ((ix<0)||(iy<0)||(ix>=xmax)||(iy>=ymax)) continue;
    uh=(iy-iy_min)*xmax_cell+(ix-ix_min);
    u=iy*xmax+ix;
    ht=(float)sqrt((double)work_array[u])+1.0; //Add 1.0 so boundaries aren't 0
    height0[uh]=(ht*ht); //Squared again

    //Keep track of spherical volume as determined from a cone. (Ie,
    //that's the 4x correction factors (2 to go from cone to half-sphere,
    //and 2 to go from half sphere to sphere)
    //(subtracting .5 for cone seems to work better.)
    if (ht>1.1){ //don't weight middle slice
      v_cone+=(4.0*(ht-.5));
    }else{
      v_cone+=(ht-.5);
    }
  }

  //For every pixels calculate its height^2 as the maximum over all other
  //pixels of R^2-dist^2 where dist is distance from center of the other
  //pixels and R is their height
  for(p1=p_interior;p1!=NULL;p1=p1->next){
    ix=(p1->i);
    iy=(p1->j);
    if ((ix<0)||(iy<0)||(ix>=xmax)||(iy>=ymax)) continue;
    uh=(iy-iy_min)*xmax_cell+(ix-ix_min);
    ht_max=-1.0;
    for(p2=p_interior;p2!=NULL;p2=p2->next){
      dx=(p2->i)-ix;
      dy=(p2->j)-iy;
      dist_cur2=dx*dx+dy*dy;
      u=((p2->j)-iy_min)*xmax_cell+((p2->i)-ix_min);

      ht=(float)(height0[u]-dist_cur2);
      if(ht<=0.0)ht=0.0;
      if (ht>ht_max)ht_max=ht;
    }
    height[uh]=((float)sqrt((double)ht_max));
  }

  //Sum up the heights to get volume of sphere in voxels
  //Also calculate effective volume in different ways.
  v_sphere=0.0;
  v_eff1=0.0;
  v_eff2=0.0;
  v_eff3=0.0;
  v_eff4=0.0;
  v_eff5=0.0;
  v_eff6=0.0;
  for(p1=p_interior;p1!=NULL;p1=p1->next){
    ix=(p1->i);
    iy=(p1->j);
    if ((ix<0)||(iy<0)||(ix>=xmax)||(iy>=ymax)) continue;
    uh=(iy-iy_min)*xmax_cell+(ix-ix_min);
    ht=height[uh];
    v_sphere+=ht;

    if (ht<2.0){
      v_eff1+=ht;
    }else{
      v_eff1+=2.0;
    }
    //if (ht<4.0){
    //  v_eff2+=ht;
    //}else{
    //  v_eff2+=4.0;
    //}
    //if (ht<6.0){
    //  v_eff3+=ht;
    //}else{
    //  v_eff3+=6.0;
    //}
    wt=integrate_gaussian_from_0_to_x(ht,sigma_a);
    v_eff2+=wt;
    wt=integrate_gaussian_from_0_to_x(ht,sigma_b);
    v_eff3+=wt;
    wt=integrate_gaussian_from_0_to_x(ht,sigma_c);
    v_eff4+=wt;
    wt=integrate_gaussian_from_0_to_x(ht,sigma_d);
    v_eff5+=wt;
    wt=integrate_gaussian_from_0_to_x(ht,sigma_e);
    v_eff6+=wt;

    //test-test-test-test
    //u=iy*xmax+ix;
    //work_array[u]=(int)(1000.0*height[uh]);
  }
  //Scale by two for lower half
  v_sphere*=2.0;
  v_eff1*=2.0;
  v_eff2*=(2.0*gauss_weight*sigma_a);
  v_eff3*=(2.0*gauss_weight*sigma_b);
  v_eff4*=(2.0*gauss_weight*sigma_c); //Make so pixel at z=0 gets weight 1
  v_eff5*=(2.0*gauss_weight*sigma_d);//not weight 1/sqrt(2*pi)/sigma which
  v_eff6*=(2.0*gauss_weight*sigma_e);//is what integrate_gaussian...() does

  free(height);
  free(height0);

  (*v_sphere_out)=v_sphere;
  (*v_cone_out)=v_cone;
  (*v_eff_1_out)=v_eff1;
  (*v_eff_2_out)=v_eff2;
  (*v_eff_3_out)=v_eff3;
  (*v_eff_4_out)=v_eff4;
  (*v_eff_5_out)=v_eff5;
  (*v_eff_6_out)=v_eff6;

  return 1.0;
}

/***********************************************************/
float integrate_gaussian_from_0_to_x(float x,float sigma){
  //Return integral of gaussian from 0 to x. (It's
  //the integral of (1/sqrt(2*pi)/sigma*exp(-x**2/2.0/sigma**2) from 0 to x.)
  //(so for x very large it goes to .5 and x very negative it goes to -.5)
  //Do it by simply tabulating sum of
  //1/sqrt(2*pi)*exp(-z**2/2) and then doing
  //z=x/sigma, d(x/sigma)=dz....

  double dz,z,sum,invnorm;
  int i;

  float absx;
  float integral;

  if (gaussian_integral==NULL){
    gaussian_integral=(double *)malloc(sizeof(double)*
				       gaussian_integral_nbins);
    dz=(max_gaussian_integral)/((double)gaussian_integral_nbins);
    gaussian_integral_bin_width=((float)dz);
    invnorm=1.0/sqrt(2.0*3.1415926535);
    z=dz/2.0; //start in middle of bin
    sum=0.0;
    for(i=0;i<gaussian_integral_nbins;i++){
      //      sum+=(dz*invnorm*exp(-z*z/2.0));
      sum+=(dz*invnorm/sqrt(1.0+z*z));
      gaussian_integral[i]=sum;
      z+=dz;
    }
  }

  if (sigma<=0.0){
    return 0.0;
  }
  if (x==0.0){
    return 0.0;
  }
  if (x<0.0){
    absx=-x;
  }else{
    absx=x;
  }
  absx/=sigma;

  i=(int)(absx/gaussian_integral_bin_width);
  if (i>=gaussian_integral_nbins){
    //    return 0.5; //Integral goes to .5 since is normalized so
    //integral from -inf to inf is 1.
    integral=((float)gaussian_integral[gaussian_integral_nbins-1]);
    return integral;
  }
  integral=((float)gaussian_integral[i]);
  if (x<0.0){
    integral=-integral;
  }
  return integral;
}

/***********************************************************/
void statistics_from_r_vs_theta(int flag,
				float *circ_out,
				float *major_out1,
				float *minor_out1,
				float *major_out2,
				float *minor_out2,
				float *vrot_out){

  //Using r_vs_theta find the circumference as well as the
  //major and minor axis lengths. The major axis is the longest
  //segment that goes through the center and is attached at the
  //boundaries, and the minor is the length of the segment that's
  //perpendicular to the major (not necessarily the smallest segment).

  //The first (major,minor) pair are calculated by finding largest
  //segment and then minor is 90degrees away. The second pair is
  //done by calculating smallest and then major is 90 degrees away.

  //Assume that r_vs_theta[] is already calculated.

  int i,j;
  int i_max,i_min;
  float major,minor;
  float tmp;

  double circ;
  double cos_dtheta;
  double ds;
  double r1,r2;

  double vrot;
  double dtmp;

  int n_180;
  int n_90;

  if ((sin_theta==NULL)||(cos_theta==NULL)){
    fill_cos_sin_arrays();
  }

  n_180=(n_points_r_vs_theta/2);
  n_90=(n_points_r_vs_theta/4);

  //Find largest and smallest segments
  minor=1.0e30;
  major=-999.0;
  i_max=-1;
  i_min=-1;
  for(i=0;i<n_180;i++){
    j=i+n_180;
    //j marks point 180 degrees away from i.
    //(Note that j should always be less than n_points_r_vs_theta.)
    tmp=r_vs_theta[i]+r_vs_theta[j];
    if (tmp>major){
      major=tmp;
      i_max=i;
    }
    if (tmp<minor){
      minor=tmp;
      i_min=i;
    }
  }

  //We're going to return two (major,minor) pairs. The first uses
  //the major axis found and above and makes the minor 90 degrees away.
  //The second uses the minor above and makes the major axis 90 degres
  //away.
  (*major_out1)=major;
  j=(i_max+n_90);
  if (j>=n_180){
    j=(i_max-n_90);
    if (j<0)j=0; //Just in case some integer division stuff leads here.
  }
  i=j+n_180;
  (*minor_out1)=r_vs_theta[i]+r_vs_theta[j];

  if (flag==0) return; //Only want minor_out1 and major_out1


  //And now do the other way
  (*minor_out2)=minor;
  j=(i_min+n_90);
  if (j>=n_180){
    j=(i_min-n_90);
    if (j<0)j=0; //Just in case some integer division stuff leads here.
  }
  i=j+n_180;
  (*major_out2)=r_vs_theta[i]+r_vs_theta[j];

  //And now the circumference and volume of rotation about major axis
  cos_dtheta=cos(dtheta);
  circ=0.0;
  vrot=0.0;
  for(i=0;i<n_points_r_vs_theta;i++){

    //For circumference
    j=i+1;
    if (j==n_points_r_vs_theta)j=0;

    r1=(double)r_vs_theta[i];
    r2=(double)r_vs_theta[j];
    ds=sqrt(r1*r1+r2*r2-2.0*r1*r2*cos_dtheta);
    circ+=ds;

    //For volume of rotation about major axis
    j=i-i_max; //Rotating about major axis.
    if (j<0)j+=n_points_r_vs_theta;
    dtmp=fabs(sin_theta[j]);
    vrot+=(dtmp*r1*r1*r1);
  }
  //vrot should be (assuming azimuthal symmetry)
  // 2pi/3*integral_0_pi{sin(theta)r**3dtheta} but above we did
  // integral_0_(twopi){|sin(theta)|r**3dtheta}.
  //Assuming that the cell is symmetric, then we just did twice
  //the actual integral. So scale here by pi/3 intead of 2pi/3
  vrot*=((3.1415926535/3.0)*dtheta);

  (*circ_out)=((float)circ);
  (*vrot_out)=((float)vrot);

  return;
}

/**********************************************/
int mark_distance_from_boundary(struct point *p_int, int nmax){

  //Fill the over[] array so that all the points contained in
  //the list p_int[] are marked so that each successive anulus
  //away from the border towards the interior has an decreasing value
  //and increasing away from the interior.
  //How far outside the border region is given by nmax.
  //Return 0 for success.

  struct point *pnew;
  struct point *pold;
  struct point *p;

  int over_start,over_cur,over_prev;

  int u0,u;
  int ix0,iy0;
  int ix,iy;

  int bound_x[]={-1, -1, -1, 0,  0, 1, 1, 1};
  int bound_y[]={-1,  0,  1, 1, -1,-1, 0, 1};
  int nbounds=8;
  int ibounds;

  int i_in,i_tot;
  int tot_in;

  int total_tries=5;
  int i_total_tries;

  if (p_int==NULL) return 1;
  if (nmax>10){
    //Put max here since we want the over[] array values to be
    //monotonically increasing as a function of annular distance.
    //But the overlap_value gets reset if its gets too high. So
    //we just kind of chose 10 as a reasonable and low value.
    printf("In mark_distance_from_boundary() nmax (%i) is too large.\n",
	   nmax);
    return 1;
  }
  if (nmax<0) return 1;

  for(i_total_tries=0;i_total_tries<total_tries;i_total_tries++){

    //Increment overlap value to be above around 25 (to reduce chance
    //that it's not high enough to get all the way to the middle of the
    //cell)
    update_overlap_value(); //Bigger than every other point in over[]
    while (overlap_value<25) update_overlap_value();

    //First mark the outer boundary
    fill_overlap_array_with_point_list(p_int);
    over_start=overlap_value;
    pnew=p_int; //Next list to loop over
    pold=NULL;
    //Loop over pold and save all neighboring points with
    //over[]<current overlap_value to the new list and mark
    //the over[] array with over_cur. Keep going until reach nmax.
    while((overlap_value-over_start)<nmax){ //

      //pnew points to list of recently added points.
      //Set pold to pnew and then setup pnew for a new list
      //Make sure to free everything appropriately.
      if (pold!=p_int) point_list_free(pold);
      pold=pnew;
      pnew=(struct point *)malloc(sizeof(struct point));
      pnew->prev=NULL;

      //Update overlap value for the next exterior boundary points
      over_prev=overlap_value;
      update_overlap_value();
      if (overlap_value<over_prev) break; //was reset, will start over

      for(p=pold;p!=NULL;p=p->next){
	ix0=(p->i);
	iy0=(p->j);
	if ((ix0<0)||(iy0<0)||(ix0>=xmax)||(iy0>=ymax)) continue;
	//Check all 8 boundary pixels
	u0=(iy0*xmax+ix0);
	for(ibounds=0;ibounds<nbounds;ibounds++){
	  ix=(ix0+bound_x[ibounds]);
	  iy=(iy0+bound_y[ibounds]);
	  if ((ix<0)||(iy<0)||(ix>=xmax)||(iy>=ymax)) continue;
	  u=(iy*xmax+ix);
	  if (over[u]<over_prev){
	    over[u]=overlap_value;
	    pnew->i=ix;
	    pnew->j=iy;
	    pnew->next=(struct point *)malloc(sizeof(struct point));
	    (pnew->next)->prev=pnew;
	    pnew=pnew->next;
	  }
	}
      }//Done with loop over pold points.
      //Terminate pnew
      if ((pnew->prev)!=NULL){
	pnew=pnew->prev;
	point_free(pnew->next);
	pnew->next=NULL;
      }else{ //Were no more points, so stop list
	pnew->next=NULL;
	break; //This last point gets free'd below
      }
    }
    if (pold!=p_int) point_list_free(pold);
    if (pnew!=p_int) point_list_free(pnew);

    if (overlap_value<over_start) continue;
    //The overlap_value<overstart is a hack. It avoids the case that
    //the overlap_value got reset doing the isearch++ loop. When that
    //happens the entire arrays gets re-initialized. It also ensures
    //that we can assume that the over[] array has interior-->exterior
    //points in increasing order.
    //The "continue" increments i_total_tries to try again


    //Now do a kind of similar type of thing going inward.
    over_prev=over_start;
    over_cur=over_prev-1;
    tot_in=0;
    do{
      for(p=p_int;p!=NULL;p=p->next){
	ix0=(p->i);
	iy0=(p->j);
	if ((ix0<0)||(iy0<0)||(ix0>=xmax)||(iy0>=ymax)) continue;
	//Check all 8 boundary pixels
	//Pixels with all 8 boundary pixels <= over_prev are _not_
	//on boundary, so reduce them to over_cur
	u0=(iy0*xmax+ix0);
	i_in=0;
	i_tot=0;
	for(ibounds=0;ibounds<nbounds;ibounds++){
	  ix=(ix0+bound_x[ibounds]);
	  iy=(iy0+bound_y[ibounds]);
	  if ((ix<0)||(iy<0)||(ix>=xmax)||(iy>=ymax)) continue;
	  u=(iy*xmax+ix);
	  i_tot++;
	  if (over[u]>over_prev) break;
	  i_in++;
	}
	if (i_in==i_tot){ //Was totally interior point
	  tot_in++;
	  over[u0]=over_cur;
	}
      }//Done with loop over p_int interior points.

      over_prev=over_cur;
      over_cur--;
      if (over_cur<=2) break; //Too small

    } while (tot_in>0) ;//Still have interior points to do

    if (over_cur>2){
      return 0; //done, success
    }

  } //End of the i_total_tries loop

  return 1; //Returning here means failure (i_total_tries too high)

}

/**********************************************/
void neighbor_interpolation(int ix0,int iy0,
			    float centerx,float centery,
			    float *array,
			    float *prediction, float *rms){
  //Do a 2d interpolation for pixels +-2 pixels around the
  //point (ix,iy). Return predicted value and rms of residuals.

  //Fit to f(x,y)=f0+a*x+b*y+m*xy (so actually "quadratic interpolation")
  //Set coordinates so (ix,iy) is in center.
  //We're going to do fit numerically.

#define ninterp 25
  int u,u0;
  float vinterp[ninterp*ninterp];
  float dvinterp[ninterp*ninterp];
  float xinterp[ninterp*ninterp];
  float yinterp[ninterp*ninterp];
  float xyinterp[ninterp*ninterp];
  float xxinterp[ninterp*ninterp];
  float yyinterp[ninterp*ninterp];
  int nuse;
  int i,j;

  int npar=6;
  float pars[6];
  float step[]={4.0,2.0,5.0,4.0,2.0,2.0};
  float eps[]={.1,.01,.01,.01,.01,.01};

  float stepcur;
  int dir;
  float sum0,sum;
  float fpred;
  float tmp;
  float vcenter;

  int iloop,nloop;
  int idx,idy;

  float res0;
  int n0;
  int again=4; //Repeat "again" residuals, repeating each time

  double dx,dy,dr,dc,ds;
  float theta,theta0;
  float r0;

  int over0,dover;
  int iconverge;
  int nconverge=10000;

  u0=(iy0*xmax+ix0);
  vcenter=array[u0];

  //Radius and theta of initial point
  dx=((double)(ix0)-((double)centerx));
  dy=((double)(iy0)-((double)centery));
  dr=sqrt(dx*dx+dy*dy);
  dc=dx/dr;
  ds=dy/dr;
  theta0=(float)acos(dc);
  if (ds<0.0) theta0=-theta0;
  over0=over[u0];
  r0=((float)over0);

  nuse=0;
  for(i=(ix0-5);i<=(ix0+5);i++){
    idx=i-ix0;
    for(j=(iy0-5);j<=(iy0+5);j++){
      if ((i==0)&&(j==0)) continue; //Skip middle point
      idy=j-iy0;
      //if ((idx<=1)&&(idx>=-1)&&(idy<=1)&&(idy>=-1)) continue; //remove middle box

      u=(j*xmax+i);
      if ((u<0)||(u>=xmax_ymax)) continue; //out of bounds
      if (over[u]==0) continue;//0 is used to mark bad points

      //Calculate theta with respect to centroid
      dx=((double)(i)-((double)centerx));
      dy=((double)(j)-((double)centery));
      dr=sqrt(dx*dx+dy*dy);
      dc=dx/dr;
      ds=dy/dr;
      theta=(float)acos(dc);
      if (ds<0.0) theta=-theta;

      //Radial location
      dover=over[u]-over0;
      if ((dover>1)||(dover<-1))continue;
      xinterp[nuse]=((float)dover);

      //delta(theta)
      tmp=((float)theta)-theta0;
      if (tmp>6.0){
	tmp-=(twopi);
      }else if(tmp<-6.0){
	tmp+=(twopi);
      }
      yinterp[nuse]=tmp;

      xyinterp[nuse]=(xinterp[nuse]*yinterp[nuse]);
      xxinterp[nuse]=(xinterp[nuse]*xinterp[nuse]);
      yyinterp[nuse]=(yinterp[nuse]*yinterp[nuse]);
      vinterp[nuse]=array[u]-vcenter; //Subtract off central point
      if (array[u]>ccd_floor){
	dvinterp[nuse]=1.0/(array[u]-ccd_floor);
      }else{
	dvinterp[nuse]=0.0;
      }
      nuse++;
    }
  }
  //Only do x-fit
  //npar=3;

 start:

  if (nuse<(npar+3)){
    (*prediction)=-999.0;
    (*rms)=-999.0;
    return;
  }
  //Now loop over the parameters and step back and forth for each
  //separately, minimizing sum of squares
  pars[0]=step[0];
  pars[1]=0.0;
  pars[2]=0.0;
  pars[3]=0.0;
  pars[4]=0.0;
  pars[5]=0.0;
  //Initial calculation of chi^2
  sum0=0.0;
  for(j=0;j<nuse;j++){
    fpred=(pars[0]+
	   pars[1]*xinterp[j]+
    	   pars[2]*xxinterp[j]+
	   pars[3]*yinterp[j]+
	   pars[4]*xyinterp[j]+
  	   pars[5]*yyinterp[j]);
    tmp=(fpred-vinterp[j]);
    sum0+=(tmp*tmp*dvinterp[j]);
    //if (tmp<0.0) tmp=-tmp;
    //sum0+=tmp;
  }


  nloop=5;
  for(iloop=0;iloop<nloop;iloop++){

    for (i=0;i<npar;i++){
      dir=1; //0=right, 1=left
      stepcur=step[i];
      for(iconverge=0;iconverge<nconverge;iconverge++){

	//Do step
	if (dir==1){
	  pars[i]-=stepcur;
	}else{
	  pars[i]+=stepcur;
	}

	sum=0.0;
	for(j=0;j<nuse;j++){
	  fpred=(pars[0]+
		 pars[1]*xinterp[j]+
		 pars[2]*xxinterp[j]+
		 pars[3]*yinterp[j]+
		 pars[4]*xyinterp[j]+
		 pars[5]*yyinterp[j]);
	  tmp=(fpred-vinterp[j]);
	  sum+=(tmp*tmp*dvinterp[j]);
	  //if (tmp<0.0) tmp=-tmp;
	  //sum+=tmp;
	}
	//      printf("%i %e %e %e %e: %e %e\n",
	//     i,pars[0],pars[1],pars[2],pars[3],sum,sum0);

	//Did we improve?
	if (sum>sum0){ //no improvement, switch direction and cut step in half
	  if (dir==1){ //We've been going to the left
	    pars[i]+=stepcur; //Undo previous step
	  }else{
	    pars[i]-=stepcur;
	  }
	  dir=1-dir; //Switch direction
	  stepcur/=2.0; //Cut step inhalf
	  if (stepcur<eps[i]) break; //Check if done
	}else{
	  sum0=sum; //Save for next step
	}

      }
      if (iconverge>=nconverge){
	printf("No convergence: %i\n",i);fflush(stdout);
	return;
      }

    }
  }

  if (again>0){
    again--;
    //Remove highest residual and repeat
    res0=-1.0;
    n0=-1;
    for(j=0;j<nuse;j++){
      fpred=(pars[0]+
	     pars[1]*xinterp[j]+
	     pars[2]*xxinterp[j]+
	     pars[3]*yinterp[j]+
	     pars[4]*xyinterp[j]+
	     pars[5]*yyinterp[j]);
      tmp=(fpred-vinterp[j]);
      //if (tmp<0.0) tmp=-tmp;
      //sum=tmp;
      sum=(tmp*tmp*dvinterp[j]);
      if (sum>res0){ //res1<res0=max
	n0=j;
	res0=sum;
      }
    }
    //Remove (n0)
    xinterp[n0]=xinterp[nuse-1];
    yinterp[n0]=yinterp[nuse-1];
    vinterp[n0]=vinterp[nuse-1];
    dvinterp[n0]=dvinterp[nuse-1];
    nuse--;
    goto start;
  }


  (*prediction)=(pars[0]+vcenter);
  //(*prediction)=(float)sqrt((double)
  //			    (pars[1]*pars[1]+pars[2]*pars[2]));
  //Calculate sum of weights to get rms
  sum=0.0;
  for(j=0;j<nuse;j++){
    sum+=(dvinterp[j]);
  }
  (*rms)=((float)sqrt(((double)sum0)/((double)sum)));
  //(*rms)=(((float)sum0)/((float)nuse));

  //printf("pars: %e %e %e %e %e %e\n",
  //	 pars[0],pars[1],pars[2],pars[3],pars[4],pars[5]);


  return;
}

/**********************************************/
int get_median(int ix, int iy,int n,float *array,float *smallest,
		 float *largest, float *median){
  //Get median of neighboring points to (ix,iy). Include the
  //8 neighboring points, but only if they're part of the same group
  //as the point (ix,iy).

  int u,u2;
  int i,j;
  int k,l;
  int iy2;

  float val;
  float sorted[50];
  int inext;
  int ind;
  int overlap_group;

  //we do +=n points
  if ((n!=1)&&(n!=2))return 1;

  inext=0;
  u=(iy*xmax+ix);

  overlap_group=over[u];

  for(j=-n;j<=n;j++){
    iy2=(iy+j);
    u2=iy2*xmax+(ix-(n+1)); //Will increment in i loop
    for(i=-n;i<=n;i++){
      u2++;
      if ((j==0)&&(i==0)) continue;
      if ((u2<0)||(u2>=xmax_ymax)) continue;
      if (over[u2]!=overlap_group) continue;
      val=array[u2];

      //Where to put it
      ind=inext;
      for(k=0;k<inext;k++){
	if (val<sorted[k]){
	  ind=k;
	  break;
	}
      }
      //Put val in the ind-th position, move everyone else up
      for(l=inext;l>ind;l--){
	sorted[l]=sorted[l-1];
      }
      sorted[ind]=val;
      inext++;
    }
  }

  //inext is how many points we found
  if (inext<=1) return 1; //flag an error

  (*smallest)=sorted[0];
  (*largest)=sorted[(inext-1)];

  //Get median
  i=(inext/2); //integer division
  j=(inext%2); //Remainder of integer division
  if (j==0){ //inext was even
    val=(sorted[i-1]+sorted[i])/2.0;
  }else{
    val=sorted[i];
  }

  (*median)=val;
  return 0;
}

/**********************************************/
int get_median_deviations(int ix, int iy,float y0,float *array,float *dev){
  //Get median of neighboring points to (ix,iy) of
  //absolute value of deviations from the value y0.

  int u,u2;
  int i,j;
  int k,l;
  int iy2;

  float val;
  float sorted[50];
  int inext;
  int ind;
  int overlap_group;

  float median;

  inext=0;
  u=(iy*xmax+ix);

  i=get_median(ix,iy,2,array,&val,&val,&median);
  if (i>0) return 1; //error
  //median=y0;

  overlap_group=over[u];
  for(j=-2;j<=2;j++){
    iy2=(iy+j);
    u2=iy2*xmax+(ix-3); //Will increment in i loop
    for(i=-2;i<=2;i++){
      u2++;
      if ((j==0)&&(i==0)) continue;
      if ((u2<0)||(u2>=xmax_ymax))continue;
      if (over[u2]!=overlap_group) continue;
      val=fabsf(array[u2]-median);

      //Where to put it
      ind=inext;
      for(k=0;k<inext;k++){
	if (val<sorted[k]){
	  ind=k;
	  break;
	}
      }
      //Put val in the ind-th position, move everyone else up
      for(l=inext;l>ind;l--){
	sorted[l]=sorted[l-1];
      }
      sorted[ind]=val;
      inext++;
    }
  }

  //inext is how many points we found
  if (inext<=1) return 1; //Flag an error

  //Get median
  i=(inext/2); //integer division
  j=(inext%2); //Remainder of integer division
  if (j==0){ //inext was even
    val=(sorted[i-1]+sorted[i])/2.0;
  }else{
    val=sorted[i];
  }
  (*dev)=val;
  return 0;
}





/***********************************************************/
void sort_list_by_fl(struct point *pin){
  //Sort the list based on the value of fl[][]
  //(So assuming that fl[][] has valid data, etc.)

  struct point *pstart;
  struct point *p0;
  struct point *p1;
  float fl0,fl1;
  float tmp;
  int i,j;

  //Go to start of list
  pstart=pin;
  while ((pstart->prev)!=NULL)pstart=(pstart->prev);

  for(p0=pstart;(p0->next)!=NULL;p0=(p0->next)){
    fl0=fl[(p0->j)*xmax+(p0->i)];
    //fl0=cur_prev[p0->i][p0->j];
    for(p1=p0->next;p1!=NULL;p1=p1->next){
      fl1=fl[(p1->j)*xmax+(p1->i)];
      //fl1=cur_prev[p1->i][p1->j];
      if (fl1<fl0){
	i=p0->i;
	j=p0->j;
	(p0->i)=(p1->i);
	(p0->j)=(p1->j);
	(p1->i)=i;
	(p1->j)=j;
	tmp=fl1;
	fl1=fl0;
	fl0=tmp;
      }
    }
  }

  return;
}


/**************************************************************/
struct point *copy_cell_for_split_regions(
					  struct point *p_in,
					  int flag){

  struct point *ptmp;
  struct point *p;
  int fret_offx;
  int fret_offy;
  float tmpx;
  struct point *p_return;

  //Assume fret_bx, etc, already calculated. Use them to move pixel
  //list by linear offset function.

  if (p_in==NULL) return NULL;

  p_return=point_malloc();
  p_return->prev=NULL;

  ptmp=p_return; //start of new list

  if (flag==0){ //Go from bottom to top
    for(p=p_in;p!=NULL;p=p->next){
      tmpx=((float)(p->i));
      fret_offx=((int)(fret_bx+fret_mx*tmpx+0.5));
      fret_offy=((int)(fret_by+fret_my*tmpx+0.5));
      (ptmp->i)=(p->i)+fret_offx; //Note it may not be in image!
      (ptmp->j)=(p->j)+fret_offy;
      ptmp->next=point_malloc();
      (ptmp->next)->prev=ptmp;
      ptmp=ptmp->next;
    }
  }else if (flag==1){ //Go from top to bottom
    //Inverting the fret_offx,y function for going from top to bottom
    for(p=p_in;p!=NULL;p=p->next){
      tmpx=(((float)(p->i))-fret_bx)/(1.0+fret_mx);
      fret_offy=((int)(fret_by+fret_my*tmpx+0.5));
      (ptmp->i)=((int)(tmpx+0.5));
      (ptmp->j)=(p->j)-fret_offy;
      ptmp->next=point_malloc();
      (ptmp->next)->prev=ptmp;
      ptmp=ptmp->next;
    }

  }else{
    //Unknown flag type
    return NULL;
  }

  //One too many points, free last
  if ((ptmp->prev)!=NULL){
    (ptmp->prev)->next=NULL;
  }else{
    printf("Pointer screw-up in copy_cell_for_split_regions().\n");
    perror(0); //exit(0); //http://r-pkgs.had.co.nz/src.html
  }
  point_free(ptmp); //One too many points

  return p_return;
}

/**************************************************************/
// First argument, array to load
// type==0          c
// type==1          fl
// type==2          third_image
// type==3          d
// type==4          fret_labels

void load_global_arrays(int type, float *v1, int *v2,int xmax_in, int ymax_in){

  xmax=xmax_in;
  ymax=ymax_in;
  xmax_ymax=xmax*ymax;
  over_array_max=xmax_ymax;
  if (work_array==NULL){
    work_array=(int *)malloc(xmax_ymax*sizeof(int));
  }
  if (over==NULL){
    update_overlap_value(); //Will malloc over() if over==NULL
  }
  if (type==0){
    c=v1;
  }else if (type==1){
    fl=v1;
  }else if (type==2){
    third_image=v1;
  }else if (type==3){
    d=v2;
  }else if (type==4){
    fret_labels=v2;
  }else{
    printf("Unknown global array type in load_global_arrays()\n");
    fflush(stdout);
    perror(0); //exit(0); //http://r-pkgs.had.co.nz/src.html
  }
}

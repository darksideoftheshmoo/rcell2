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
#include <string.h>
#include <math.h>

#include "fl_dist.h"
#include "contiguous.h"

struct radius_offsets {
  float radius;
  int *off_x;
  int *off_y;
  int n_off;
  struct radius_offsets *prev;
  struct radius_offsets *next;
};

struct radius_offsets *radius_list_start=NULL;

int *work=NULL;
int work_xmax=0;
int work_ymax=0;

#define max_points_in_list 50000
int list_x[max_points_in_list];
int list_y[max_points_in_list];

int label_cur=0;
#define max_label 30000

/****************************************************************/
void maximum_pixels_within_fixed_radius(float *array,
					int xmax,
					int ymax,
					struct point *p_in,
					float radius,
					int **list_x_out,
					int **list_y_out,
					int *n_points_out,
					struct point **p_center){
  //Vary the center of a circle of size radius and find the center that
  //has the highest total fluorescence per pixel. Discard points that
  //are outside the cell boundaries. It will return the center point that
  //it found in p_center.

  //However, I added the feature that sometimes it doesn't search for
  //a new center, but just uses the center that gets passed in.
  //If p_center==NULL, then it will do a search for a new center. If
  //p_center!=NULL, it will skip the search and go ahead and use that
  //point as the center.

  //If the radius value is <= 0.0, then consider it as an infinite radius,
  //and simply return all the points. Ignore p_center in this case.

  int i,j,k;
  int u;
  float tmp;
  float r2;

  int *off_x;
  int *off_y;
  int n_off;

  int ix,iy;
  int ix_center,iy_center;
  float array_max;
  float sum;
  int n_sum;

  struct radius_offsets *r_use;

  struct point *p;
  struct point *pstart;
  struct point *p_max;

  (*list_x_out)=list_x;
  (*list_y_out)=list_y;
  (*n_points_out)=0; //In case any errors, return 0 for this value

  if (radius<=0.0){ //In this case, simply return the list of pixels
    for(pstart=p_in;(pstart->prev)!=NULL;pstart=pstart->prev); //Start of list
    n_sum=0;
    for(p=pstart;p!=NULL;p=(p->next)){
      ix=(p->i);
      iy=(p->j);
      if ((ix>=0)&&(ix<xmax)&&(iy>=0)&&(iy<ymax)){
      //Add a new point
	list_x[n_sum]=ix;
	list_y[n_sum]=iy;
	n_sum++;
	if (n_sum>=max_points_in_list){
	  //	  printf("Too many points in fl_dist! (%i >= %i).\n",
	  // n_sum,max_points_in_list);
	  n_sum--;
	}
      }
    }
    (*n_points_out)=n_sum;
    return;
  }

  //If radius is positive, then continue with normal stuff.

  label_cur++; //To label the current set of pixels.
  if (label_cur>max_label){
    if (work!=NULL){
      memset(work,0,(xmax*ymax)*sizeof(int));
      label_cur=1;
    }else{
      printf("Too many cell labels (%i), but no work array!!!\n",label_cur);
      perror(0); //exit(0); //http://r-pkgs.had.co.nz/src.html
    }
  }

  //First check if we have already made a set of pixel offsets of the
  //current radius (which we'll use to loop in circle of size radius quickly)
  for(r_use=radius_list_start;r_use!=NULL;r_use=(r_use->next)){
    if ((r_use->radius)==radius) break;
  }
  if (r_use==NULL){ //Didn't find any location, have to make a new
      //set of offsets
    if (radius_list_start!=NULL){
      r_use=radius_list_start;
      while((r_use->next)!=NULL) r_use=(r_use->next);
      //Now at last position, make a new one
      (r_use->next)=(struct radius_offsets *)
	malloc(sizeof(struct radius_offsets));
      ((r_use->next)->prev)=r_use;
      r_use=(r_use->next);
      r_use->next=NULL;
    }else{ //First time
      radius_list_start=(struct radius_offsets *)
	malloc(sizeof(struct radius_offsets));
      r_use=radius_list_start;
      (r_use->next)=NULL;
      (r_use->prev)=NULL;
    }
    //Now set up the lists
    r2=radius*radius;
    //First count up how many points, then make a list of the points
    n_off=0;
    for(i=-250;i<250;i++){
      for(j=-250;j<250;j++){
	tmp=(float)(i*i+j*j);
	if (tmp<=r2) n_off++;
      }
    }
    //Make an array of the right size
    off_x=(int *)malloc(n_off*sizeof(int));
    off_y=(int *)malloc(n_off*sizeof(int));
    k=0;
    for(i=-250;i<250;i++){
      for(j=-250;j<250;j++){
	tmp=(float)(i*i+j*j);
	if (tmp<=r2){
	  if (k<n_off){
	    off_x[k]=i;
	    off_y[k]=j;
	  }
	  k++;
	}
      }
    }
    if (k!=n_off){
      printf("------> Inconsistency in setting up arrays in fl_dist.\n");
      if ((r_use->prev)!=NULL){
	r_use=(r_use->prev);
	free(r_use->next);
	r_use->next=NULL;
	return; //no points
      }
    }else{
      (r_use->radius)=radius;
      (r_use->n_off)=n_off;
      (r_use->off_x)=off_x;
      (r_use->off_y)=off_y;
    }
  }else{ //if we already have some offsets setup
    off_x=(r_use->off_x);
    off_y=(r_use->off_y);
    n_off=(r_use->n_off);
  }

  //Go to the start of the points that we're looking at
  for(pstart=p_in;(pstart->prev)!=NULL;pstart=pstart->prev);

  //Set up the work array to label where the cell is
  if ((work_xmax!=xmax)||(work_ymax!=ymax)){
    free(work);
    work=(int *)malloc(xmax*ymax*sizeof(int));
    work_xmax=xmax;
    work_ymax=ymax;
    memset(work,0,(xmax*ymax)*sizeof(int));
  }

  //Label work array with the current points
  for(p=pstart;(p!=NULL);p=p->next){
    ix=p->i;
    iy=p->j;
    if ((ix>=0)&&(ix<xmax)&&(iy>=0)&&(iy<ymax)){
      work[ iy*xmax + ix ]=label_cur;
    }
  }

  if ((*p_center)==NULL){ //Search for new center
    //Begin searching for center of cell which has maximum for the offsets
    p_max=NULL;
    array_max=-999.0;
    for (p=pstart;(p!=NULL);p=p->next){
      ix_center=(p->i);
      iy_center=(p->j);
      sum=0.0;
      n_sum=0;
      for(k=0;k<n_off;k++){
	ix=ix_center+off_x[k];
	iy=iy_center+off_y[k];
	if ((ix>=0)&&(ix<xmax)&&(iy>=0)&&(iy<ymax)){
	  u=(iy*xmax+ix);
	  if (work[u]==label_cur){
	    sum+=(array[u]);
	    n_sum++;
	  }
	}
      }
      if (n_sum>0){
	tmp=(sum/((float)n_sum));
	if (tmp>array_max){
	  array_max=tmp;
	  p_max=p; //Save location
	}
      }
    }

    if (p_max==NULL){ //Didn't find any points
      return;
    }
    (*p_center)=(struct point *)malloc(sizeof(struct point));
    (*p_center)->next=NULL;
    (*p_center)->prev=NULL;
    (*p_center)->i=p_max->i; //For returning center location
    (*p_center)->j=p_max->j; //For returning center location
  }else{
    p_max=(*p_center);
  }

  //Now make a list of all the points starting at this point of this radius
  //which are within the cell.
  n_sum=0;
  ix_center=(p_max->i);
  iy_center=(p_max->j);
  for(k=0;k<n_off;k++){
    ix=ix_center+off_x[k];
    iy=iy_center+off_y[k];
    if ((ix>=0)&&(ix<xmax)&&(iy>=0)&&(iy<ymax)){
      u=(iy*xmax+ix);
      if (work[u]==label_cur){
	//Add a new point
	list_x[n_sum]=ix;
	list_y[n_sum]=iy;
	n_sum++;
	if (n_sum>=max_points_in_list){
	  printf("Too many points in fl_dist! (%i >= %i).\n",
		 n_sum,max_points_in_list);
	  n_sum--;
	}
      }
    }
  }

  (*n_points_out)=n_sum;
  return;

}

/************************************************************/
void maximum_contiguous_pixels(float *array,
			       int xmax,
			       int ymax,
			       struct point *p_in,
			       int **list_x_out,
			       int **list_y_out,
			       int *n_points_out){

  //Find regions of contiguous pixels above a cut within the
  //interior points pointed at by p_in. Return list in
  //list_x,list_y

  struct contiguous_search c_contiguous;
  struct contiguous_search *csearch;
  struct point *p;
  struct point *pstart;

  int i,j;
  int u;
  float array_max;
  float array_mean;
  int n_points;
  int n_regions;
  int n_contiguous,n_contiguous_max;
  float cut;
  float cut_step;
  int offset;
  int ix,iy;

  (*list_x_out)=list_x;
  (*list_y_out)=list_y;
  (*n_points_out)=0; //In case any errors, return 0 for this value

  label_cur++; //To label the current set of pixels.
  if (label_cur>max_label){
    if (work!=NULL){
      memset(work,0,(xmax*ymax)*sizeof(int));
      label_cur=1;
    }else{
      printf("Too many cell labels (%i), but no work array!!!\n",label_cur);
      perror(0); //exit(0); //http://r-pkgs.had.co.nz/src.html
    }
  }

  *(n_points_out)=0; //In case return with no points

  csearch=&(c_contiguous);

  (csearch->data_array)=array;
  (csearch->xmax)=xmax;
  (csearch->ymax)=ymax;
  //Setup the work array to label where the interior points are
  if ((work_xmax!=xmax)||(work_ymax!=ymax)){
    free(work);
    work=(int *)malloc(xmax*ymax*sizeof(int));
    work_xmax=xmax;
    work_ymax=ymax;
    memset(work,0,(xmax*ymax)*sizeof(int));
    label_cur=1;
  }

  //Rewind point list
  for(pstart=p_in;(pstart->prev)!=NULL;pstart=(pstart->prev));

  //Label work array with the current points (also find minimum value
  //of array)
  array_max=-999.0;
  array_mean=0.0;
  n_points=0;
  for(p=pstart;(p!=NULL);p=p->next){
    ix=p->i;
    iy=p->j;
    if ((ix>=0)&&(ix<xmax)&&(iy>=0)&&(iy<ymax)){
      u=iy*xmax+ix;
      work[u]=label_cur;
      array_mean+=array[u];
      n_points++;
      if (array[u]>array_max){
	array_max=array[u];
      }
    }
  }
  if (n_points>0){
    array_mean/=((float)n_points);
  }else{
    return; //No points
  }


  (csearch->p)=pstart; //Will only check these points
  (csearch->label_array)=work;
  (csearch->label_value)=label_cur;

  //Vary cut and find how many contiguous regions we have above the
  //current cut value. (We're going to find the cut that maximizes the
  //_number_ of contiguous regions below the cut).

  //Set the step size to (mean-min)/npixels
  n_points=0;
  for(p=pstart;(p!=NULL);p=p->next){
    ix=(p->i);
    iy=(p->j);
    if ((ix>=0)&&(ix<xmax)&&(iy>=0)&&(iy<ymax)){
      u=iy*xmax+ix;
      if ((array[u]>=array_mean)&&(array[u]<=array_max)){
	n_points++;
      }
    }
  }
  if (n_points==0){ //Distribution must be uniform, probably don't
    //want to do anything
    return;
  }

  cut_step=(array_max-array_mean)/((float)n_points);

  (csearch->cut_behavior)=equal_to_labels_cut_below;
  cut=array_max;
  do{
    cut-=cut_step;
    (csearch->cut_low)=cut;
    do_contiguous_search(csearch);
    n_regions=(csearch->n_lists_found);
    //printf("%e-->%i\n",cut,n_regions);fflush(stdout);
    if (cut<0.0) break;
  }while(n_regions==1) ;

  cut+=(cut_step);
  (csearch->cut_behavior)=equal_to_labels_cut_above;
  (csearch->cut_high)=cut;
  do_contiguous_search(csearch);

  //Find largest region (since we've gone 1 beyond n_regions=1 cut)
  n_contiguous_max=-999;
  j=-999;
  for(i=0;i<(csearch->n_lists_found);i++){
    n_contiguous=(csearch->npoints_in_list)[i];
    if(n_contiguous>n_contiguous_max){
      n_contiguous_max=n_contiguous;
      j=i;
    }
  }
  if (n_contiguous_max>max_points_in_list){
    printf("Too many contiguous points: %i > %i (setting to %i).\n",
	   n_contiguous_max,max_points_in_list,n_contiguous_max);
    n_contiguous_max=max_points_in_list;
  }
  if (n_contiguous_max>0){
    offset=(csearch->list_start)[j];
    for(i=0;i<n_contiguous_max;i++){
      list_x[i]=*((csearch->list_found_x)+offset+i);
      list_y[i]=*((csearch->list_found_y)+offset+i);
    }
  }
  (*n_points_out)=n_contiguous_max;
  printf("%i\n",n_contiguous_max);

  return;
}


/************************************************************/
void get_gauss_2d_parameters(float *array,int xmax,int ymax,
			     float back,
                             struct point *p_in,
                             float *mean_x,float *mean_y,
			     float *sig,
			     int **list_x_out,
			     int **list_y_out,
			     int *n_points_out){
  //We're going to do a 2d-gaussian+background "fit" to the data.
  //array, xmax, ymax are the data array we're looking at.
  //p_in is a list of points in the array that we're going to use.
  //p_center is a first guess at mean_x and mean_y
  //back is the assumed level of the background. (We're just going to
  //go ahead and use whatever is passed in).

  //For the calculation we just do mean_x=average x, mean_y=average y
  //and sig=rms(r) for all points above the background

  double back_cut;
  double mu_x,mu_y;
  double f;
  int n;
  double wt;
  struct point *p;

  int ix,iy,u;

  double dr2,half_inv_sig2;
  double mu2_x,mu2_y;
  double dtmp1,dtmp2;
  double sig2;

  double sig_r,sig2_r;
  double sig_x,sig_y,sig2_x,sig2_y;

  *(mean_x)=0.0; //just in case return with error
  *(mean_y)=0.0;
  *(sig)=1.0e4; //Essentially flat

  (*list_x_out)=list_x;
  (*list_y_out)=list_y;
  (*n_points_out)=0; //In case any errors, return 0 for this value

  if (p_in==NULL){
    return;
  }

  if (back<=0.0){
    back_cut=0.0;
  }else{
    //back_cut=((double)back)+2.0*(sqrt((double)back)); //back+2 sigma
    back_cut=((double)back);
  }

  //Now do weighted average+rms for points above cut
  n=0; //number of points in averages
  wt=0.0; //sum of weights in averages
  mu_x=0.0;
  mu_y=0.0;
  sig_x=0.0;
  sig_y=0.0;
  for(p=p_in;p!=NULL;p=(p->next)){
    ix=(p->i);
    iy=(p->j);
    if ((ix>=0)&&(ix<xmax)&&(iy>=0)&&(iy<ymax)){
      u=(iy*xmax+ix);
      f=((double)array[u]);
      if (f>back_cut){
	f-=back_cut;
	n++;
	wt+=f;
	mu_x+=(((double)ix)*f);
	mu_y+=(((double)iy)*f);
	sig_x+=(((double)(ix*ix))*f);
	sig_y+=(((double)(iy*iy))*f);
      }
    }
  }
  if (n<3){ //Some minimum number of points
    return;
  }
  mu_x/=wt;
  mu_y/=wt;
  sig_x=(sig_x/wt-(mu_x*mu_x));
  sig_y=(sig_y/wt-(mu_y*mu_y));
  if ((sig_x<=0.0)||(sig_y<=0.0)){
    return;
  }else{
    sig_x=sqrt(sig_x);
    sig_y=sqrt(sig_y);
    //Assume it's circular so sig=sig_x=sig_y. But do an average
    //since it won't be perfect.
    sig_r=(sig_x+sig_y)/2.0;
  }
  //printf("First:  ---------------> (%e,%e,%e)\n",
  // mu_x,mu_y,sig_r);fflush(stdout);

  //Now as a second pass, do a weighted average based on first estimate
  //of gaussian.

  half_inv_sig2=0.5/(sig_r*sig_r);

  n=0; //number of points in averages
  wt=0.0; //sum of weights in averages
  mu2_x=0.0;
  mu2_y=0.0;
  sig2_x=0.0;
  sig2_y=0.0;
  for(p=p_in;p!=NULL;p=(p->next)){
    ix=(p->i);
    iy=(p->j);
    if ((ix>=0)&&(ix<xmax)&&(iy>=0)&&(iy<ymax)){
      u=(iy*xmax+ix);
      f=((double)array[u]);
      if (f>back_cut){
	dtmp1=((double)ix)-mu_x;
	dtmp2=((double)iy)-mu_y;
	dr2=(dtmp1*dtmp1+dtmp2*dtmp2);
	f*=exp(-dr2*half_inv_sig2);
	n++;
	wt+=f;
	mu2_x+=(((double)ix)*f);
	mu2_y+=(((double)iy)*f);
	sig2_x+=(((double)(ix*ix))*f);
	sig2_y+=(((double)(iy*iy))*f);
      }
    }
  }
  if (n<3){ //Some minimum number of points
    return;
  }
  mu2_x/=wt;
  mu2_y/=wt;
  sig2_x=(sig2_x/wt-(mu2_x*mu2_x));
  sig2_y=(sig2_y/wt-(mu2_y*mu2_y));
  if ((sig2_x<=0.0)||(sig2_y<=0.0)){
    return;
  }else{
    sig2_x=sqrt(sig2_x);
    sig2_y=sqrt(sig2_y);
    sig2_r=(sig2_x+sig2_y)/2.0;
  }


  //Finally, scale by 1.5 to get 67% probability roughly
  sig2_r*=1.5;

  //printf("second:    ------->(%e,%e,%e)\n",
  //mu2_x,mu2_y,sig2_r);fflush(stdout);


  *(mean_x)=((float)mu2_x);
  *(mean_y)=((float)mu2_y);
  *(sig)=((float)sig2_r);

  //Now make a list of points within one sigma of mean.
  sig2=sig2_r*sig2_r*2.0;
  n=0;
  for(p=p_in;p!=NULL;p=(p->next)){
    ix=(p->i);
    iy=(p->j);
    if ((ix>=0)&&(ix<xmax)&&(iy>=0)&&(iy<ymax)){
      u=(iy*xmax+ix);
      dtmp1=((double)ix)-mu2_x;
      dtmp2=((double)iy)-mu2_y;
      dr2=(dtmp1*dtmp1+dtmp2*dtmp2);
      if (dr2<=sig2){
	list_x[n]=ix;
	list_y[n]=iy;
	n++;
      }
    }
  }
  (*n_points_out)=n;


  return;

}


/************************************************************/
void find_vacuole(struct point *p_in,
		  float *data,
		  int xmax, int ymax,
		  int **list_x_out,
		  int **list_y_out,
		  int *n_points_out){
  //Search through points p_interior to find a contiguous set
  //of low points that we're going to identify as the "vacuole."
  //The data is passed in data and is size (xmax,ymax).


  //For contiguous searches
  struct contiguous_search c_contiguous;
  struct contiguous_search *csearch;

  struct point *p;
  struct point *p1;
  struct point *p2;
  float val1,val2;
  float cut_min;
  struct point *pmin;

  int ix,iy,u;
  int imax,offset;
  int istep;

  int i,j,n_contiguous,nmax,n_total;
  int nmax1,nmax2;
  float frac,frac_min;

  csearch=&(c_contiguous);

  (csearch->data_array)=data;
  (csearch->xmax)=xmax;
  (csearch->ymax)=ymax;

  (*list_x_out)=list_x;
  (*list_y_out)=list_y;
  (*n_points_out)=0; //In case any errors, return 0 for this value

  //Mark interior points
  //Set up the work array to label where the cell is
  if ((work_xmax!=xmax)||(work_ymax!=ymax)){
    free(work);
    work=(int *)malloc(xmax*ymax*sizeof(int));
    work_xmax=xmax;
    work_ymax=ymax;
    memset(work,0,(xmax*ymax)*sizeof(int));
  }
  label_cur++; //To label the current set of pixels.
  if (label_cur>max_label){
    if (work!=NULL){
      memset(work,0,(xmax*ymax)*sizeof(int));
      label_cur=1;
    }else{
      printf("Too many cell labels (%i), but no work array!!!\n",label_cur);
      perror(0); //exit(0); //http://r-pkgs.had.co.nz/src.html
    }
  }
  //Label work array with the current points
  n_total=0;
  for(p=p_in;(p!=NULL);p=p->next){
    ix=p->i;
    iy=p->j;
    if ((ix>=0)&&(ix<xmax)&&(iy>=0)&&(iy<ymax)){
      work[ iy*xmax + ix ]=label_cur;
      n_total++;
    }
  }

  if (n_total<10){ //Some small cut off
    return;
  }

  //Sort the data
  for(p1=p_in;p1->next!=NULL;p1=p1->next){
    ix=p1->i;
    iy=p1->j;
    if ((ix<0)||(ix>=xmax)||(iy<0)||(iy>=ymax)) continue;
    u=(iy*xmax+ix);
    val1=data[u];
    for(p2=p1->next;p2!=NULL;p2=p2->next){
      ix=p2->i;
      iy=p2->j;
      if ((ix<0)||(ix>=xmax)||(iy<0)||(iy>=ymax)) continue;
      u=(iy*xmax+ix);
      val2=data[u];
      if (val2<val1){
	p2->i=p1->i;
	p2->j=p1->j;
	p1->i=ix;
	p1->j=iy;
	val1=val2;
      }
    }
  }

  (csearch->p)=p_in; //Will only check these points
  (csearch->label_array)=work;
  (csearch->label_value)=label_cur;

  cut_min=1.0e30;
  frac_min=1.0e30;
  //Vary each cut and find max size contiguous above and below cut
  pmin=NULL;
  p1=p_in;
  p2=NULL;
  p=p_in;


  istep=128;
  i=n_total/4;
  while (istep>i)istep=istep/2;
  if (istep==0){
    printf("Very few points in vacuole search: %i\n",n_total);
  }
  while(p!=NULL){
    ix=p->i;
    iy=p->j;
    if ((ix<0)||(ix>=xmax)||(iy<0)||(iy>=ymax)) goto next_p;
    u=(iy*xmax+ix);
    val1=data[u];
    (csearch->cut_low)=val1;
    (csearch->cut_high)=val1;
    (csearch->cut_behavior)=equal_to_labels_cut_above;
    do_contiguous_search(csearch);
    //Find maximum sized list
    nmax=-1;
    for(i=0;i<(csearch->n_lists_found);i++){
      n_contiguous=(csearch->npoints_in_list)[i];
      if (n_contiguous>nmax){
	nmax=n_contiguous;
      }
    }
    nmax1=nmax;

    (csearch->cut_behavior)=equal_to_labels_cut_below;
    do_contiguous_search(csearch);
    //Find maximum sized list
    nmax=-1;
    imax=-1;
    for(i=0;i<(csearch->n_lists_found);i++){
      n_contiguous=(csearch->npoints_in_list)[i];
      if (n_contiguous>nmax){
	imax=i;
	nmax=n_contiguous;
      }
    }
    nmax2=nmax;

    if ((nmax1<0)||(nmax2<0)) goto next_p;
    frac=((float)(nmax1+nmax2))/((float)(n_total));
    if (frac<frac_min){
      frac_min=frac;
      cut_min=(csearch->cut_low);
      pmin=p;
    }

  next_p:
    //printf("---->%i  %i\n",(int)p,istep);fflush(stdout);
    //Step forward by istep--A faster way than to get global minimum
    //than calculating each point
    for(i=0;i<istep;i++){
      p=p->next;
      if (p==p2){
	if (istep==1) goto done;
	if (pmin==NULL){
	  printf("Found nothing in vacuole search.\n");
	  return;
	}
	//We're starting a new loop. Set bounds on p and
	//start again.
	frac_min=1.0e30;
	cut_min=1.0e30;
	j=0;
	for(p1=pmin;(p1->prev)!=NULL;p1=p1->prev){//p1 should never be NULL
	  j++;
	  if (j==istep) break;
	}
	j=0;
	for(p2=pmin;p2!=NULL;p2=p2->next){
	  j++;
	  if (j==istep) break;
	}
	p=p1;
	pmin=NULL;
	istep=istep/2;
	break;
      }
    }

  }

 done:
  //Find point list for best case
  ix=pmin->i;
  iy=pmin->j;
  u=(iy*xmax+ix);
  (csearch->cut_low)=data[u];
  (csearch->cut_high)=data[u];
  (csearch->cut_behavior)=equal_to_labels_cut_below;
  do_contiguous_search(csearch);
  //Find maximum sized list
  nmax=-1;
  imax=-1;
  for(i=0;i<(csearch->n_lists_found);i++){
    n_contiguous=(csearch->npoints_in_list)[i];
    if (n_contiguous>nmax){
      imax=i;
      nmax=n_contiguous;
    }
  }

  offset=(csearch->list_start)[imax];
  for(i=0;i<nmax;i++){
    list_x[i]=*((csearch->list_found_x)+offset+i);
    list_y[i]=*((csearch->list_found_y)+offset+i);
  }

  (*n_points_out)=nmax;
  return;

}

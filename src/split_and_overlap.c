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

#include "tif_routines.h"
#include "split_and_overlap.h"
#include "contiguous.h"

float diff_stat_reg1_reg2(float,float,float,float,float,int *,
			  float *,int,int);

int *find_split_regions(float *fl, int xmax, int ymax){
  //Assume we have a split image. Find a cut value that divides the
  //two images.
  //We're going to return a pointer to an array of integers dimensioned
  //xmax by ymax. 0 means between images, 1 in lower, and 2 in upper.

  int i,j,k;
  int n;
  //int n_large;
  int n_boundary;
  int boundary_low,boundary_high;

  float cut,cut_low,cut_high;
  float cut_step;

  int u;

  int *clist_x;
  int *clist_y;
  int list_cur;

  int *split=NULL;
  int xmax_ymax;

  int n_contiguous;

  float sum;

  struct contiguous_search csearch_memory;
  struct contiguous_search *csearch;

  csearch=&csearch_memory;
  xmax_ymax=xmax*ymax;

  if (split==NULL){
    split=(int *)malloc(xmax_ymax*sizeof(int)); //Will be returned
  }
  memset(split,0,(xmax_ymax)*sizeof(int));

  //To talk to contiguous search routine
  (csearch->p)=NULL; //Will search entire image
  (csearch->data_array)=fl;
  (csearch->xmax)=xmax;
  (csearch->ymax)=ymax;

  //The boundary between upper and lower part of image might not extend
  //entirely across the image. First find these boundary regions.
  //Only use the center region
  //Max and minimum pixels
  boundary_low=(ymax/2-100);
  boundary_high=(ymax/2+100);
  cut_high=0.0;
  cut_low=1.0e30;
  for(i=0;i<xmax;i++){
    for(j=boundary_low;j<boundary_high;j++){
      u=(j*xmax+i);
      split[u]=5;
      if (fl[u]<cut_low){
	cut_low=fl[u];
      }
      if (fl[u]>cut_high){
	cut_high=fl[u];
      }
    }
  }
  (csearch->label_array)=split;
  (csearch->label_value)=5; //Middle region from above

  cut_step=cut_low/50.0;

  (csearch->cut_behavior)=(equal_to_labels_cut_below);
  if (cut_step<=0.0)cut_step=10.0;

  //Now vary the cut until we have a region at least 100 pixels
  n_boundary=100; //Look for a set at least 100 pixels
  cut=cut_low;
  while(cut<cut_high){
    (csearch->cut_low)=cut;
    do_contiguous_search(csearch);
    //printf("%e-->%i",cut,csearch->n_lists_found);fflush(stdout);
    n=0;
    for(i=0;i<(csearch->n_lists_found);i++){
      n_contiguous=(csearch->npoints_in_list)[i];
      if (n_contiguous>=n_boundary)n++;
    }
    //printf("--->%i\n",n);fflush(stdout);
    if(n==1) break;
    cut+=cut_step;
  }

  if (n!=1){ //Couldn't find a good cut
    printf("Couldn't find good cut, dividing in half.\n");
    sum=(float)(ymax/2);
  }else{
    //Average the contiguous region
    sum=0.0;
    clist_y=(csearch->list_found_y);
    for(i=0;i<(csearch->n_lists_found);i++){
      n_contiguous=(csearch->npoints_in_list)[i];
      if (n_contiguous>=n_boundary){
	list_cur=(csearch->list_start)[i];
	for(k=list_cur;k<(list_cur+n_contiguous);k++){
	  sum+=((float)clist_y[k]);
	}
	break;
      }
    }
    sum/=((float)n_contiguous);
    if (sum==0.0){ //Couldn't find region even though was a region
      printf("Middle region was found and then lost.\n");
      error(0); //exit(0); //http://r-pkgs.had.co.nz/src.html
    }
  }

  //Sum is now the average y-value of cut. Extend it across the middle
  //since sometimes the cut region doesn't go all the way across
  k=(int)(sum+0.5);
  for(j=0;j<k;j++){
    u=j*xmax;
    for(i=0;i<xmax;i++){
      split[u]=2;
      u++;
    }
  }
  for(j=k;j<ymax;j++){
    u=j*xmax;
    for(i=0;i<xmax;i++){
      split[u]=1;
      u++;
    }
  }

  //Now search each region separately to try to get rid of more of
  //the boundary region area.
  (csearch->cut_low)=cut; //Same as above
  for(j=1;j<=2;j++){
    (csearch->label_value)=j; //for split[]=j region

    do_contiguous_search(csearch);
    //Remove all pixels found
    clist_x=(csearch->list_found_x);
    clist_y=(csearch->list_found_y);
    for(i=0;i<(csearch->n_lists_found);i++){
      n_contiguous=(csearch->npoints_in_list)[i];
      list_cur=(csearch->list_start)[i];
      for(k=list_cur;k<(list_cur+n_contiguous);k++){
	u=(clist_y[k]*xmax)+clist_x[k];
	split[u]=0;
      }
    }
  }

  return split;
}

/*****************************************************/
void calculate_split_offset(int *split,
			    float *fl,
			    float *mx_out,
			    float *bx_out,
			    float *my_out,
			    float *by_out,
			    int xmax, int ymax){
  //We now have two regions labelled a 1 and 2 in the split[] array.
  //Search for a vector-offset that will overlap them the best.

  //For each value of offset we're going to sum the squares of the
  //differences pixel by pixel. However, first scale the two images
  //since they have different filters, etc.

  int jstart;
  float scale,mean1,mean2;
  float tmp1,tmp2;
  float sum0,sum1;

  float mx,bx,my,by;
  float mx_start,bx_start,my_start,by_start;
  float step_b,step_m;
  float eps_b,eps_m;
  float comp_b,comp_m;
  float scale_b,scale_m;

  int i,j,k;
  int u;
  int xmax_ymax;

  int count_down;

  xmax_ymax=xmax*ymax;
  //Calculate means to get a scale
  tmp1=0.0;
  tmp2=0.0;
  mean1=0.0;
  mean2=0.0;
  for(u=0;u<xmax_ymax;u++){
    if (split[u]==1){
      mean1+=fl[u];
      tmp1+=1.0;
    }else if (split[u]==2){
      mean2+=fl[u];
      tmp2+=1.0;
    }
  }
  mean1/=tmp1;
  mean2/=tmp2;
  scale=mean1/mean2;

  //printf("scale=%e (%e/%e)\n",scale,mean1,mean2);

  //Do search for minimum of "diff_stat_reg1_reg2(offsets)"
  //We're stepping one pixel at a time, so we're probably
  //sensitive to local minima, but probably (0,256) is a good
  //enough first guess to be ok.
  jstart=(ymax/2);
  i=0;
  j=jstart;

  //by=-((float)jstart);
  //my=0.0;
  //bx=0.0;
  //mx=0.0;
  by=*(by_out);
  bx=*(bx_out);
  mx=*(mx_out);
  my=*(my_out);

  scale_b=0.5;
  scale_m=0.5;
  eps_b=scale_b/2.0;
  eps_m=eps_b/((float)xmax);
  comp_b=0.6;
  comp_m=2.5/((float)xmax);
  //printf("bcomp=%e and mcomp=%e\n",comp_b,comp_m);

  sum1=diff_stat_reg1_reg2(mx,bx,my,by,scale,split,fl,xmax,ymax);

  count_down=100;
  do{

    mx_start=mx;
    bx_start=bx;

    my_start=my;
    by_start=by;

    //First we fit for by, then my, then bx, then my.
    //This is because we expect biggest effects in by and my and
    //much less in bx and mx.

    //Search back and forth for by:
    step_b=1.0;
    k=1; //forward or reverse direction
    while(step_b>eps_b){
      do{
	sum0=sum1;
	if (k==0){
	  by+=step_b;
	}else{
	  by-=step_b;
	}
	sum1=diff_stat_reg1_reg2(mx,bx,my,by,scale,split,fl,xmax,ymax);
      }while(sum1<sum0);
      //No longer improving, switch direction, cut step size
      k=1-k;
      step_b*=scale_b;
    }
    //Went 1 too far, undo step (note we just switched direction)
    step_b/=scale_b;
    if (k==0){
      by+=step_b;
    }else{
      by-=step_b;
    }
    sum1=sum0; //Current minimum

    //Search back and forth for my:
    step_m=5.0/((float)ymax);
    k=1; //forward or reverse direction
    while(step_m>eps_m){
      do{
	sum0=sum1;
	if (k==0){
	  my+=step_m;
	}else{
	  my-=step_m;
	}
	sum1=diff_stat_reg1_reg2(mx,bx,my,by,scale,split,fl,xmax,ymax);
      }while(sum1<sum0);
      //No longer improving, switch direction, cut step size
      k=1-k;
      step_m*=scale_m;
    }
    //Went 1 too far, undo step (note we just switched direction)
    step_m/=scale_m;
    if (k==0){
      my+=step_m;
    }else{
      my-=step_m;
    }
    sum1=sum0; //Current minimum

    //Search back and forth for bx:
    step_b=1.0;
    k=0; //forward or reverse direction
    while(step_b>eps_b){
      do{
	sum0=sum1;
	if (k==0){
	  bx+=step_b;
	}else{
	  bx-=step_b;
	}
	sum1=diff_stat_reg1_reg2(mx,bx,my,by,scale,split,fl,xmax,ymax);
      }while(sum1<sum0);
      //No longer improving, switch direction, cut step size
      k=1-k;
      step_b*=scale_b;
    }
    //Went 1 too far (note we just switched directions)
    step_b/=scale_b;
    if (k==0){
      bx+=step_b;
    }else{
      bx-=step_b;
    }
    sum1=sum0; //Current minimum

    //Search back and forth for mx:
    step_m=5.0/((float)xmax);
    k=0; //forward or reverse direction
    while(step_m>eps_m){
      do{
	sum0=sum1;
	if (k==0){
	  mx+=step_m;
	}else{
	  mx-=step_m;
	}
	sum1=diff_stat_reg1_reg2(mx,bx,my,by,scale,split,fl,xmax,ymax);
      }while(sum1<sum0);
      //No longer improving, switch direction, cut step size
      k=1-k;
      step_m*=scale_m;
    }
    //Went 1 too far (note we just switched directions)
    step_m/=scale_m;
   if (k==0){
      mx+=step_m;
    }else{
      mx-=step_m;
    }
    sum1=sum0; //Current minimum

    /*
    printf("----->delbx=%e, delby=%e, delmx=%e delmy=%e\n",
	   ((float)(fabs((double)(bx-bx_start)))),
	   ((float)(fabs((double)(by-by_start)))),
	   ((float)(fabs((double)(mx-mx_start)))),
	   ((float)(fabs((double)(my-my_start))))
	   );
    */

    //printf("offx=(%e)x+(%e), offy=(%e)x+(%e): %e\n",mx,bx,my,by,sum1);

    count_down--;
    if (count_down<=0){
      printf("Too many tries getting slopes\n");
      printf("Going with what we have.\n");
      break;
    }
  }while(
	 (((float)(fabs((double)(bx-bx_start))))>comp_b)||
	 (((float)(fabs((double)(by-by_start))))>comp_b)||
	 (((float)(fabs((double)(mx-mx_start))))>comp_m)||
	 (((float)(fabs((double)(my-my_start))))>comp_m)
	 );

  printf("offx=(%e)x+(%e)\n",mx,bx);
  printf("offy=(%e)x+(%e)\n",my,by);
  (*mx_out)=mx;
  (*bx_out)=bx;
  (*my_out)=my;
  (*by_out)=by;

  return;
}

/*********************************************************/
float diff_stat_reg1_reg2(float mx, float bx,
			  float my, float by,
			  float scale12, int *split,
			  float *fl,
			  int xmax, int ymax){
  //Calculate a statistic that compares regions 1 and 2 as
  //labelled in split[].
  //Offset in x and y direction are calculated as lines with
  //respect to x.
  //The function is, for (x',y') in region 2
  // x' = b_x + m_x*x
  // y' = b_y + m_y*x
  //(Note the slope is respect to x) where (x,y) are in region 1.


  int ui,uj;
  int i,j;
  int u1,u2;
  float sum;
  float total;
  float tmp;

  int offx,offy;

  sum=0.0;
  total=0.0;
  //for(i=0;i<xmax;i++){
  for(i=50;i<xmax-50;i++){
    offx=(int)(bx+mx*((float)i)+0.5);
    offy=(int)(by+my*((float)i)+0.5);
    for(j=50;j<ymax-50;j++){
      u1=j*xmax+i;
      if (split[u1]==1){
	ui=i+offx;
	uj=j+offy;
	if ((ui>=0)&&(ui<xmax)&&(uj>=0)&&(uj<ymax)){
	  u2=uj*xmax+ui;
	  if (split[u2]==2){
	    tmp=((fl[u2]*scale12)-fl[u1]);
	    sum+=(tmp*tmp);
	    total+=1.0;
	  }
	}else{
	  //total=-999.0;
	  //goto finish;
	}
      }
    }
  }

  //finish:
  if (total<=0.0){
    //printf("(%i,%i)=infinite!\n",offx,offy);fflush(stdout);
    return 1.0e30;
  }

  //printf("(%i,%i)=%e\n",offx,offy,(sum/total));fflush(stdout);
  return (sum/total);


}

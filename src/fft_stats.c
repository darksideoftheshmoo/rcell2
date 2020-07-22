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
#include <math.h>
#include <stdlib.h>

#define twopi 3.1415926535*2.0

#include "segment.h"
#include "fft_stats.h"
#include "fft.h"

#define max_log_2 50

int npoints=128;

/**********************************************************************/
double FFT_ratio(struct point *p0){
  //Given a sequence of points making a polygon, calculate FFT of
  //r vs theta (centered at polygon centroid), and return ratio of the
  //energy in all frequencies above 0 to the energy in 0.

  struct complex *ft;

  int k;

  double tmpr,tmpi,f0,f1;
  double return_value;

  if(p0==NULL){
    return -1.0;
  }
  ft=FFT_of_r_vs_theta(p0);
  if(ft!=NULL){
    f0=0.0;
    k=0;
    tmpr=(ft+k)->r;
    tmpi=(ft+k)->i;
    f0 += ( tmpr*tmpr + tmpi*tmpi);
    f1=0.0;
    for(k=1;k<npoints;k+=1){
      tmpr=(ft+k)->r;
      tmpi=(ft+k)->i;
      f1 += (tmpr*tmpr + tmpi*tmpi);
    }
    f0=sqrt(f0);
    f1=sqrt(f1);
    if(f0>0.0){
      tmpr=f1/f0;
    }else{
      tmpr=1.0e30;
    }
    return_value=tmpr;
  }else{
    //    printf("Error in FFT!!!!!!!!!\n");
    return_value=-1.0;
  }

  free(ft);
  return return_value;

}

/**********************************************************************/
struct complex *FFT_of_r_vs_theta(struct point *p){
  //Take a list of (x,y) points as a polygon.  Using the centroid, calculate
  //r vs theta and then calculate FFT.

  struct complex *r=NULL;
  struct complex *ft;
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

  double dtheta,theta;

  int i;

  //Find start of list
  while((p->prev)!=NULL) p=p->prev;
  start=p;

  //Calculate centroid.  The x-component of the centroid is the sum of
  //the average x positions weighted by the length of each segment, and
  //similarly for the y-component.
  p0=start;
  p1=start->next;
  if(p1==NULL){  //Only two points, not a polygon
    return NULL;
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

  r=(struct complex *)malloc(npoints*sizeof(struct complex));
  p0=start;
  p1=start->next;
  //Now transform the polygon into npoints of (r,theta) for theta from 0 to
  //twopi
  dtheta=twopi/((double)npoints);
  for(i=0;i<npoints;i++){
    theta=dtheta*((double)i);
    x1=x0+cos(theta)*1024.0; //1024 to get out of image region
    y1=y0+sin(theta)*1024.0;

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
	//	printf("No intersection, returning NULL\n");
	free(r);
	return NULL;
      }
    }
    //printf("%i %e, %e\n",i,tmpx,tmpy);
    x1=((double)tmpx)-x0;
    y1=((double)tmpy)-y0; 
    (r+i)->r=sqrt( x1*x1 + y1*y1 );
    (r+i)->i=0.0;
  }

  //Calculate FFT
  ft=FFT_1d(r,npoints,1);

  free(r);
  			       
  return ft;
}



















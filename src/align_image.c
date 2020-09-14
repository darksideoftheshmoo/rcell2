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
#include "align_image.h"

float *f_save=NULL;
int xdim=-999;
int ydim=-999;
int xdim_ydim;
int xuse_low=0;
int xuse_high=0;
int yuse_low=0;
int yuse_high=0;
int use_boundary=50;

int offx_prev;
int offy_prev;

double diff_stat(int,int,float*);

/*****************************************************************/
void align_image(float *f, int xdim_in, int ydim_in, int flag){
  //Compare f to f_save to find an overlap offset. Change the values
  //of f to the new offset positions.

  //If flag=0, then just copy f into f_save and do nothing else.
  //If (xdim_in,ydim_in)!=(xdim,ydim) then treat it as if flag=0.

  int i,j,k,l;
  int u1,u2;
  int offx,offy;
  double sum0,sum1;
  float *work;
  int step;

  if ((flag==0)||(xdim_in!=xdim)||(ydim_in!=ydim)){
    free(f_save);
    xdim=xdim_in;
    ydim=ydim_in;
    xdim_ydim=xdim*ydim;

    offx_prev=0;
    offy_prev=0;
    xuse_low=use_boundary;
    if (xuse_low>xdim) xuse_low=0;
    xuse_high=(xdim-use_boundary);
    if (xuse_high<=xuse_low) xuse_high=xdim;
    yuse_low=use_boundary;
    if (yuse_low>ydim) yuse_low=0;
    yuse_high=(ydim-use_boundary);
    if (yuse_high<=yuse_low) yuse_high=ydim;

    f_save=(float *)malloc(xdim_ydim*sizeof(float));
    for (i=0;i<xdim_ydim;i++){
      f_save[i]=f[i];
    }
    return;
  }

  //Vary an offset to minimize diff_stat
  //Calculate an offset from this image to the previous
  //i=offx_prev;
  //j=offy_prev;
  i=0;
  j=0;
  
  for(step=1;step>=1;step/=2){
    do{
      
      offx=i;
      offy=j;

      sum1=diff_stat(i,j,f);
      do{
	      sum0=sum1;
	      i+=step;
	      sum1=diff_stat(i,j,f);
      }while(sum1<sum0);
      i-=step; //Went one too far.
      if (i==offx){ //Didn't move, try other direction
	     sum1=sum0;
	     do{
	       sum0=sum1;
	       i-=step;
	       sum1=diff_stat(i,j,f);
	     }while(sum1<sum0);
	     i+=step; //Went one too far.
      }

      sum1=sum0;
      do{
	     sum0=sum1;
	     j+=step;
	     sum1=diff_stat(i,j,f);
      }while(sum1<sum0);
      j-=step; //Went one too far.

      if (j==offy){ //Didn't move, try other direction
	     sum1=sum0;
	     do{
	       sum0=sum1;
	       j-=step;
	       sum1=diff_stat(i,j,f);
	     }while(sum1<sum0);
	     j+=step; //Went one too far.
      }
      
      
    }while((i!=offx)||(j!=offy)); //Keep going until no changes

  }

  printf("Offset for new fluorescence image=(%i,%i)\n",offx,offy);
  offx_prev=offx;
  offy_prev=offy;
  
  //Correct image for offset
  if ((offx!=0)||(offy!=0)){
    work=(float *)malloc(xdim_ydim*sizeof(float));
    for(i=0;i<xdim_ydim;i++)work[i]=f[i];

    //offx=-offx;
    //offy=-offy;
    for(i=0;i<xdim;i++){
      k=i+offx;
      for(j=0;j<ydim;j++){
	l=j+offy;
	u2=(j*xdim+i);
	if((k>=0)&&(k<xdim)&&(l>=0)&&(l<ydim)){
	  u1=(l*xdim+k); //u1=u2+offset
	  f[u2]=work[u1];
	}else{
	  f[u2]=work[u2]; //Something kind of random on boundaries....
	}
      }
    }
    free(work);
  }

  return;
  
}

/**********************************************************************/
/*
double diff_stat(int offx,int offy, float *f1){
  //Calculate a statistic to look to see if should move images
  //by overall offset. (Comparing fluorescence images f1 and f2
  //for this purpose).
  int i,j,k,l;
  int u1,u2;
  double sum;
  double tmp;
  double total;
  float *f2;

  f2=f_save;

  sum=0.0;
  total=0.0;
  for(i=xuse_low;i<xuse_high;i++){
    k=i+offx;
    for(j=yuse_low;j<yuse_high;j++){
      l=j+offy;
      if((k>=0)&&(k<xdim)&&(l>=0)&&(l<ydim)){
        u1=(l*xdim+k); //u1=u2+offset
	u2=(j*xdim+i);
        tmp=(double)(f2[u2]-f1[u1]);
        sum+=(tmp*tmp);
        total+=1.0;
      }
    }
  }

  sum/=total;
  printf("(%i,%i)=(%e,%e)\n",offx,offy,total,sum);fflush(stdout);
  return sum;
}
*/


/**********************************************************************/
double diff_stat(int offx,int offy, float *f1){
  //Calculate a statistic to look to see if should move images
  //by overall offset. (Comparing fluorescence images f1 and f2
  //for this purpose).
  int i,j,k,l;
  int u1,u2;
  double sum1,sum2;
  double tmp;
  float *f2;

  f2=f_save;

  sum1=0.0;
  sum2=0.0;
  for(i=xuse_low;i<xuse_high;i++){
    k=i+offx;
    for(j=yuse_low;j<yuse_high;j++){
      l=j+offy;
      if((k>=0)&&(k<xdim)&&(l>=0)&&(l<ydim)){
        u1=(l*xdim+k); //u1=u2+offset
	u2=(j*xdim+i);
        tmp=(double)(f2[u2]*f1[u1]);
        sum1+=tmp;
	sum2+=((double)f1[u1]);
      }
    }
  }

  sum1/=sum2;
  //printf("(%i,%i)=(%e,%e)\n",offx,offy,sum2,sum1);fflush(stdout);
  return (-sum1);
}



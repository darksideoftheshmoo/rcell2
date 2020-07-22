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

#include "flatten.h"

float get_slope(float *,int,int,int,int,int,int,int);
void flatten_image_linear(float *,int,int,enum flatten_action_xy);
void flatten_image_piecewise_linear(float *,int,int,
				    enum flatten_action_xy);

/**************************************************************/
float *flatten_image(float *b, int xdim, int ydim,
		     enum flatten_action_overwrite overw,
                     enum flatten_action_xy xy,
		     enum flatten_action_linear lin){

  int i;
  float *cors;
  float av;

  if (overw==return_correction_array){
    //We're going to return a pointer to some allocated memory.
    cors=(float *)malloc(sizeof(float)*xdim*ydim);
    //Start by copying initial image there
    for(i=0;i<(xdim*ydim);i++)cors[i]=b[i];
  }else{ //For overwrite case
    cors=b;
  }
  
  if (lin==linear){
    flatten_image_linear(cors,xdim,ydim,xy);
  }else{
    flatten_image_piecewise_linear(cors,xdim,ydim,xy);
    flatten_image_linear(cors,xdim,ydim,xy);
    flatten_image_piecewise_linear(cors,xdim,ydim,xy);
    flatten_image_linear(cors,xdim,ydim,xy);
    flatten_image_piecewise_linear(cors,xdim,ydim,xy);
    flatten_image_linear(cors,xdim,ydim,xy);  
  }

  if (overw==return_correction_array){
    //Return a scale factor
    //Average correction should be zero by design, but just to make sure:
    av=0.0;
    for(i=0;i<(xdim*ydim);i++){
      av+=(b[i]-cors[i]);
    }
    av/=((float)(xdim*ydim));

    for(i=0;i<(xdim*ydim);i++){
      cors[i]-=av;
      if (b[i]>0.0){
	cors[i]/=b[i];
      }else{
	cors[i]=1.0;
      }
    }
  }
  


  return cors;
}


/**************************************************************/
void flatten_image_linear(float *b, int xdim, int ydim,
			  enum flatten_action_xy xy){
  //Take input image b(xdim,ydim) and flatten based on linear
  //slope.
  //We're going to actually _replace_ b with the flattenned image
  
  //To do this, do the x- and y- projections separately
  
  float mx,my;
  
  int i,j,k;
  float midx,midy;
  
  float cor;
  int count;
  
  printf("Doing linear-correction to flatten image.\n");
  
  //Correct around middle of image
  midx=((float)xdim)/2.0;
  midy=((float)ydim)/2.0;
  
  for (count=0;count<10;count++){
    //x-slope and y-slope
    if (xy!=y_only){
      mx=get_slope(b,xdim,ydim,0,xdim-1,0,ydim-1,0);
    }else{
      mx=0.0;
    }
    if (xy!=x_only){
      my=get_slope(b,xdim,ydim,0,xdim-1,0,ydim-1,1);
    }else{
      my=0.0;
    }
    //printf("Using slope per pixel: mx=%e, my=%e\n",mx,my);
    for(i=0;i<xdim;i++){
      for(j=0;j<ydim;j++){
	k=(j*xdim+i);
	cor=mx*(((float)i)-midx)+my*(((float)j)-midy);
	b[k]-=cor;
      }
    }
    
  }
  
  
  return;
}

/**************************************************************/
void flatten_image_piecewise_linear(float *b, int xdim, int ydim,
				    enum flatten_action_xy xy){

  //Take input image b(xdim,ydim) and flatten based on non-linear
  //slope.
  //We're going to actually _replace_ b with the flattenned image
  //To do the "nonlinear" thing, we're just going to divide the image
  //into different regions and do a separate correction for each.

  //To do this, do the x- and y- projections separately

  float mx,my;

  int i,j,k;
  float midx,midy;

  float cor;
  int count;

  int ndivx=4;
  int ndivy=4;
  int divx,divy;
  int i0,i1,j0,j1;

  int iuse,juse;

  float corav;
  
  printf("Doing piece-wise linear-correction to flatten image.\n");
  
  for (divx=0;divx<ndivx;divx++){
    i0=(xdim/ndivx)*divx;
    i1=i0+(xdim/ndivx)-1;
    for (divy=0;divy<ndivy;divy++){
      j0=(ydim/ndivy)*divy;
      j1=j0+(ydim/ndivy)-1;

      //Correct around middle of image
      midx=((float)i0)+((float)(i1-i0+1))/2.0;
      midy=((float)j0)+((float)(j1-j0+1))/2.0;
      
      for (count=0;count<10;count++){
	//x-slope and y-slope
	if (xy!=y_only){
	  mx=get_slope(b,xdim,ydim,i0,i1,j0,j1,0);
	}else{
	  mx=0.0;
	}
	if (xy!=x_only){
	  my=get_slope(b,xdim,ydim,i0,i1,j0,j1,1);
	}else{
	  my=0.0;
	}
	
	corav=0.0;
	for(i=0;i<xdim;i++){
	  for(j=0;j<ydim;j++){
	    k=(j*xdim+i);

	    if (i<i0){
	      iuse=i0;
	      if (j<j0){
		juse=j0;
	      }else if (j>j1){
		juse=j1;
	      }else{
		juse=j;
	      }
	    }else if (i>i1){
	      iuse=i1;
	      if (j<j0){
		juse=j0;
	      }else if (j>j1){
		juse=j1;
	      }else{
		juse=j;
	      }
	    }else{
	      iuse=i;
	      if (j<j0){
		juse=j0;
	      }else if (j>j1){
		juse=j1;
	      }else{
		juse=j;
	      }
	    }
	    cor=mx*(((float)iuse)-midx)+my*(((float)juse)-midy);

	    b[k]-=cor;
	    corav+=cor;
	  }
	}
	
	//Keep average correction around zero.
	corav/=((float)(xdim*ydim));
	for(i=0;i<(xdim*ydim);i++){
	  b[i]+=corav;
	}
	
      }
      
      //printf("Region x(%i to %i) and y(%i to %i):\n",i0,i1,j0,j1);
      //printf("     Final slope per pixel is mx=%e, my=%e\n",mx,my);
    }
  }
  
  
  return;
}

/**************************************************************/
float get_slope(float *b, int xdim, int ydim,
		int i0, int i1, int j0, int j1,
		int flag){

  //Project onto x- and y- axes and get slopes.
  //(xdim,ydim) are the dimensions of the data array and 
  //(i0,j0)-(i1,j1) mark corners of region we're intersted in.
  //flag=0 means project onto x-axis, 1 means y axis.
  
  int i,j,k,inc;
  double sum;
  int nbins;

  double I2,I,F,FI,N;
  double m;

  if ((i0<0)||(i1>=xdim)||(j0<0)||(j0>=ydim)){
    printf("Unexpected values in get_slope: (%i,%i)-(%i,%i)\n",
	   i0,j0,i1,j1);
    i0=0;
    j0=0;
    i1=xdim;
    j1=ydim;
  }
  
  
  I=0.0;
  I2=0.0;
  F=0.0;
  FI=0.0;
  if (flag==0){ //Project onto x-axis
    nbins=(i1-i0+1);
    for(i=i0;i<=i1;i++){
      sum=0.0;
      k=(j0*xdim+i);
      inc=xdim;
      for(j=j0;j<=j1;j++){
	sum+=((double)b[k]);
	k+=inc;
      }
      I+=((double)i);
      I2+=((double)(i*i));
      F+=sum;
      FI+=((double)i)*sum;
    }
  }else{ //Project onto y-axis
    nbins=(j1-j0+1);
    for(j=j0;j<=j1;j++){
      sum=0.0;
      k=j*xdim+i0;
      inc=1;
      for(i=i0;i<=i1;i++){
	sum+=((double)b[k]);
	k+=inc;
      }
      I+=((double)j);
      I2+=((double)(j*j));
      F+=sum;
      FI+=((double)j)*sum;
    }
  }  
  
  N=(double)nbins;
  
  m=((I*F/N)-FI)/((I2/N)-(I*I));
  
  return ((float)m);
  
}






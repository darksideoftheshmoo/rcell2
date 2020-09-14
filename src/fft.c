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
#include "fft.h"

#define max_log_2 50

struct complex *FFT_recursive(int, int, int);

int current_fft_length=-1;
struct complex *w=NULL;
struct complex *fft_in=NULL;
int inverse_recursive=1;

/**********************************************************************/
struct complex *FFT_2d(struct complex *data, int N, int inverse){
  //Do a full 2d FFT.  This involves an FFT in one direction, followed
  //by an FFT in the other.
  //*data points to an an array of N*N points
  //N is the length of the arrays and the number of them.
  //inverse should be 1 or -1.  If -1, then do the inverse transform.

  //Output points to start of output

  int i,j;
  int step;
  int start;
  struct complex *out;
  struct complex *tmp;

  double dtmp;

  //Make sure number of points N is power of 2
  for(i=N;i>1;i=i/2){
    if((i%2)!=0){
      printf("N not a power of 2: %i.\n",N);
      return NULL;
    }
  }

  //Make sure the exp(twopi*i/N) values are defined correctly
  if(N!=current_fft_length){
    free(w);
    w=(struct complex *)malloc(N*sizeof(struct complex));
    dtmp=twopi/( (double) N );
    for(i=0;i<N;i++){
      (w+i)->r = cos( dtmp*( (double) i ) );
      (w+i)->i = sin( dtmp*( (double) i ) );
    }
    current_fft_length=N;
  }

  out=(struct complex *)malloc(sizeof(struct complex)*N*N);
  inverse_recursive=inverse;

  //Do first direction
  step=1;
  fft_in=data;
  for(i=0;i<N;i++){
    start=i*N;
    tmp=FFT_recursive(N,start,step);
    //Copy onto the output array
    if(tmp==NULL){
      printf("Failure in 2d FFT.\n");
      free(out);
      return NULL;
    }
    for(j=0;j<N;j++){
      (out+start+j)->r = (tmp+j)->r;
      (out+start+j)->i = (tmp+j)->i;
    }
    free(tmp);
  }

  //Now do the other direction on the fft of the first directions
  step=N;
  fft_in=out;
  for(i=0;i<N;i++){
    start=i;
    tmp=FFT_recursive(N,start,step);
    //Copy onto the output array
    if(tmp==NULL){
      printf("Failure in 2d FFT.\n");
      free(out);
      return NULL;
    }
    for(j=0;j<N;j++){
      (out+start+step*j)->r = (tmp+j)->r;
      (out+start+step*j)->i = (tmp+j)->i;
    }
    free(tmp);
  }

  //  if(inverse_recursive==-1){
  //  dtmp=1.0/( (double) (N*N) ); //N*N? or N
  //  for(i=0;i<N*N;i++){
  //    ((out+i)->r)*=dtmp;
  //    ((out+i)->i)*=dtmp;
  //  }
  // }

  dtmp=1.0/( (double) N );
  for(i=0;i<N*N;i++){
    ((out+i)->r)*=dtmp;
    ((out+i)->i)*=dtmp;
  }


  return out;
}

/**********************************************************************/
struct complex *FFT_1d(struct complex *data, int N, int inverse){

  //Return pointer to array of complex fourier coefficients.  Just pass
  //it onto the recursive routine, but do some checks here.
  
  int i;
  int start=0;
  int step=1;

  double dtmp;
  struct complex *ftmp;

  //Make sure number of points N is power of 2
  for(i=N;i>1;i=i/2){
    if((i%2)!=0){
      printf("N not a power of 2: %i.\n",N);
      return NULL;
    }
  }

  //Make sure the exp(twopi*i/N) values are defined correctly
  if(N!=current_fft_length){
    free(w);
    w=(struct complex *)malloc(N*sizeof(struct complex));
    dtmp=twopi/( (double) N );
    for(i=0;i<N;i++){
      (w+i)->r = cos( dtmp*( (double) i ) );
      (w+i)->i = sin( dtmp*( (double) i ) );
    }
    current_fft_length=N;
  }

  //The recursive routine assumes all the data is in the array *in.
  fft_in=data;

  //Define this global variable so recursive routine doesn't have to
  //keep passing it
  inverse_recursive=inverse;

  ftmp=FFT_recursive(N,start,step);
  if(ftmp==NULL){
    printf("Error doing FFT 1d.\n");
    return NULL;
  }

  if(inverse_recursive==-1){
    dtmp=1.0/( (double) N );
    for(i=0;i<N;i++){
      (ftmp+i)->r = ((ftmp+i)->r)*dtmp;
      (ftmp+i)->i = ((ftmp+i)->i)*dtmp;
    }
  }

  return ftmp;


}

/**********************************************************************/
struct complex *FFT_recursive(int N, int start, int step){
  //N is the length of the arrays and must be a power of 2
  //start is location in data array that we consider the start of
  //the array, and step is how to get to the next point.  The data array
  //is the in[] which is defined for this file.

  //inverse is 1 for forward FFT and -1 for reverse.

  //Return NULL if something went wrong, otherwise a pointer to the
  //start of an array of the N fourier coefficients

  //We're going to do this recursively.  There might be faster ways to
  //code this up, but this is the easiest I think.

  struct complex *fft_even;
  struct complex *fft_odd;
  struct complex *fft_N;

  int i,j,k,l;
  int N_over_2;

  double dtmp_even_real,dtmp_even_im;
  double dtmp_odd_real,dtmp_odd_im;
  double w_real,w_im;

  //Define output N-long output array
  fft_N=(struct complex *)malloc(N*sizeof(struct complex));

  //If N=1, then no more recursion, just return the current value
  if(N==1){
    fft_N->r=(fft_in+start)->r;
    fft_N->i=(fft_in+start)->i;
    return fft_N;
  }

  //Get the 2 N/2 FFT transforms of the odd and even components of the
  //current array.  (We're not checking here that N is even.  It should
  //already be checked by the routine calling FFT_recursive (that N is
  //a power of 2).)  Also, we're not going to worry that we're going to
  //step outside of the bounds of in[].  This should also be checked by
  //the calling routine.
  N_over_2=N/2;
  fft_even=FFT_recursive(N_over_2,start,(step*2));
  fft_odd=FFT_recursive(N_over_2,start+step,(step*2));

  l=current_fft_length/N;
  //current_fft_length is the length of the FFT from the original call.
  //It's also the N we used to define the w[] trigonometric terms.  We
  //want here exp(i*twopi/N*k), but w[k] is exp(i*twopi/current_fft_length*k).
  //Since N is current_fft_length divided by 2**x for some integer x, we
  //can just do the above division without worrying about truncation, etc.
  for(i=0;i<N_over_2;i++){
    dtmp_even_real=(fft_even+i)->r;
    dtmp_even_im=(fft_even+i)->i;
    dtmp_odd_real=(fft_odd+i)->r;
    dtmp_odd_im=(fft_odd+i)->i;

    //First do the 0 to N/2 terms
    j=((inverse_recursive)*i*l)%current_fft_length;  //Modulo N division
    if(j<0)j+=current_fft_length;
    w_real=(w+j)->r; //The trigonometric terms
    w_im=(w+j)->i;
    (fft_N+i)->r = dtmp_even_real +
      w_real*dtmp_odd_real - w_im*dtmp_odd_im;
    (fft_N+i)->i = dtmp_even_im +
      w_real*dtmp_odd_im + w_im*dtmp_odd_real;

    //And now the N/2-1 to N-1 terms
    k=i+N_over_2;
    j=((inverse_recursive)*k*l)%current_fft_length;  //Modulo N division
    if(j<0)j+=current_fft_length;
    w_real=(w+j)->r;
    w_im=(w+j)->i;
    (fft_N+k)->r = dtmp_even_real +
      w_real*dtmp_odd_real - w_im*dtmp_odd_im;
    (fft_N+k)->i = dtmp_even_im +
      w_real*dtmp_odd_im + w_im*dtmp_odd_real;

  }
  free(fft_even);  //Keep memory clean!
  free(fft_odd);

  return fft_N;
}


  




















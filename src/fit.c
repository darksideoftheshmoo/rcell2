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
#include "fit.h"

//We're trying to find where the gradient is zero.  Here are the limits
//to define how close to zero we have to be to stop.
#define epsilon_mu .1
#define epsilon_sig .1
#define epsilon_n .1

double gauss_gradients(float *, int, double, double, double,
		       int);

//Fit a gaussian to the data y(i).  We're going to minimize the curve
//s=sum((y_i - gaussian(x_i))^2.  However, we're going to take x_i=i,
//so the fitted mean and sigma should be adjusted for this by the calling
//routine if there's some scale and offset that should be applied.
int gauss_fit(float *y,
	      int n_points,
	      float *mu_fit,
	      float *sig_fit,
	      float *scale_fit){
  //mu_fit, sig_fit, n_fit are starting values for gaussian parameters, and
  //also the returned values.
  //Function returns 1 if all went ok, and 0 otherwise

  //Instead of following the gradient to find the minimimum, I'm going
  //to do vary each variable separately.  This is easier to do, although
  //it's not necessarily more robust.  Also speed isn't really a concern
  //here since we're not going to do this that often.

  double x[3];
  double dx[]={5.0, 1.0, .5};
  double smallest_grad[]={1.0,1.0,1.0};
  double step,grad;
  double s;

  int old_direction, new_direction;

  int i,j;
  int count_down=200; //To force a return

  x[0]=(double)(*mu_fit);
  x[1]=(double)(*sig_fit);
  x[2]=(double)(*scale_fit);

  do{

    s=gauss_gradients(y,n_points,x[0],x[1],x[2],3);
    //    printf("Gauss:(mu,sig,scale)=(%e,%e,%e) (S=%e)\n",
    // 	   x[0],x[1],x[2],s);
    if (s<-1.0e20){
      return 0;
    }
    
    for(i=0;i<3;i++){

      //Initialization for loops below.  Get initial values
      grad=gauss_gradients(y,n_points,x[0],x[1],x[2],i);
      if (grad<-1.0e20){
	return 0;
      }
      if(grad==0.0){
	break;
      }else if(grad>0.0){
	new_direction=1;
	step=-dx[i]; //Stepping against gradient
      }else{
	new_direction=-1;
	step=dx[i];
      }
      
      while(fabs(grad)>smallest_grad[i]){
	old_direction=new_direction;

	//Now step in direction of decreasing s until slope direction
	//switches (we're going to focus in on the slope zero)
	while(new_direction==old_direction){
	  x[i]+=step;

	  //if ((i==0)&&(x[0]<0.0)) { //Consider us done if mean is negative
	  //  grad=0.0;
	  //  goto done_fit;
	  //}
	  grad=gauss_gradients(y,n_points,x[0],x[1],x[2],i);
	  if (grad<-1.0e20){
	    return 0;
	  }	  
	  //printf("%i: grad=%le, step=%le\n",i,grad,step);	
	  if(grad==0.0){
	    break;
	  }else if(grad>0.0){
	    new_direction=1;
	  }else{
	    new_direction=-1;
	  }
	}
	
	//Switch step direction since we crossed zero and scale down so
	//focus in on zero
	step*=(-.5);
      }
      
    }

    if(
       (fabs(gauss_gradients(y,n_points,x[0],x[1],x[2],0)<
	     smallest_grad[0]))&&
       (fabs(gauss_gradients(y,n_points,x[0],x[1],x[2],1)<
	     smallest_grad[1]))&&
       (fabs(gauss_gradients(y,n_points,x[0],x[1],x[2],2)<
	     smallest_grad[2]))
       ){
      j=1;
    }else{
      j=0;
    }
    
    count_down--;
    if (count_down==0){
      printf("Fit taking too long.\n");
      return 0;
    }
    
  }while(j==0);
  
  //	      done_fit:
  s=gauss_gradients(y,n_points,x[0],x[1],x[2],3);
  printf("Final Params: mu=%e, sig=%e, scale=%e (S=%e)\n",
	 x[0],x[1],x[2],s);
  (*mu_fit)=(float)x[0];
  (*sig_fit)=(float) x[1];
  (*scale_fit)=(float)x[2];
  return 1;
}

/************************************************/
double gauss_gradients(float *y, int n_points, 
		       double mu, double sig, double scale,
		       int which){
  //Depending on the value of "which" calculate the gradient with respect
  //to mu, sig, or the normalization factor scale. (This is the gradient
  //s=sum((y_i - gaussian(x_i))^2 where we take x_i=i.)
  
  int i;

  double xi_mu;
  double fi,gi;
  
  double one_over_two_sig_squared;
  double grad;

  if(sig<=0.0){
    printf("Tried to fit Gaussian with non-positive sigma: %e\n",sig);
    return -1.0e30;
  }
  if(scale<=0.0){
    printf("Tried to fit Gaussian with non-positive scale factor: %e\n",
	   scale);
    return -1.0e30;
  }
  
  //Calculate gradients for current parameter settings
  grad=0.0;
    
  //Make a bunch of if checks.  The program would be smaller if we put
  //the if checks in the loops, but it would be slower.
  one_over_two_sig_squared=1.0/2.0/sig/sig;
  if(which==0){ //Gradient with respect to mu
    for(i=0;i<n_points;i++){
      xi_mu=((double)i)-mu;
      fi=scale*exp(-(xi_mu)*(xi_mu)*one_over_two_sig_squared);
      gi=((double)(y[i]))-fi;

      grad+=((xi_mu)*fi*gi);
    }
    grad*=(-2.0/sig/sig);

  }else if(which==1){ //Gradient with respect to sig
    for(i=0;i<n_points;i++){
      xi_mu=((double)i)-mu;
      fi=scale*exp(-(xi_mu)*(xi_mu)*one_over_two_sig_squared);
      gi=((double)(y[i]))-fi;

      grad+=((xi_mu)*(xi_mu)*fi*gi);

    }
    grad*=(-2.0/sig/sig/sig);

  } else if(which==2){ //Gradient with respect to scale
    for(i=0;i<n_points;i++){
      xi_mu=((double)i)-mu;
      fi=scale*exp(-(xi_mu)*(xi_mu)*one_over_two_sig_squared);
      gi=((double)(y[i]))-fi;

      grad+=(fi*gi);

    }
    grad*=(-2.0/scale);

  } else if(which==3){ //Return the value of S
    for(i=0;i<n_points;i++){
      xi_mu=((double)i)-mu;
      fi=scale*exp(-(xi_mu)*(xi_mu)*one_over_two_sig_squared);
      gi=((double)(y[i]))-fi;

      grad+=(gi*gi); //Is really S, not gradient

    }

  } else {
    printf("Didn't choose any of the gradient options in\n");
    printf("gauss_gradients: %i.\n",which);
    return 1.0e30;
  }

  return grad;

}





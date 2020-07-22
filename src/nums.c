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
#include <string.h>
#include "nums.h"
#include "tif_routines.h"

#define xdim 4
#define ydim 6
char zero[]={
  "    "
  " ###"
  " # #"
  " # #"
  " # #"
  " ###"
};
char one[]={
  "    "
  " ## "
  "  # "
  "  # "
  "  # "
  " ###"
};
char two[]={
  "    "
  " ###"
  "   #"
  " ###"
  " #  "
  " ###"
};
char three[]={
  "    "
  " ###"
  "   #"
  " ###"
  "   #"
  " ###"
};
char four[]={
  "    "
  " # #"
  " # #"
  " ###"
  "   #"
  "   #"
};
char five[]={
  "    "
  " ###"
  " #  "
  " ###"
  "   #"
  " ###"
};
char six[]={
  "    "
  " ###"
  " #  "
  " ###"
  " # #"
  " ###"
};
char seven[]={
  "    "
  " ###"
  "   #"
  "  # "
  " #  "
  " #  "
};
char eight[]={
  "    "
  " ###"
  " # #"
  " ###"
  " # #"
  " ###"
};
char nine[]={
  "    "
  " ###"
  " # #"
  " ###"
  "   #"
  " ###"
};

char *numbers[10]={zero,one,two,three,four,five,six,seven,eight,nine};

void add_numbers_to_data(int n, int x, int y, int *d, int xmax, int ymax){

  char *p;

  int u;
  int i,j,k;
  int digit;
  int ix,iy;

  //  y-=3; //So center of number is in middle
  //  x-=2;

  k=n; //Which number we want to use
  do{
    digit=(k%10);
    k=k/10;
    
    p=numbers[digit]; //Which number to add

    for(i=xdim;i>=0;i--){
      for(j=ydim-1;j>=0;j--){
	       u=j*xdim+i;
	       //u=i*ydim+j;
	       if(*(p+u)=='#'){
	         ix=x+i;
	         iy=y+j;
	         if((ix>=0)&&(ix<xmax)&&(iy>=0)&&(iy<ymax)){
	           // if(d[ix][iy]!=found_border){ //Don't overwrite border
	           d[(iy*xmax)+ix]=cell_label;
	           //}
	         }
	       }
      }
    }
    
    x-=xdim; //Move location to left for next digit since we do first 
    //digits of n first
  } while(k>0);


return;
}

/************************************************************/
int search_data_for_number(int n, int *x0, int *y0,
			    int a, int *bfi, int xmax, int ymax){
  //Search the data for the number n. Assume that the numbers are
  //marked with the number (int a). The location is returned as
  //(x0,y0).
  //bf=image of dimensions (xmax,ymax).
  //Return 0 if we found a location, 1 otherwise.
  //This search is pretty slow....

  char *p;

  int u;
  int i,j,k;
  int digit;
  int ix,iy;
  int x,y;
  int xstart,ystart;
  int check_lr;

  (*x0)=0;
  (*y0)=0;

  for(xstart=0;xstart<xmax;xstart++){
    for(ystart=0;ystart<ymax;ystart++){

      x=xstart;
      y=ystart;
      k=n; //Which number we want to use
      do{
	digit=(k%10);
	k=k/10;
	
	p=numbers[digit]; //Which number to look for
 
	for(i=xdim;i>=0;i--){
	  for(j=ydim-1;j>=0;j--){
	    u=j*xdim+i;
	    ix=x+i;
	    iy=y+j;
	    if((ix>=0)&&(ix<xmax)&&(iy>=0)&&(iy<ymax)){
	      if (*(p+u)=='#'){
		if (bfi[(iy*xmax)+ix]!=a){
		  goto next_start;
		}
	      }else{
		if (bfi[(iy*xmax)+ix]==a){
		  goto next_start;
		}
	      }
	    }else{
	      goto next_start; //Numbers don't go over bounds
	    }
	  } 
	} 
	
	x-=xdim; //Move location to left for next digit since we do first 
	// digits of n first
      } while(k>0);

      //printf("Found %i at (%i,%i)\n",n,xstart,ystart);

      check_lr=0;
    check_left_right: ;
      //Before returning check that positions to the left and right
      //don't contain any numbers.
      //We're all primed to check to the left:
      for (digit=0;digit<10;digit++){
	p=numbers[digit]; //Which number to look for
 	for(i=xdim;i>=0;i--){
	  for(j=ydim-1;j>=0;j--){
	    u=j*xdim+i;
	    if (*(p+u)=='#'){
	      ix=x+i;
	      iy=y+j;
	      if((ix>=0)&&(ix<xmax)&&(iy>=0)&&(iy<ymax)){
		if (bfi[(iy*xmax)+ix]!=a)goto next_digit;
	      } 
	    } 
	  } 
	}
	//If we're here then we found a number. Reject this position
	goto next_start;
      next_digit: ;
      }

      if (check_lr==0){
	check_lr=1;
	x=xstart+xdim;
	goto check_left_right;
      }

      //We're done if we're here
      (*x0)=xstart;
      (*y0)=ystart;
      return 0;


    next_start: ;
    }
  }

  //Didn't find number
  return 1;
}


/****************************************************/
void digits_to_string(char *s, int n, int max){
  //Change a number n to a string.  Do it so that all the numbers from
  //0 to max contain the same number of digits.  Ie, if max is 102, then
  //the number 7 would become 007, etc.
  //append the result to s

  int i,j;
  int ten;
  char s2[]=" "; //Initialize like this so it has a '\0' at the end

  *s='\0';
  i=max;
  ten=10;
  while((i/ten)>0){
    ten=ten*10;
  }
  ten=ten/10;

  for(j=ten;j>0;j=j/10){
    s2[0]=(n/j)+'0';
    strcat(s,s2);
    n=(n%j);
  }

return;
}

/*************************************************************/
void date_stamp(int year, int month, int day){
  //just printf() a year/month/day format.  Do it so using
  //digits_to_string() so that say 1 appears as 01 for the month, etc.

  char stamp[10];

  stamp[0]='\0'; //digits_to_string() appends digits to the passed string
  //so be sure to clear it by setting first element to '\0'
  digits_to_string(stamp,year,9999);
  printf("%s/",stamp);
  stamp[0]='\0';
  digits_to_string(stamp,month,12);
  printf("%s/",stamp);
  stamp[0]='\0';
  digits_to_string(stamp,day,31);
  printf("%s",stamp);
  
return;
}

/*************************************************************/
void time_stamp(int hours, int minutes, int seconds){
  //just printf() a hh:mm:ss format.  Do it so using
  //digits_to_string() so that say 1 appears as 01, etc.

  char stamp[10];

  stamp[0]='\0'; //digits_to_string() appends digits to the passed string
  //so be sure to clear it by setting first element to '\0'
  digits_to_string(stamp,hours,24);
  printf("%s:",stamp);
  stamp[0]='\0';
  digits_to_string(stamp,minutes,60);
  printf("%s:",stamp);
  stamp[0]='\0';
  digits_to_string(stamp,seconds,60);
  printf("%s",stamp);
  
return;
}





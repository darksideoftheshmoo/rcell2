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
#include "contiguous.h"

int list_cur;
int *history=NULL;
float *c_array;
int *d_array;

int *clist_x=NULL;
int *clist_y=NULL;
int *list_start=NULL;
int *npoints_in_list=NULL;
int xmax_save=-999;
int ymax_save=-999;
static int xmax;       // made static because of multiple definitions error, see: https://stackoverflow.com/a/7190020/11524079
static int ymax;       // made static
static int xmax_ymax;  // made static

static float cut_low;  // made static
static float cut_high;  // made static
int label_value;

int list_cur;
int n_contiguous;
int cur_history_value=0;
int max_history_value=10000;

void make_contiguous_list(int, int, int (*)(int));

//void make_contiguous_list_cut_low(int,int);
//void make_contiguous_list_cut_high(int,int);
//void make_contiguous_list_cut_between(int,int);
//void make_contiguous_list_label(int,int);
//void make_contiguous_list_label_cut_low(int,int);
//void make_contiguous_list_label_cut_high(int,int);
void update_cur_history_value();

//Simulate a stack (got stack-overflows I think when doing
//recursively)

struct stack_xy {
  unsigned short int i;
  unsigned short int j;
  unsigned char dir;
};
struct stack_xy *stack=NULL;
struct stack_xy *sp;

//The different functions that are used to check in the search
int check_cut_low(int);
int check_cut_high(int);
int check_cut_between(int);
int check_label(int);
int check_label_cut_low(int);
int check_label_cut_high(int);


/*********************************************************/
void do_contiguous_search(struct contiguous_search *csearch){

  int behavior;
  int u;
  int i,j;
  int n_lists;

  int (*check_function)(int); //Pointer to a function (will be which
  //function to use for criterion on which points are to be included
  //in the contiguous list. (The syntax says that "check_function" is
  //a pointer to an int function and takes an int as an argument.)

  struct point *p;

  if (((xmax_save)!=(csearch->xmax))||
      ((ymax_save)!=(csearch->ymax))){
    xmax=(csearch->xmax);
    ymax=(csearch->ymax);
    xmax_save=xmax;
    ymax_save=ymax;
    xmax_ymax=xmax*ymax;
    free(history);
    free(clist_x);
    free(clist_y);
    free(list_start);
    free(npoints_in_list);
    history=NULL;
  }
  if (history==NULL){ //For first call also
    history=(int *)malloc(xmax_ymax*sizeof(int));
    memset(history,cur_history_value,(xmax_ymax)*sizeof(int));
    clist_x=(int *)malloc(xmax_ymax*sizeof(int));
    clist_y=(int *)malloc(xmax_ymax*sizeof(int));
    list_start=(int *)malloc(xmax_ymax*sizeof(int));
    npoints_in_list=(int *)malloc(xmax_ymax*sizeof(int));
  }

  //Always allocate this here (we free it below)
  stack=(struct stack_xy *)malloc((xmax_ymax+1)*sizeof(struct stack_xy));

  c_array=(csearch->data_array);
  d_array=(csearch->label_array);

  list_cur=0;
  n_lists=0;

  update_cur_history_value();

  check_function=NULL;
  //Set up cut-parameters for the contiguous list search
  behavior=(csearch->cut_behavior);
  if (behavior==cut_below){
    cut_low=(csearch->cut_low);
    check_function=check_cut_low;

  }else if (behavior==cut_above){

    cut_high=(csearch->cut_high);
    check_function=check_cut_high;

  }else if (behavior==cut_between){

    cut_low=(csearch->cut_low);
    cut_high=(csearch->cut_high);
    check_function=check_cut_between;

  }else if (behavior==equal_to_labels){

    label_value=(csearch->label_value);
    check_function=check_label;

  }else if (behavior==equal_to_labels_cut_below){

    label_value=(csearch->label_value);
    cut_low=(csearch->cut_low);
    check_function=check_label_cut_low;

  }else if (behavior==equal_to_labels_cut_above){

    label_value=(csearch->label_value);
    cut_high=(csearch->cut_high);
    check_function=check_label_cut_high;

  }else{

    printf("Unknown behavior for contiguous search: %i\n",
	   behavior);fflush(stdout);
    goto done;
  }

  //Call the contiguous-search function for every point depending
  //on whether we want the entire image or the linked list of points
  //pointed to by (csearch->p).
  if ((csearch->p)!=NULL){
    for(p=(csearch->p);p!=NULL;p=(p->next)){
      i=(p->i);
      j=(p->j);
      if ((i>=0)&&(i<xmax)&&(j>=0)&&(j<ymax)){
	      u=(j*xmax+i);
	      if (history[u]!=cur_history_value){
	        list_start[n_lists]=list_cur;
	        make_contiguous_list(i,j,check_function);
	        if (n_contiguous>0){ //Set by make_contiguous_list()
	          npoints_in_list[n_lists]=n_contiguous;
	          n_lists++;
	        }
	      }
      }
    }
  }else{
    for(i=0;i<xmax;i++){
      for(j=0;j<ymax;j++){
	      u=(j*xmax+i);
	      if (history[u]!=cur_history_value){
	        list_start[n_lists]=list_cur;
	        make_contiguous_list(i,j,check_function);
	        if (n_contiguous>0){ //Set by make_contiguous_list()
	          npoints_in_list[n_lists]=n_contiguous;
	          n_lists++;
	        }
	      }
      }
    }
  }

  done:

  (csearch->n_lists_found)=n_lists;
  (csearch->list_start)=list_start;
  (csearch->npoints_in_list)=npoints_in_list;
  (csearch->list_found_x)=clist_x;
  (csearch->list_found_y)=clist_y;

  free(stack);

  return;

}

/*************************************************************/
void update_cur_history_value(){
  cur_history_value++;
  if(cur_history_value>max_history_value){
    cur_history_value=1;
  }
  return;
}

/*************************************************************/
void make_contiguous_list(int istart, int jstart, int (*check)(int)){
  //Do the contiguous search in a non-recursive manner.
  //Start search at (istart,jstart).
  //Last argument says that "check" is a pointer to a function that has
  //one integer argument
  //(I used to do this recursively but on the macintosh it would only do
  //like 5500 recursions before giving a segmentation fault. (On the
  //linux machines it seems to be able to do like 250000 recursions or
  //so and so probably never has that problem. Regardless, I just do a
  //simulated recursion here using stack[] to store the information I
  //need.

  //Return values:
  //clist_x[] and clist[] contain a list of contiguous points, contiguous
  //to (istart,jstart) that pass the cut requirements.
  //At start, list_cur tells where in clist_x[] and clist_y[] to put the
  //next point. On return list_cur points to next free location.
  //n_contiguous is returned as the total number of contiguous points
  //found.

  int u;
  int i,j;
  int dir;

  n_contiguous=0; //Count of how many points we find

  sp=stack; //Start stack pointer at beginning
  sp++;//Note we don't use the sp=stack position, but we use
  //it as a marker below for when we're done
  (sp->i)=(unsigned short int) istart;
  (sp->j)=(unsigned short int) jstart;
  (sp->dir)=0; //All new points start this at zero (will record
  //whether we've checked this point already right, up, down, left, etc)

  do{
    //Check current point
    i=(sp->i);
    j=(sp->j);
    u=(j*xmax+i);

    dir=(int)(sp->dir);
    (sp->dir)++; //Increment for the next time we see this point

    if (dir==0){ //This is a new point, check if it passes

      //Check if we've been here before for this search
      if(history[u]==cur_history_value){
	      sp--; //pop() (so we're done with this point)
	      continue;
      }
      history[u]=cur_history_value; //Save that we've already been here

      //Check if this points passes cut criterion
      if ( ( (*check)(u) )==1 ){ //failed
        sp--; //pop()
      	continue;
      }
      //Good point, add to the list
      clist_x[list_cur]=i;
      clist_y[list_cur]=j;
      list_cur++;
      n_contiguous++;
    }

    //Now determine whether we want to check right, left, up, or down
    if (dir==0){ //Check to the right
      if(i<(xmax-1)){
	       sp++;   //These four lines could constitute a "push()"
	               //routine but I only do them a few (and one
	               //above) and it's probably a little faster to write it
	               //explicitly here.
	       (sp->i)=(unsigned short int )(i+1);
	       (sp->j)=(unsigned short int)(j);
	       (sp->dir)=0; //All new points start this at zero
      }
    }else if (dir==1){ //Check down
      if (j>0){
	       sp++;
	       (sp->i)=(unsigned short int )(i);
	       (sp->j)=(unsigned short int)(j-1);
	       (sp->dir)=0; //All new points start this at zero
      }
    }else if (dir==2){ //Check left
      if (i>0){
	       sp++;
	       (sp->i)=(unsigned short int )(i-1);
	       (sp->j)=(unsigned short int)j;
	       (sp->dir)=0; //All new points start this at zero
      }
    }else if (dir==3){ //Check up
      if (j<(ymax-1)){
	       sp++;
	       (sp->i)=(unsigned short int )(i);
	       (sp->j)=(unsigned short int)(j+1);
	       (sp->dir)=0; //All new points start this at zero
      }
    }else{ //Checked all directions, done with this point
      sp--; //pop()
    }

  }while(sp!=stack); //Keep going until popped all the way back

  return;
}

/*************************************************************/
  int check_cut_low(int u){
  if(c_array[u]>=cut_low){
    return 1;  //If this point fails the criterion

  }
  return 0;
}
/*************************************************************/
int check_cut_high(int u){
  if(c_array[u]<=cut_high){
    return 1;  //If this point fails the criterion
  }
  return 0;
}
/*************************************************************/
int check_cut_between(int u){
  if((c_array[u]>=cut_high)||(c_array[u]<=cut_low)){
    return 1;  //If this point fails the criterion
  }
  return 0;
}
/*************************************************************/
int check_label(int u){
  if (d_array[u]!=label_value){
    return 1;  //If this point fails the criterion
  }
  return 0;
}
/*************************************************************/
int check_label_cut_low(int u){
  if((d_array[u]!=label_value)||(c_array[u]>=cut_low)){
    return 1;  //If this point fails the criterion
  }
  return 0;
}
/*************************************************************/
int check_label_cut_high(int u){
  if((d_array[u]!=label_value)||(c_array[u]<=cut_high)){
    return 1;  //If this point fails the criterion
  }
  return 0;
}


/*
void make_contiguous_list_cut_low(int x, int y){
  //Search through the array c[] for a set of contiguous points that
  //are all less than cut_low.  Use a recursive method.
  //The array history[] is defined externally and is assumed initialized to
  //0 for the first entry.  It tells whether we've checked this point
  //already or not.

  int u;

  u=(y*xmax+x);

  //Check if this point fails the criterion
  if(c_array[u]>=cut_low) return;

  //Check if we've been here before for this search
  if(history[u]==cur_history_value){
    return;
  }
  history[u]=cur_history_value;

  //Add to the list
  clist_x[list_cur]=x;
  clist_y[list_cur]=y;
  list_cur++;
  n_contiguous++;
  //Check to the right
  if(x<xmax-1){
    make_contiguous_list_cut_low(x+1,y);
  }
  //Check down
  if(y>0){
    make_contiguous_list_cut_low(x,y-1);
  }
  //Check to the left
  if(x>0){
    make_contiguous_list_cut_low(x-1,y);
  }
  //Check up
  if(y<ymax-1){
    make_contiguous_list_cut_low(x,y+1);
  }
  return;
}

void make_contiguous_list_cut_high(int x, int y){
  //Search through the array c[] for a set of contiguous points that
  //are all less than cut_low.  Use a recursive method.
  //The array history[] is defined externally and is assumed initialized to
  //0 for the first entry.  It tells whether we've checked this point
  //already or not.

  int u;

  u=(y*xmax+x);

  //Check if this point fails the criterion
  if(c_array[u]<=cut_high) return;

  //Check if we've been here before for this search
  if(history[u]==cur_history_value){
    return;
  }
  history[u]=cur_history_value;

  //Add to the list
  clist_x[list_cur]=x;
  clist_y[list_cur]=y;
  list_cur++;
  n_contiguous++;
  //Check to the right
  if(x<xmax-1){
    make_contiguous_list_cut_high(x+1,y);
  }
  //Check down
  if(y>0){
    make_contiguous_list_cut_high(x,y-1);
  }
  //Check to the left
  if(x>0){
    make_contiguous_list_cut_high(x-1,y);
  }
  //Check up
  if(y<ymax-1){
    make_contiguous_list_cut_high(x,y+1);
  }
  return;
}

void make_contiguous_list_cut_between(int x, int y){
  //Search through the array c[] for a set of contiguous points that
  //are all less than cut_low.  Use a recursive method.
  //The array history[] is defined externally and is assumed initialized to
  //0 for the first entry.  It tells whether we've checked this point
  //already or not.

  int u;

  u=(y*xmax+x);

  //Check if this point fails the criterion
  if (c_array[u]<=cut_low) return;
  if (c_array[u]>=cut_high) return;

  //Check if we've been here before for this search
  if(history[u]==cur_history_value){
    return;
  }
  history[u]=cur_history_value;

  //Add to the list
  clist_x[list_cur]=x;
  clist_y[list_cur]=y;
  list_cur++;
  n_contiguous++;
  //Check to the right
  if(x<xmax-1){
    make_contiguous_list_cut_between(x+1,y);
  }
  //Check down
  if(y>0){
    make_contiguous_list_cut_between(x,y-1);
  }
  //Check to the left
  if(x>0){
    make_contiguous_list_cut_between(x-1,y);
  }
  //Check up
  if(y<ymax-1){
    make_contiguous_list_cut_between(x,y+1);
  }
  return;
}


void make_contiguous_list_label(int x, int y){
  //Search through the array c[] for a set of contiguous points that
  //are all less than cut_low.  Use a recursive method.
  //The array history[] is defined externally and is assumed initialized to
  //0 for the first entry.  It tells whether we've checked this point
  //already or not.

  int u;

  u=(y*xmax+x);

  //Check if this point fails the criterion
  if (d_array[u]!=label_value) return;

  //Check if we've been here before for this search
  if(history[u]==cur_history_value){
    return;
  }
  history[u]=cur_history_value;

  //Add to the list
  clist_x[list_cur]=x;
  clist_y[list_cur]=y;
  list_cur++;
  n_contiguous++;
  //Check to the right
  if(x<xmax-1){
    make_contiguous_list_label(x+1,y);
  }
  //Check down
  if(y>0){
    make_contiguous_list_label(x,y-1);
  }
  //Check to the left
  if(x>0){
    make_contiguous_list_label(x-1,y);
  }
  //Check up
  if(y<ymax-1){
    make_contiguous_list_label(x,y+1);
  }
  return;
}

void make_contiguous_list_label_cut_low(int x, int y){
  //Search through the array c[] for a set of contiguous points that
  //are all less than cut_low.  Use a recursive method.
  //The array history[] is defined externally and is assumed initialized to
  //0 for the first entry.  It tells whether we've checked this point
  //already or not.

  int u;

  u=(y*xmax+x);

  //Check if this point fails the criterion
  if (d_array[u]!=label_value) return;
  if(c_array[u]>=cut_low) return;

  //Check if we've been here before for this search
  if(history[u]==cur_history_value){
    return;
  }
  history[u]=cur_history_value;

  //Add to the list
  clist_x[list_cur]=x;
  clist_y[list_cur]=y;
  list_cur++;
  n_contiguous++;
  //Check to the right
  if(x<xmax-1){
    make_contiguous_list_label_cut_low(x+1,y);
  }
  //Check down
  if(y>0){
    make_contiguous_list_label_cut_low(x,y-1);
  }
  //Check to the left
  if(x>0){
    make_contiguous_list_label_cut_low(x-1,y);
  }
  //Check up
  if(y<ymax-1){
    make_contiguous_list_label_cut_low(x,y+1);
  }
  return;
}

void make_contiguous_list_label_cut_high(int x, int y){
  //Search through the array c[] for a set of contiguous points that
  //are all less than cut_low.  Use a recursive method.
  //The array history[] is defined externally and is assumed initialized to
  //0 for the first entry.  It tells whether we've checked this point
  //already or not.

  int u;

  u=(y*xmax+x);

  //Check if this point fails the criterion
  if (d_array[u]!=label_value) return;
  if(c_array[u]<=cut_high) return;

  //Check if we've been here before for this search
  if(history[u]==cur_history_value){
    return;
  }
  history[u]=cur_history_value;

  //Add to the list
  clist_x[list_cur]=x;
  clist_y[list_cur]=y;
  list_cur++;
  n_contiguous++;
  //Check to the right
  if(x<xmax-1){
    make_contiguous_list_label_cut_high(x+1,y);
  }
  //Check down
  if(y>0){
    make_contiguous_list_label_cut_high(x,y-1);
  }
  //Check to the left
  if(x>0){
    make_contiguous_list_label_cut_high(x-1,y);
  }
  //Check up
  if(y<ymax-1){
    make_contiguous_list_label_cut_high(x,y+1);
  }
  return;
}
*/

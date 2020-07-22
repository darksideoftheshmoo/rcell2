#ifndef _point_
   #define _point_
   #include "point.h"
#endif

//These are the functions for different ways to search
//for contiguous pixels.
//The routine do_contiguous_search() will determine which
//of the different contiguous searches to do.  It's a little
//ugly, but I'm allowing a variable to determine whether to use
//the global c[][] array, the global fl[][] array, etc, etc.
//

enum {cut_below, cut_above, cut_between, equal_to_labels,
      equal_to_labels_cut_below, equal_to_labels_cut_above};

struct contiguous_search {
  struct point *p; //List of points to start at
  float *data_array;
  int xmax;
  int ymax;
  float cut_low;
  float cut_high;
  int cut_behavior;
  int *label_array;
  int label_value;
  int *list_found_x; //Each contiguous set
  int *list_found_y;
  int *list_start;  //start in list_found_x() array
  int *npoints_in_list; //How many points in the list
  int n_lists_found; //How many separate lists
};

void do_contiguous_search(struct contiguous_search *);



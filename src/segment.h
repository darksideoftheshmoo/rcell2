#ifndef _point_
   #define _point_
   #include "point.h"
#endif

int get_data_from_text_file(char *);
int find_cells(struct point ***,struct point ***);
void free_pixels_from_earlier_time_points();
int calculate_fluorescence(void);
int calculate_fluorescence_with_r_info(void);
void add_cell_number_to_the_data(int);
void remove_overlaps(void);
int overlap(struct point *,int,int);
float circumference_squared_over_area(int);
void all_free(void);
void update_list_of_found_cells(int,int,int);
int output_cells(char *,char *,int *);

//V1.2a
int output_cells_single_file(char *,char *,int *);

float integrate_data_along_cell_boundary(struct point *);
int do_segments_intersect(float, float,
                          float, float,
                          float, float,
                          float, float,
                          float *, float *);
void background_level(int);
void internal_structure(int,int);
void next_prev_fl_comparison(void);
void load_global_arrays(int, float *, int *, int, int);
void add_boundary_points_to_data(struct point *, int);
void align_found_cells_to_fl(int);
int recombination_check(int,int,int *,int *,float *);
int  output_individual_cells_to_file(int,char *,
				     float *,
				     int,int,
				     int,int,int);
void debug_test(int);

extern int new_phase;

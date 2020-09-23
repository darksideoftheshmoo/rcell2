
#define cell_in 1 
#define cell_border 2
#define cell_out 3
#define pixel_removed 4
#define found_border 5
#define cell_label 6
#define cell_nucleus 7
//-->Make it possible to have different shades for borders
#define found_border_a 8
#define found_border_b 9
#define found_border_c 10
#define found_border_d 11
#define found_border_e 12
#define found_border_f 13
#define found_border_g 14
#define delete_pixel 15
float *get_data_from_tif_file(char *,int,float *,int *,int *);
int output_data_to_tif_file(char *,float *, int, int,
			    int *,int,int,int);


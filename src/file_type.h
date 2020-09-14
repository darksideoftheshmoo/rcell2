//V1.2a 
#define list_file 0  //List of phase and fuorescence images passed as arguments
                     //default in version 1.1
#define oif_file 1  //Olymus Image File
#define oem_file 2  //Open Enviroment for Microscopy file

extern int file_type;
#ifndef _file_type_
#define _file_type_
int file_type;
#endif


//V1.2a Mapping bright field an fluorescence images
#define bf_fl_mapping_time 0
#define bf_fl_mapping_list 1

extern int bf_fl_mapping;
#ifndef _bf_fl_mapping_
#define _bf_fl_mapping_
int bf_fl_mapping;
#endif


#define bright_field 0
#define metamorph_deconvolution 1
#define hexagonal_grid 2
#define membrane_tagged_fl 3
#define no_third_image 4
#define nuclear_label 5
#define vacuole_label 6
#define fret_bf_bottom_only 7
#define fret_bf_top_only 8
#define fret_bf_bottom_and_top 9
#define unknown_image_type 10
#define confocal_transmission 11

extern int image_type;
extern int third_image_type;
extern int overall_id_offset;
extern int recalculate_internal;
#ifndef _image_type_
#define _image_type_
//int image_type;
//int third_image_type;
//int overall_id_offset;
//int recalculate_internal;
#endif

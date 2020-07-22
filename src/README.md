# Makevars

## Makevars

    CC=ccache clang -Qunused-arguments
    CXX=ccache clang++ -Qunused-arguments
    CCACHE_CPP2=yes
    PKG_LIBS= -ltiff

## Makevars.win

No funciona, probé de todo pero Windows es realmente un bardo.

    CC=ccache clang -Qunused-arguments
    CXX=ccache clang++ -Qunused-arguments
    CCACHE_CPP2=yes
    PKG_CPPFLAGS = -I./inst/libtiff/include
    PKG_LIBS = ../inst/libtiff/lib/x64/libtiff.a ../inst/libtiff/lib/x64/libjpeg.a ../inst/libtiff/lib/x64/libz.a ../inst/libtiff/lib/x64/liblzma.a

## Makevars (old)

    CC=ccache clang -Qunused-arguments
    CXX=ccache clang++ -Qunused-arguments
    CCACHE_CPP2=yes
    PKG_CFLAGS = -I/usr/include/glib-2.0 -I/usr/lib/glib-2.0/include
    PKG_LIBS=-lglib-2.0 -ltiff

# libtiff

Para ver si la puedo incluir dentro del paquete:

  * http://www.libtiff.org/build.html
  * http://mazamascience.com/WorkingWithData/?p=1151
  * http://r-pkgs.had.co.nz/src.html#clang
  * https://stackoverflow.com/a/43599233
  * http://kbroman.org/minimal_make/
  * https://cran.r-project.org/web/packages/roxygen2/vignettes/roxygen2.html
  * https://github.com/r-lib/roxygen2/issues/1063

## para incluirla en windows

https://github.com/ropensci/ijtiff/blob/master/tools/winlibs.R



# CellID original parameters
```
*** Cell_ID Version 1.4.6 *** 

Usage:
  cell [OPTION?] - Cell ID options

Application Options:
  -b, --bright                   Text file list of brightfield tif images
  -f, --fluor                    Text file of fluorescent tif images
  -3, --third                    List of third images to be used
  -d, --dark=dark.txt            List of dark images to be used
  -t, --flat=flat.txt            List of flat images to be used
  -p, --param=parameters.txt     Parameters file 
  -o, --output                   Basename of output files (including dirs)# 

Help Options:
  -h, --help                     Show help options
  --help-all                     Show all help options
  --help-cell                    Parameters used to find cells
  --help-image-type              Specifications for the input images
  --help-misc                    Miscelaneus options

Cell Segmenting
  --brf                          background_reject_factor segmentation parameter
  --max-ppc                      max_pixels_per_cell segmentation parameter
  --min-ppc                      min_pixels_per_cell segmentation parameter
  --max-som                      max_split_over_minor cell splitting parameter
  --max-dow                      max_dist_over_waist cell splitting parameter
  --max-som-t0                   max_split_over_minor_t0 cell splitting parameter
  --max-dow-t0                   max_dist_over_waist_t0 cell splitting parameter
  --track                        tracking_comparison parameter
  --align-fl                     align_fl_to_first/bf option (first, bf)
  --nucl1                        nucleus_radius_1 parameter 
  --nucl2                        nucleus_radius_2 parameter 
  --nucl3                        nucleus_radius_3 parameter 
  --nucl4                        nucleus_radius_4 parameter 
  --nucl5                        nucleus_radius_5 parameter 
  --nucl6                        nucleus_radius_6 parameter 

Image Type options
  --fret-bf                      fret_bf options for split images (top, bottom or both)
  --fret-nuclear                 fret_nuclear options for split images (top, bottom )
  --image-type=bright            image_type options (bright, decon, hex, confocal)
  --3-img-label=none             third_image label option (nuclear, vacuole)

Miscelaneus options
  --force-nuc                    force_nucleus_in_center option
  --list-mapping                 bf_fl_mapping list
  --bf-as-fl                     treat_brightfield_as_fluorescence option
  --align-ind                    align_individual_cells option
  --align-ind-boundary           align_individual_cells_boundary option
  --output-ind                   output_individual_cells option
  --append_output                append_output option, argument is the id-offset (>=0)
  --output_paw                   Make output for PAW.
```
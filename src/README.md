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

# Getting mask data out from cellid

Clues:

## output_data_to_tif_file

En `cell.c` dice:

```C
output_data_to_tif_file(line,  // nombre del archivo ".out.tif"
                        bf,    // numero tipo puntero a float definido en cell.c: float *bf=NULL;
                        xmax,
                        ymax,
                        bf_fl_labels,
                        0,
                        8,
                        0
                        )
```

Está definida en `tif.c` donde dice:

```C
int output_data_to_tif_file(
          char *file,
          float *output_data, // output_data es un array 
          int xmax_data,      //(xmax_data,ymax_data)=size of input array (ie, "output_data")
          int ymax_data,
          int *labels,        // labels es un array: tells where to add boundaries, etc. (NULL for none)
          int type,           // type determines what set of labels to write out
          int bit_size,
          int invert          // if 1: Flip values back from array_max-c[][]
          )
```

## add_boundary_points_to_data

`add_boundary_points_to_data` usage in `cell.c` defined in `segment.c`

Usage passes `NULL` as an argument, and definition says: `//if p_in==NULL then do all n_found borders.`

Takes as input `struct point *p_in`, a struct defined in `point.h` as:

```C
struct point {
  int i,j;
  struct point *next;
  struct point *prev;
};
```

Que parece ser **algo recursivo** : https://stackoverflow.com/questions/588623/self-referential-struct-definition

`int i,j;` quizas vendrian a ser la posición del pixel, pero no hay valor de intensidad ni nada más en la definición...

La cuestion es que next y prev son punteros al elemento siguiente y anterior, respectivamente, de una "lista".

En algun lado del codigo de cellid esta struct se usa para guardar una lista de pares de numeros enteros `i` y `j`.

## Notas sobre este tipo de cosas recursivas en C

En este archivo se define el struct pero no se declara.

Podria hacerse todo a la vez así:

```C
struct point {
  int i,j;
  struct point *next;
  struct point *prev;
} structure_decaration_too;
```

Or only define:

```C
struct {
  int i,j;
  struct point *next;
  struct point *prev;
} structure_decaration_only;
```

### Uso en fft_stats.c

Toma de input un argumento llamado `p0` que es un puntero `*` a un `struct point` (una lista de puntos?).

Del comentario de texto, sabemos que `p0` es una secuencia de puntos que define un polígono.

Como es el `fft.stat`, ese polígono es la frontera o "*boundary*" de la célula.

```C
double FFT_ratio(struct point *p0){
  //Given a sequence of points making a polygon, calculate FFT of
  //r vs theta (centered at polygon centroid), and return ratio of the
  //energy in all frequencies above 0 to the energy in 0.
```

### Uso en contiguous.h

```C
struct contiguous_search {
  struct point *p; //List of points to start at
```

### Ejemplo de uso en `segment.c`

```C
struct point *ptmp;  // definicion de "ptmp" como puntero a la struct "point"

// stuff ...

    for(ptmp=(b0->boundary);ptmp!=NULL;ptmp=(ptmp->next)){
      ix=ptmp->i;
      iy=ptmp->j;
      
      // mode stuff ...
    }
```

### Arrow operator

> (pointer_name)->(variable_name)
>
>  The -> operator in C or C++ gives the value held by "variable_name" to structure variable "pointer_name".

Ver: https://www.geeksforgeeks.org/arrow-operator-in-c-c-with-examples/

Otra explicación:

> foo->bar gets "the member called bar" from "the struct that foo points to".

Entonces:

  * `ptmp=(b0->boundary)` siginifica traeme `boundary` del struct `b0` y asignalo a `ptmp`.
  * en el for loop de `segment.c`, el incrementador es `ptmp=(ptmp->next)` y se interpreta como traeme el `next` de la struct, que en este caso también es un puntero, pero a la siguiente struct de la "lista". 

https://stackoverflow.com/questions/2575048/arrow-operator-usage-in-c

## Conclusion

Hay objetos globales en el programa que guardan una lista de cuales son los puntos que corresponden a un borde o al interior de la célula.

Por ejemplo, `boundary` es un array de estas listas (o sea, un array de `struct point`s).

Lo más probables es que este `boundary` sea un array con un elemento por cada célula en **una** imagen.

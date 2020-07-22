# Including copies 

  cp /usr/include/glib-2.0/glib.h .
  cp -r /usr/include/glib-2.0/glib glib
  cp /usr/lib/glib-2.0/include/glibconfig.h ./

# Makevars

    CC=ccache clang -Qunused-arguments
    CXX=ccache clang++ -Qunused-arguments
    CCACHE_CPP2=yes
    PKG_CFLAGS = -I/usr/include/glib-2.0 -I/usr/lib/glib-2.0/include
    PKG_LIBS=-lglib-2.0 -ltiff

## After glib2.0 include

    CC=ccache clang -Qunused-arguments
    CXX=ccache clang++ -Qunused-arguments
    CCACHE_CPP2=yes
    PKG_CFLAGS = -Iglib_includes
    PKG_LIBS= -ltiff -lcellMagick

Error al cargar paquete con devtools::check():

    Error in dyn.load(dllfile) : 
      unable to load shared object '/home/nicomic/Projects/Rdevel/cellMagick/src/cellMagick.so':
      /home/nicomic/Projects/Rdevel/cellMagick/src/cellMagick.so: undefined symbol: g_option_context_add_main_entries

Clean rebuild:

    /opt/R/3.6.3/bin/R CMD INSTALL --preclean --no-multiarch --with-keep.source cellMagick

No ayuda.

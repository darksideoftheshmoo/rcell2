/*
Function library to handle Oympus Image File (oif) files.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
//#include <tiffio.h>
//#include <time.h>
#include "oif.h"


int oif_load(char *fileName, struct oifData *oif){

  if( (oif->file=fopen(fileName,"r"))!=NULL ){
    printf("Reading %s\n", fileName);
    oif->name=fileName;
    return 1;
  }else {
    printf("File %s not found.\n",fileName);
    return 0;
  }

}
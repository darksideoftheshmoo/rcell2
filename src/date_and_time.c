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
//This program will take an output from /usr/bin/tiffdump and search for
//the special metamorph tag 33628.  The 16th value that tiffdump returns
//is the offset into the tif file of the date and time integers.

//04/16/03-->Changed it to read in the tags from the metamorph file
//and look for tag 33628. If that fails then go to tiffdump and get
//time stamp there. The time stamp is the time the _file_ was made (ie
//if the tif file originated in a stack, it might be the time the
//tif file was taken from the stack, etc.)
//Getting time with 33628 is better since that tag records the actual
//time at which the image was taken.

//(See http://support.universal-imaging.com/docs/T10243.pdf) for
//links to how to read Metamorph tags.


#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "date_and_time.h"

#define UIC1Tag 33628
#define STACK_TAG 33629
#define UIC3Tag 33630
#define UIC4Tag 33631
#define dateTag 306

unsigned long YMDtoJulian(unsigned int, unsigned int, unsigned int);

struct tag_list {
  unsigned short tag;
  unsigned short type;
  unsigned int N; //Number of values for this type
  unsigned int offset; //Contains the value if value fits into 4 bytes
  //Otherwise offset to the data, all offsets from beginning of file
  struct tag_list *next;
  struct tag_list *prev;
};

struct tag_list *read_tags(FILE *);

//Below are used to read characters into bytes and integers,
//defined just before reading in tags
unsigned char *t1=NULL;
unsigned char *t2=NULL;
unsigned char *f1=NULL;
unsigned char *f2=NULL;
unsigned char *f3=NULL;
unsigned char *f4=NULL;
unsigned short two_bytes;
unsigned int four_bytes;


/*********************************************************************/
int get_date_and_time(char *name, int *d, int *t,
				float *xstage_pos, float *ystage_pos){

  FILE *fp;

  int i,j;
  int offset;

  int a1,a2,a3,a4,a5,a6;
  char c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12;
  char c13,c14;
  char zero='0';
  char dummy;

  struct tag_list *tags;
  struct tag_list *tg;

  int offset16,offset28;
  int xnum,xden,ynum,yden;

  int good=1;

  (*xstage_pos)=-1.0;
  (*ystage_pos)=-1.0;
  (*d)=0;
  (*t)=0;

  //Open tif file
  if((fp=fopen(name,"rb"))==NULL){
    printf("Couldn't open %s file.\n",name);
    return 0;
  }

  //Get list of all the tags (sorted numerically)
  tags=read_tags(fp);
  if (tags==NULL){
    printf("Failed to read in TIFF tags.\n");
    fclose(fp);
    return 0;
  }
  //Make sure tags is at start of list (it should be, but check anyway)
  while((tags->prev)!=NULL)tags=(tags->prev);

  //Look for tag 33628
  for(tg=tags;(tg!=NULL);tg=tg->next){
    if((tg->tag)==UIC1Tag) break;
  }
  if (tg==NULL) goto get_standard_time_stamp;

  //We found the Metamorph 33268 tag.
  offset=tg->offset;
  fseek(fp,offset,SEEK_SET); //Sets file pointer to offset bytes  

  //The values consist of an ID and a value. We want to get the
  //value corresponding to ID=16
  offset16=-1;
  offset28=-1;
  for(i=0;i<(tg->N);i++){
    *f1=fgetc(fp);
    *f2=fgetc(fp);
    *f3=fgetc(fp);
    *f4=fgetc(fp);
    j=four_bytes;
    
    *f1=fgetc(fp);
    *f2=fgetc(fp);
    *f3=fgetc(fp);
    *f4=fgetc(fp);
    offset=four_bytes;
    
    if (j==16){ //ID=16 is creation time
      offset16=offset;
    }else if (j==28){
      offset28=offset; //ID=28 is stage position
    }
    
  }
  
  if (offset28>0){ //For stage position in absolute coordinates
    //Value points to an offset that contains a table of long integers
    //where values are Xposition-numerator, Xposition-denominator,
    //Yposition-numerator, Yposition-denominator, for each plane in the
    //image.
    offset=offset28;
    fseek(fp,offset,SEEK_SET);
    *f1=fgetc(fp);
    *f2=fgetc(fp);
    *f3=fgetc(fp);
    *f4=fgetc(fp);
    xnum=four_bytes;
    *f1=fgetc(fp);
    *f2=fgetc(fp);
    *f3=fgetc(fp);
    *f4=fgetc(fp);
    xden=four_bytes;
    *f1=fgetc(fp);
    *f2=fgetc(fp);
    *f3=fgetc(fp);
    *f4=fgetc(fp);
    ynum=four_bytes;
    *f1=fgetc(fp);
    *f2=fgetc(fp);
    *f3=fgetc(fp);
    *f4=fgetc(fp);
    yden=four_bytes;

    if ((xden>0)&&(yden>0)){
      (*xstage_pos)=((float)xnum)/((float)xden);
      (*ystage_pos)=((float)ynum)/((float)yden);
    }
  }
  
  
  
  if (offset16<0){
    printf("17th value in 33628 tag should have been 16. It was %i.\n",j);
    fclose(fp);
    goto get_standard_time_stamp;
  }
  
  offset=offset16;
  //The value points to the offset that contains the time and date
  fseek(fp,offset,SEEK_SET);
  *f1=fgetc(fp);
  *f2=fgetc(fp);
  *f3=fgetc(fp);
  *f4=fgetc(fp);
  *d=four_bytes; //"Julian date (days since 1/1/4713 B.C.)"
  //see http://www.universal-imaging.com/ftp/support/stack/STK.doc
  //which describes metamorph tiff tag number 33628
  //They have the comment: "For more information about this method of numbering
  //days, the best explanation we've found is “Peter Meyer's Julian Day Numbers” 
  //(http://www.magnet.ch/serendipity/hermetic/cal_stud/jdn.htm)."
  //(This is no longer there but a Google search for
  //"Peter Meyer's Julian Day Numbers" gives good results)

  *f1=fgetc(fp);
  *f2=fgetc(fp);
  *f3=fgetc(fp);
  *f4=fgetc(fp);
  *t=four_bytes; //time
  //My comparisons to what the metamorph program itself gives predict that
  //t is the time in milliseconds since 12:00 midnight, ie, the number of
  //milliseconds in that day.  And date is an integer that counts the
  //number of days.  To convert to day/month/year, pass date to the
  //JulianToYMD(unsigned long, unsigned short *, unsigned char *,
  //unsigned char *)
  //routine below.
  
  goto close_everything;
  
 get_standard_time_stamp:
  //Didn't find metamorph encoded date, or not using
  //Look for the datetime stamp instead  
  
  //Look for tag dateTag
  while((tags->prev)!=NULL)tags=(tags->prev);
  for(tg=tags;(tg!=NULL);tg=tg->next){
    if((tg->tag)==dateTag) break;
  }
  if (tg==NULL){
    printf("No time stamp found\n");
    good=0;
    goto close_everything;
  }
  offset=tg->offset;
  fseek(fp,offset,SEEK_SET); //Sets file pointer to offset bytes  
  
  //Format is "YYYY:MM:DD HH:MM:SS" for creation time of file
  fscanf(fp,
	 "%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c%c",
	 &c1,&c2,&c3,&c4,&dummy,&c5,&c6,&dummy,
	 &c7,&c8,&dummy,&c9,&c10,&dummy,&c11,&c12,&dummy,&c13,&c14);
  a1=(c1-zero)*1000+(c2-zero)*100+(c3-zero)*10+(c4-zero);
  a2=(c5-zero)*10+(c6-zero)+1; //Does months as 0 to 11, but
  //it seems Metamorph does 1-12 and we want to keep compatible
  a3=(c7-zero)*10+(c8-zero);
  a4=(c9-zero)*10+(c10-zero);
  a5=(c11-zero)*10+(c12-zero);
  a6=(c13-zero)*10+(c14-zero);
  //The below sscan() won't work because it reads in say
  //"04" as 0 and then 4, instead of just 4.
  //sscanf(s,"%i:%i:%i %i:%i:%i",&a1,&a2,&a3,&a4,&a5,&a6);
  //printf("%i %i %i %i %i %i\n",a1,a2,a3,a4,a5,a6);
  //fflush(stdout);

  //Convert this to the Julian thing that Metamorph does. This is good
  //because it keeps things in same format for the rest of the code and
  //also, I'm sorting based on the time stamps, so keeping it like this
  //avoids having to check changes in the year, etc (maybe user is running
  //on New Year's eve....)
  (*d)=YMDtoJulian(a1,a2,a3);
  *t=((a4*3600+a5*60+a6)*1000); //milliseconds since midnight
  
 close_everything:
  
  fclose(fp);
  //Free tags linked list
  if ((tags->next)==NULL){
    free(tags);
  }else{
    for(tg=(tags->next);(tg->next)!=NULL;tg=tg->next){
      free(tg->prev);
    }
    free(tg);
  }

  return good; //Done
}




/****************************************************/
void JulianToYMD(unsigned int julian,
                 unsigned short *year,
		 unsigned char *month,
                 unsigned char *day){
  //This program is from http://www.universal-imaging.com/ftp/support/stack/STK.doc
  //which describes the secret metamorph tiff tag number 33628

  int a, b, c, d, e, alpha;
  int z = julian+1;

  // dealing with Gregorian calendar reform

  if (z < 2299161L) {
           a = z;

  } else {

           alpha = (long) ((z - 1867216.25) / 36524.25);

	          a = z + 1 + alpha - alpha / 4;

  }
  b = ( a > 1721423L ? a + 1524 : a + 1158 );

  c = (long) ((b - 122.1) / 365.25);

  d = (long) (365.25 * c);

  e = (long) ((b - d) / 30.6001);

  *day = (unsigned char)(b - d - (long)(30.6001 * e));

  *month = (unsigned char)((e < 13.5) ? e - 1 : e - 13);

  *year = (short)((*month > 2.5 ) ? (c - 4716) : c - 4715);

return;
}

/****************************************************/
unsigned long YMDtoJulian(unsigned int year,
			  unsigned int month,
			  unsigned int day){
  short a,b = 0;
  short work_year = year;
  short work_month = month;
  short work_day = day;

  // correct for negative year
  if (work_year < 0) {
    work_year++;
  }

  if (work_month <= 2) {
    work_year--;
    work_month += (short)12;
  }
  
  // deal with Gregorian calendar
  if (work_year*10000. + work_month*100. + work_day >= 15821015.){
    a = (short)(work_year/100.);
    b = (short)(2 - a + a/4);
  }
  
  return(long) (365.25*work_year) +
    (long) (30.6001 * (work_month+1)) +
    work_day + 1720994L + b;
}


/*********************************************************************/
void Int_To_Hours_Minutes_Seconds(unsigned int t, 
				  int *hours,
				  int *minutes,
				  int *seconds){
  //My own routine to take the total integer t which is seconds from midnight
  //and convert to h:m:s. (Note it's _seconds_, not milliseconds)
  //--asg
  
  double tmp;

  tmp=((double)t)/3600.0; //3600 is number of seconds per hour
  *hours=(int)tmp;
  tmp-=floor(tmp); //Fractional part of an hour
  tmp*=60.0; //Number of minutes
  *minutes=(int)tmp;
  tmp-=floor(tmp);
  tmp*=60.0; //Number of seconds
  *seconds=(int)tmp;

  return;
}


/*****************************************************************/
struct tag_list *read_tags(FILE *fp){
  //Read in all the tags from the file pointed at by fp.

  unsigned char c1,c2;
  unsigned short n_entries;
  unsigned int IFD;

  int error;
  int i;

  struct tag_list *tags=NULL;
  struct tag_list *cur;

  int little_big_endian=0;
  int flip_endian=0;

  //First check whether the cpu we're on is big endian or little
  //endian (which is different than the file itself, and will tell
  //us how to read interpret the characters we read)
  t1=(unsigned char *)&two_bytes;
  t2=(t1+1);
  (*t1)=42;(*t2)=0;
  if (two_bytes==42){
    little_big_endian=0; //little endian
  }else{
    little_big_endian=1; //big endian
  }
  flip_endian=0; //Assume we match the tif-file type so no flipping

  //Read in header information
  c1=fgetc(fp); //First two bytes tell whether big- or little-endian
  c2=fgetc(fp);
  if ((c1=='I')&&(c2=='I')){ //File is little endian
    if (little_big_endian==1){ //but cpu is big endian
      flip_endian=1;
    }
  }else if ((c1!='I')||(c2!='I')){
    if ((c1=='M')&&(c2=='M')){
      //File itself is big endian
      if (little_big_endian==0){ //but cpu is little
	flip_endian=1; //Flip stuff
      }
    }
  }else{
    printf("File is not a tiff file.\n");
    return NULL;
  }

  if (flip_endian==0){
    t1=(unsigned char *)&two_bytes;
    t2=(t1+1);
    f1=(unsigned char *)&four_bytes; //Next two bytes the number 42 for tiff
    f2=(f1+1);
    f3=(f2+1);
    f4=(f3+1);
  }else{
    t2=(unsigned char *)&two_bytes;
    t1=(t2+1);
    f4=(unsigned char *)&four_bytes; //Next two bytes the number 42 for tiff
    f3=(f4+1);
    f2=(f3+1);
    f1=(f2+1);
  }   

  (*t1)=fgetc(fp);(*t2)=fgetc(fp);
  if (two_bytes!=42){ //next two bytes 42 for tiff files
    printf("File is not a tiff file: %i\n",two_bytes);
    return NULL;
  }

  //Start the list
  tags=(struct tag_list *)malloc(sizeof(struct tag_list));
  tags->next=NULL;
  tags->prev=NULL;
  tags->tag=0;
  tags->type=0;
  tags->N=0;
  tags->offset=0;

  //Next four bytes offset to first IFD;
  (*f1)=fgetc(fp);(*f2)=fgetc(fp);(*f3)=fgetc(fp);(*f4)=fgetc(fp);
  IFD=four_bytes;
  while(IFD!=0){

    if ((error=fseek(fp,IFD,SEEK_SET))!=0){
      //SEEK_SET means offset is relative to start of file
      //SEEK_CUR would mean relative to current position
      //SEEK_END would mean relative to the end of the file
      printf("Error finding next IFD: %i\n",error);
      return NULL;
    }

    //First two bytes are number of entries
    (*t1)=fgetc(fp);(*t2)=fgetc(fp);
    n_entries=two_bytes;
    for (i=0;i<n_entries;i++){
      //Make a new tag structure
      cur=(struct tag_list *)malloc(sizeof(struct tag_list));
      (*t1)=fgetc(fp);(*t2)=fgetc(fp);
      (cur->tag)=two_bytes;

      (*t1)=fgetc(fp);(*t2)=fgetc(fp);
      (cur->type)=two_bytes;

      (*f1)=fgetc(fp);(*f2)=fgetc(fp);(*f3)=fgetc(fp);(*f4)=fgetc(fp);
      (cur->N)=four_bytes;

      (*f1)=fgetc(fp);(*f2)=fgetc(fp);(*f3)=fgetc(fp);(*f4)=fgetc(fp);
      (cur->offset)=four_bytes;

      //      printf("Tag %i, field_type=%i, count=%i, offset=%x\n",
      //             cur->tag,cur->type,cur->N,cur->offset);fflush(stdout);

      //Place in structure sorted according to tag number.
      while((tags->tag)<(cur->tag)){
        if((tags->next)==NULL)break;
        tags=(tags->next);
      }
      while((tags->tag)>(cur->tag)){
        if((tags->prev)==NULL)break;
        tags=(tags->prev);
      }
      //Want to place ours right after the current location.
      (cur->next)=(tags->next);
      if ((tags->next)!=NULL){
        (tags->next)->prev=cur;
      }
      (cur->prev)=tags;
      tags->next=cur;

    }//End loop over number of entries in this IFD

    //Offset to next IFD if there is one
    (*f1)=fgetc(fp);(*f2)=fgetc(fp);(*f3)=fgetc(fp);(*f4)=fgetc(fp);
    IFD=four_bytes;

  }

  //Finally, remove the first element (which is just a marker) from
  //the list and return the start of the list.
  while ((tags->prev)!=NULL) tags=(tags->prev);

  if (tags!=NULL){
    tags=(tags->next);
    free(tags->prev);
    tags->prev=NULL;
  }

  return tags;
}



float get_exposure(char *name){
  //Look for the expression "Exposure: " in the binary file. Metamorph
  //seems to record the exposure time like that.


  FILE *fp;
  char exposure[]="Exposure: ";
  unsigned char c;
  int fc;
  int i;
  int last;
  char line[500];
  float time;

  //Open tif file
  if((fp=fopen(name,"rb"))==NULL){
    printf("Couldn't open %s file.\n",name);
    return -999.0;
  }
  
  last=strlen(exposure);

  time=-999.0;
  i=0;
  while ((fc=fgetc(fp))!=EOF){
    c=(unsigned char)fc;
    if (c!=exposure[i]){
      i=0;
    }else{
      i=i+1;
    }
    if (i==last){
      //Found word "Exposure"--Read characters until we get to a space
      i=0;
      do {
	fc=fgetc(fp);
	c=(unsigned char)fc;
	line[i]=c;
	i++;
      } while (c!=' ');
      line[i-1]='\0'; //Replace final ' ' with '\0'
      sscanf(line,"%e",&time);
      break;
    }
  }
  fclose(fp);
  return time;
}

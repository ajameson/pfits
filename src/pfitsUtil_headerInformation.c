#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include "fitsio.h"
#include <stdint.h>

#define VERSION 1.0

int main(int argc,char *argv[])
{
  fitsfile *infptr;
  int status=0;
  int intVal=0;
  int colnum=0;
  float floatVal=0;
  char str[1024];
  
  if ( !fits_open_file(&infptr, argv[1], READONLY, &status) )
    {
      fits_movabs_hdu(infptr, 1, NULL, &status);
      fits_read_key(infptr,TSTRING,"OBSERVER",str,NULL,&status); printf("[OBSERVER] %s\n",str);     
      fits_read_key(infptr,TSTRING,"PROJID",str,NULL,&status); printf("[PROJID] %s\n",str);
      fits_read_key(infptr,TSTRING,"FRONTEND",str,NULL,&status); printf("[FRONTEND] %s\n",str);
      fits_read_key(infptr,TSTRING,"BACKEND",str,NULL,&status); printf("[BACKEND] %s\n",str);
      fits_read_key(infptr,TSTRING,"BECONFIG",str,NULL,&status); printf("[BECONFIG] %s\n",str);
      fits_read_key(infptr,TSTRING,"OBS_MODE",str,NULL,&status); printf("[OBS_MODE] %s\n",str);
      fits_read_key(infptr,TSTRING,"DATE-OBS",str,NULL,&status); printf("[DATE-OBS] %s\n",str);
      fits_read_key(infptr,TSTRING,"SRC_NAME",str,NULL,&status); printf("[SRC_NAME] %s\n",str);
      fits_read_key(infptr,TSTRING,"RA",str,NULL,&status); printf("[RA] %s\n",str);
      fits_read_key(infptr,TSTRING,"DEC",str,NULL,&status); printf("[DEC] %s\n",str);     
      printf("\n");
      fits_movnam_hdu(infptr,BINARY_TBL,(char *)"SUBINT",0,&status);
      fits_read_key(infptr,TINT,"NAXIS2",&intVal,NULL,&status); printf("[SUBINT:NSUB] %d\n",intVal);
      fits_read_key(infptr,TINT,"NCHAN",&intVal,NULL,&status); printf("[SUBINT:NCHAN] %d\n",intVal);
      fits_read_key(infptr,TSTRING,"CHAN_BW",str,NULL,&status); printf("[SUBINT:CHAN_BW] %s\n",str);     
      fits_read_key(infptr,TINT,"NPOL",&intVal,NULL,&status); printf("[SUBINT:NPOL] %d\n",intVal);
      fits_read_key(infptr,TSTRING,"POL_TYPE",str,NULL,&status); printf("[SUBINT:POL_TYPE] %s\n",str);     
      fits_read_key(infptr,TINT,"NBIN",&intVal,NULL,&status); printf("[SUBINT:NBIN] %d\n",intVal);
      fits_read_key(infptr,TSTRING,"TBIN",str,NULL,&status); printf("[SUBINT:TBIN] %s\n",str);
      fits_read_key(infptr,TSTRING,"NBITS",str,NULL,&status); printf("[SUBINT:NBITS] %s\n",str);
      fits_read_key(infptr,TSTRING,"NSBLK",str,NULL,&status); printf("[SUBINT:NSBLK] %s\n",str);

      if (status) fits_report_error(stderr, status);
      fits_close_file(infptr, &status);
    }

  

}

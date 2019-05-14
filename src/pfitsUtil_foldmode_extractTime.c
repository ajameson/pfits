// Code to extract one or more sub-integrations from a specified PSRFITS file

// Currently copy the entire file and then delete
// Could also just copy the subintegrations of interest *** TO DO

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include "fitsio.h"
#include <stdint.h>

#define VERSION 1.0

int main(int argc,char *argv[])
{
  fitsfile *infptr,*outfptr;
  int status=0;
  int i,j,k;
  int hdu=1;
  char inname[1024];
  char outname[1024];
  char usename[1024];
  int nchan,npol,nbin,nsubint;
  int outSub1 = 20;
  int outSub2 = 25;
  int newNchan;
  char nullVal = 0;
  float nullVal_f = 0;
  int initflag=0;
  
  int colnum_data_in=0;
  int colnum_data_out=0;
  int colnum_out_datFreq=0;
  int colnum_in_datFreq=0;
  int colnum_out_datWts=0;
  int colnum_in_datWts=0;
  int colnum_out_datScl=0;
  int colnum_in_datScl=0;
  int colnum_out_datOffs=0;
  int colnum_in_datOffs=0;
  
  int16_t *dataIn;
  int16_t *dataOut;
  float *colArrayIn;
  float *colArrayOut;
  
  int writePos;
  int readPos;

  long naxes[4];
  int naxis=3;
  char tdim[16];

  printf("Running pfitsUtil_foldmode_extractTime\n");
  
  for (i=1;i<argc;i++)
    {
      if (strcmp(argv[i],"-f")==0) // input file
	strcpy(inname,argv[++i]);
      else if (strcmp(argv[i],"-o")==0) // output file
	strcpy(outname,argv[++i]);
    }

  sprintf(usename,"!%s",outname);
  printf("Using files: %s and %s\n",inname,outname);
  
  // Copy the first file to the output file
  /* Open the input file */

  if ( !fits_open_file(&infptr, inname, READONLY, &status) )
    {
      /* Create the output file */
      if ( !fits_create_file(&outfptr, usename, &status) )
	{
	  /* Copy every HDU until we get an error */
	  while( !fits_movabs_hdu(infptr, hdu++, NULL, &status) )
	    fits_copy_hdu(infptr, outfptr, 0, &status);
	  /* Reset status after normal error */
	  if (status == END_OF_FILE) status = 0;	  
	}
      fits_movnam_hdu(infptr,BINARY_TBL,(char *)"SUBINT",0,&status);
      fits_read_key(infptr,TINT,"NCHAN",&nchan,NULL,&status);
      fits_read_key(infptr,TINT,"NPOL",&npol,NULL,&status);
      fits_read_key(infptr,TINT,"NBIN",&nbin,NULL,&status);
      fits_read_key(infptr,TINT,"NAXIS2",&nsubint,NULL,&status);
      
    }
    fits_close_file(infptr, &status);	      
    /* if error occured, print out error message */
    if (status) fits_report_error(stderr, status);

    fits_movnam_hdu(outfptr,BINARY_TBL,(char *)"SUBINT",0,&status);
    if (outSub2 < nsubint)
      {
	printf("Deleting from subint #%d with %d subints (%d) (%d)\n",outSub2+1,nsubint-outSub2-1,outSub2+1+nsubint-outSub2-1,nsubint);
	fits_delete_rows(outfptr,outSub2+1,nsubint-outSub2-1,&status);
      }

    if (outSub1 > 1)
      fits_delete_rows(outfptr,1,outSub1-1,&status);

    fits_close_file(outfptr, &status);
    if (status) fits_report_error(stderr, status);
    printf("Completed\n");
}

// Code to extract one or more frequency bands from a specified PSRFITS file

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
  int outChan1 = 600;
  int outChan2 = 700;
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

  printf("Running pfitsUtil_foldmode_extractFreq\n");
  
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
    
    /* if error occured, print out error message */
    if (status) fits_report_error(stderr, status);

    newNchan= outChan2-outChan1+1; // SHOULD CHECK FOR BAND INVERSION
  
    dataIn = (int16_t *)malloc(sizeof(int16_t)*nchan*nbin*npol);
    dataOut = (int16_t *)malloc(sizeof(int16_t)*newNchan*nbin*npol);
    colArrayIn = (float *)malloc(sizeof(float)*nchan*npol);
    colArrayOut = (float *)malloc(sizeof(float)*newNchan*npol);

    printf("Re-sizing file size\n");
    
    fits_movnam_hdu(infptr,BINARY_TBL,(char *)"SUBINT",0,&status);
    fits_movnam_hdu(outfptr,BINARY_TBL,(char *)"SUBINT",0,&status);
    fits_get_colnum(infptr, CASEINSEN, "DATA", &colnum_data_in, &status);
    fits_get_colnum(outfptr, CASEINSEN, "DATA", &colnum_data_out, &status);  
    if (status) {fits_report_error(stderr, status); exit(1);}

    fits_update_key(outfptr, TINT, (char *)"NCHAN", &newNchan, NULL, &status );
    fits_modify_vector_len(outfptr,colnum_data_out,(long)((long)newNchan*(long)npol*(long)nbin),&status);

    fits_get_colnum(outfptr, CASEINSEN, "DAT_FREQ", &colnum_out_datFreq, &status);
    fits_get_colnum(infptr, CASEINSEN, "DAT_FREQ", &colnum_in_datFreq, &status);  
    fits_modify_vector_len(outfptr,colnum_out_datFreq,newNchan,&status);
    fits_get_colnum(outfptr, CASEINSEN, "DAT_WTS", &colnum_out_datWts, &status);
    fits_get_colnum(infptr, CASEINSEN, "DAT_WTS", &colnum_in_datWts, &status);  
    fits_modify_vector_len(outfptr,colnum_out_datWts,newNchan,&status);
    fits_get_colnum(outfptr, CASEINSEN, "DAT_SCL", &colnum_out_datScl, &status);
    fits_get_colnum(infptr, CASEINSEN, "DAT_SCL", &colnum_in_datScl, &status);  
    fits_modify_vector_len(outfptr,colnum_out_datScl,newNchan*npol,&status);
    fits_get_colnum(outfptr, CASEINSEN, "DAT_OFFS", &colnum_out_datOffs, &status);
    fits_get_colnum(infptr, CASEINSEN, "DAT_OFFS", &colnum_in_datOffs, &status);  
    fits_modify_vector_len(outfptr,colnum_out_datOffs,newNchan*npol,&status);
    if (status) {fits_report_error(stderr, status); exit(1);}
    
    naxes[0] = nbin;
    naxes[1] = newNchan;
    naxes[2] = npol;

    sprintf(tdim,"TDIM%d",colnum_data_out);
    fits_delete_key(outfptr, tdim, &status);
    fits_write_tdim(outfptr, colnum_data_out, naxis, naxes, &status);
    if (status) {fits_report_error(stderr, status); exit(1);}
    printf("Complete re-sizeing file\n");
    
    for (i=0;i<nsubint;i++)
      {
	printf("Processing sub-integration %d\n",i+1);
	fits_read_col(infptr,TSHORT,colnum_data_in,i+1,1,nbin*nchan*npol,&nullVal,dataIn,&initflag,&status);
	for (j=0;j<npol;j++)
	  {
	    writePos = j*newNchan*nbin;
	    readPos  = j*nchan*nbin+outChan1*nbin; // SHOULD CHECK FOR BAND INVERSION
	    memcpy(dataOut+writePos,dataIn+readPos,newNchan*nbin*2); // *2 for TSHORT
	  }
	fits_write_col(outfptr,TSHORT,colnum_data_out,i+1,1,newNchan*nbin*npol,dataOut,&status);
	if (status) fits_report_error(stderr, status);

	fits_read_col(infptr,TFLOAT,colnum_in_datFreq,i+1,1,nchan,&nullVal_f,colArrayIn,&initflag,&status);
	fits_write_col(outfptr,TFLOAT,colnum_out_datFreq,i+1,1,newNchan,colArrayIn+outChan1,&status);
	fits_read_col(infptr,TFLOAT,colnum_in_datWts,i+1,1,nchan,&nullVal_f,colArrayIn,&initflag,&status);
	fits_write_col(outfptr,TFLOAT,colnum_out_datWts,i+1,1,newNchan,colArrayIn+outChan1,&status);

	fits_read_col(infptr,TFLOAT,colnum_in_datScl,i+1,1,nchan*npol,&nullVal_f,colArrayIn,&initflag,&status);
	for (j=0;j<npol;j++)
	  {
	    for (k=0;k<newNchan;k++)
	      colArrayOut[j*newNchan+k] = colArrayIn[j*nchan+k];	    
	  }
	fits_write_col(outfptr,TFLOAT,colnum_out_datScl,i+1,1,newNchan*npol,colArrayOut,&status);

	fits_read_col(infptr,TFLOAT,colnum_in_datOffs,i+1,1,nchan*npol,&nullVal_f,colArrayIn,&initflag,&status);
	for (j=0;j<npol;j++)
	  {
	    for (k=0;k<newNchan;k++)
	      colArrayOut[j*newNchan+k] = colArrayIn[j*nchan+k];	    
	  }
	fits_write_col(outfptr,TFLOAT,colnum_out_datOffs,i+1,1,newNchan*npol,colArrayOut,&status);

      }

    fits_close_file(infptr, &status);	  
    fits_close_file(outfptr, &status);
    free(dataIn);
    free(dataOut);
    free(colArrayIn);
    free(colArrayOut);
    printf("Completed\n");
}

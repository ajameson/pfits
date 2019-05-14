// Code to combine two or more PSRFITS files in the frequency direction to
// produce a single output file.
//
// Usage:
// ./pfitsUtil_foldmode_combineFreq -o <output name> file1 file2 file3 file4 file5
//
// This software will work with PSRFITS fold mode files and PSRFITS search mode files
//
// Version 1: 3rd May 2018, G. Hobbs
//
// gcc -lm -o pfitsUtil_foldmode_combineFreq pfitsUtil_foldmode_combineFreq.c -lcfitsio
//

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
  char outname[1024];
  char usename[1024];
  char inname[26][1024];
  int i,j,k,p;
  int nIn=0;
  int hdu=1;
  int colnum=0;
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
  
  int nchan,nbin,npol;
  long naxes[4];
  int naxis=3;
  int newNchan;
  int16_t *dataVals;
  int16_t nullVal = 0;
  float nullVal_f = 0;
  float *dataArray;
  int subint_in;
  int subint_out;
  int initflag=0;
  int bytespersample = 2;
  unsigned long writePos;
  char tdim[16];
  
  for (i=1;i<argc;i++)
    {
      if (strcmp(argv[i],"-o")==0)
	strcpy(outname,argv[++i]);
      else
	{
	  strcpy(inname[nIn],argv[i]);
	  nIn++;
	  printf("nIn = %d\n",nIn);
	}
    }
  sprintf(usename,"!%s",outname);
  
  printf("fitsFreqCombine version %d\n",(int)(VERSION));

  printf("Output file is %s\n",outname);
  
  // Copy the first file to the output file
  /* Open the input file */
  if ( !fits_open_file(&infptr, inname[0], READONLY, &status) )
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

      fits_close_file(infptr, &status);
    }
    
    /* if error occured, print out error message */
    if (status) fits_report_error(stderr, status);

    // Fix some header parameters
    fits_movabs_hdu(outfptr, 1, NULL, &status);
    fits_update_key(outfptr, TSTRING, (char *)"BACKEND", (char *)"Medusa", NULL, &status );
    fits_update_key(outfptr, TSTRING, (char *)"FRONTEND", (char *)"UWL", NULL, &status );
    //    fits_update_key(outfptr, TSTRING, (char *)"SRC_NAME", (char *)"J1022+1001", NULL, &status );
    if (status) {fits_report_error(stderr, status); exit(1);}
    printf("Completed fixing header parameters\n");
    

    // Allocate memory
    dataVals = (int16_t *)malloc(sizeof(int16_t)*nchan*nbin);
    dataArray = (float *)malloc(sizeof(float)*nchan);
    
    newNchan = nchan*nIn;

    fits_movnam_hdu(outfptr,BINARY_TBL,(char *)"SUBINT",0,&status);
    fits_get_colnum(outfptr, CASEINSEN, "DATA", &colnum_data_out, &status);  
    fits_modify_vector_len(outfptr,colnum_data_out,(nchan*npol*nbin*nIn),&status);

    fits_update_key(outfptr, TINT, (char *)"NCHAN", &newNchan, NULL, &status );
	
    naxes[0] = nbin;
    naxes[1] = nchan*nIn;
    naxes[2] = npol;
    sprintf(tdim,"TDIM%d",colnum_data_out);
    fits_delete_key(outfptr, tdim, &status);
    fits_write_tdim(outfptr, colnum_data_out, naxis, naxes, &status);

    if (status) {fits_report_error(stderr, status); exit(1);}
    printf("Completed updating output file size\n");

    
    // Now increase sizes for DAT_FREQ, DAT_WTS, DAT_SCL, DAT_OFFS
    fits_get_colnum(outfptr, CASEINSEN, "DAT_FREQ", &colnum_out_datFreq, &status);  
    fits_modify_vector_len(outfptr,colnum_out_datFreq,newNchan,&status);
    fits_get_colnum(outfptr, CASEINSEN, "DAT_WTS", &colnum_out_datWts, &status);  
    fits_modify_vector_len(outfptr,colnum_out_datWts,newNchan,&status);
    fits_get_colnum(outfptr, CASEINSEN, "DAT_SCL", &colnum_out_datScl, &status);  
    fits_modify_vector_len(outfptr,colnum_out_datScl,newNchan*npol,&status);
    fits_get_colnum(outfptr, CASEINSEN, "DAT_OFFS", &colnum_out_datOffs, &status);  
    fits_modify_vector_len(outfptr,colnum_out_datOffs,newNchan*npol,&status);
    if (status) {fits_report_error(stderr, status); exit(1);}
    printf("Reading input files\n");
    
    
    // Copy DAT_FREQ from the other files directly to ensure no confusion if frequency bandwidths etc. have changed
    for (i=0;i<nIn;i++)
      {
	fits_open_file(&infptr, inname[i], READONLY, &status);
	if (status) fits_report_error(stderr, status);
	printf("Processing input file %d %d\n",i+1,colnum_data_out);
	
	subint_in = 1; // Should check if there is only 1 subint
	subint_out = 1; // Should check this
	
	fits_movnam_hdu(infptr,BINARY_TBL,(char *)"SUBINT",0,&status);
	fits_get_colnum(infptr, CASEINSEN, "DATA", &colnum_data_in, &status);  
	fits_get_colnum(infptr, CASEINSEN, "DAT_FREQ", &colnum_in_datFreq, &status);  
	fits_get_colnum(infptr, CASEINSEN, "DAT_WTS", &colnum_in_datWts, &status);  
	fits_get_colnum(infptr, CASEINSEN, "DAT_SCL", &colnum_in_datScl, &status);  
	fits_get_colnum(infptr, CASEINSEN, "DAT_OFFS", &colnum_in_datOffs, &status);  

	if (status) fits_report_error(stderr, status);
	
	for (p=0;p<npol;p++)
	  {
	    fits_read_col(infptr,TSHORT,colnum_data_in,subint_in,1+p*nchan*nbin,nbin*nchan,&nullVal,dataVals,&initflag,&status);
	    if (status) {printf("POS 0\n"); fits_report_error(stderr, status); exit(1);}
	    writePos = 1+p*newNchan*nbin + i*nchan*nbin;
	    fits_write_col(outfptr,TSHORT,colnum_data_out,subint_out,writePos,nchan*nbin,dataVals,&status);
	    if (status) {printf("POS 1\n"); fits_report_error(stderr, status); exit(1);}

	    // Copy the DAT_FREQ, DAT_OFFS, DAT_SCL and DAT_WTS columns
	    if (p==0) // DAT_FREQ and DAT_WTS just for channels not polarisations
	      {
		writePos = 1+i*nchan;

		printf("Writing col\n");
		fits_read_col(infptr,TFLOAT,colnum_in_datFreq,subint_in,1,nchan,&nullVal_f,dataArray,&initflag,&status);
		fits_write_col(outfptr,TFLOAT,colnum_out_datFreq,subint_out,writePos,nchan,dataArray,&status);
		//
		fits_read_col(infptr,TFLOAT,colnum_in_datWts,subint_in,1,nchan,&nullVal_f,dataArray,&initflag,&status);
		fits_write_col(outfptr,TFLOAT,colnum_out_datWts,subint_out,writePos,nchan,dataArray,&status);
		if (status) {printf("POS 3\n"); fits_report_error(stderr, status); exit(1);}
	      }
	    //
	    // DAT_OFFS and DAT_SCL also change with polarisation
	    writePos = 1+p*newNchan+i*nchan;
	    
	    fits_read_col(infptr,TFLOAT,colnum_in_datScl,subint_in,1+p*nchan,nchan,&nullVal_f,dataArray,&initflag,&status);
	    fits_write_col(outfptr,TFLOAT,colnum_out_datScl,subint_out,writePos,nchan,dataArray,&status);
	    //
	    fits_read_col(infptr,TFLOAT,colnum_in_datOffs,subint_in,1+p*nchan,nchan,&nullVal_f,dataArray,&initflag,&status);
	    fits_write_col(outfptr,TFLOAT,colnum_out_datOffs,subint_out,writePos,nchan,dataArray,&status);
	    
	    // 
	    printf("Next pol\n");
	  }
	
	fits_close_file(infptr, &status);	
      }


    // Deallocate memory
    free(dataVals);
    free(dataArray);
    
    printf("Closing the output file\n");
    // Now close the file
    fits_close_file(outfptr,  &status);
    if (status) fits_report_error(stderr, status);
    
}

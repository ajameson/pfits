// Code to combine two or more PSRFITS search mode files in the frequency direction to
// produce a single output file.
//
// Usage:
// ./pfitsUtil_searchmode_combineFreq -o <output name> file1 file2 file3 file4 file5
//
// This software will work with PSRFITS fold mode files and PSRFITS search mode files
//
// Version 1: 3rd May 2018, G. Hobbs
//
// gcc -lm -o pfitsUtil_searchmode_combineFreq pfitsUtil_searchmode_combineFreq.c -lcfitsio
// gcc -lm -o pfitsUtil_searchmode_combineFreq pfitsUtil_searchmode_combineFreq.c -I/pulsar/psr/software/20170525/src/util/anaconda2/include -L../cfitsio/ -L/pulsar/psr/software/20170525/src/util/anaconda2/lib/ -lcfitsio -lcurl -lssl -lcrypto

#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include "fitsio.h"
#include <stdint.h>

#define VERSION 1.0
#define MAX_SUBBANDS 30

int main(int argc,char *argv[])
{
  fitsfile *infptr,*outfptr;
  int status=0;
  char outname[1024];
  char usename[1024];
  char inname[MAX_SUBBANDS][1024];
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
  
  int nchan,nsblk,npol,nbit;
  long naxes[4];
  int naxis=3;
  int newNchan;
  char *dataVals;
  char *writeVals;
  char nullVal = 0;
  float nullVal_f = 0;
  float *dataArray;
  int subint_in;
  int subint_out;
  int initflag=0;
  float bytespersample;
  float newObsFreq;
  long long sizeWriteVals = 0;
  long long writePos;
  char cval[16];
  int nsubint=0;
  char tdim[16];
  int data_colnum;
  float freqFirst,freqLast;
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
  
  printf("pfitsUtil_searchmode_cmobineFreq version %d\n",(int)(VERSION));
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
      fits_read_key(infptr,TINT,"NSBLK",&nsblk,NULL,&status);
      fits_read_key(infptr,TINT,"NBITS",&nbit,NULL,&status);
      printf("Nchan = %d, npol = %d nsblk = %d nbit = %d\n",nchan,npol,nsblk,nbit);
      fits_read_key(infptr,TINT,"NAXIS2",&nsubint,NULL,&status);
      fits_close_file(infptr, &status);
    }

  // SET NSUBINT TO 1
  //  nsubint = 2;
  // REMOVE THAT

  
  bytespersample = nbit/8.0;
  printf("Created initial file\n");
  /* if error occured, print out error message */
  if (status) fits_report_error(stderr, status);

  //  fits_close_file(outfptr,  &status);
  //  if (status) fits_report_error(stderr, status);


  newNchan = nchan*nIn;

  // Fix some header parameters
  fits_movabs_hdu(outfptr, 1, NULL, &status);
  if (status) {fits_report_error(stderr, status); exit(1);}
  fits_update_key(outfptr, TINT, (char *)"OBSNCHAN", &newNchan, NULL, &status );

  if (status) {fits_report_error(stderr, status); exit(1);}

  printf("Completed fixing header parameters\n");

  // Increase the number of channels
  
  // Allocate memory
  dataVals = (char *)malloc(sizeof(char)*nchan*nsblk*bytespersample);
  sizeWriteVals = (long long)(((double)newNchan*(double)nsblk*(double)nsubint*(double)npol)*bytespersample);
  writeVals = (char *)malloc(sizeof(char)*sizeWriteVals);
  dataArray = (float *)malloc(sizeof(float)*nchan); 
  printf("Size of writeVals = %Ld\n",sizeWriteVals);

  fits_movnam_hdu(outfptr,BINARY_TBL,(char *)"SUBINT",0,&status);
  fits_get_colnum(outfptr, CASEINSEN, "DATA", &colnum_data_out, &status);  
  if (status) {fits_report_error(stderr, status); exit(1);}
  
  printf("Starting increase %d %d %d %d %d %d\n",nchan,npol,nsblk,nIn,colnum_data_out,nchan*npol*nsblk*nIn);
  fits_report_error(stderr,status);
 
  fits_modify_vector_len(outfptr,colnum_data_out,(long)((long)nchan*(long)npol*(long)nsblk*nIn*bytespersample),&status);
  if (status) {fits_report_error(stderr, status); exit(1);}
  //  printf("Increasing NCHAN\n");
  fits_update_key(outfptr, TINT, (char *)"NCHAN", &newNchan, NULL, &status );
  if (status) {fits_report_error(stderr, status); exit(1);}
  printf("Increased vector len %g\n",bytespersample);
  sprintf(tdim,"TDIM%d",colnum_data_out);

  //  fits_close_file(outfptr,  &status);
  //  exit(1);
  
    //
  naxes[0] = (int)(nchan*nIn); 
  naxes[1] = npol;
  naxes[2] = nsblk*bytespersample;

  fits_delete_key(outfptr, tdim, &status);
  fits_write_tdim(outfptr, colnum_data_out, naxis, naxes, &status);
  if (status) {fits_report_error(stderr, status); exit(1);}
  printf("Completed updating output file size\n");
  
  // Now increase sizes for DAT_FREQ, DAT_WTS, DAT_SCL, DAT_OFFS
  fits_get_colnum(outfptr, CASEINSEN, "DAT_FREQ", &colnum_out_datFreq, &status);  
  if (status) {fits_report_error(stderr, status); exit(1);}
  printf("Colnum = %d %d\n",colnum_out_datFreq,newNchan);
  fits_modify_vector_len(outfptr,colnum_out_datFreq,newNchan*npol,&status);
  if (status) {fits_report_error(stderr, status); exit(1);}
  printf("GOT HERE\n");
  
  fits_get_colnum(outfptr, CASEINSEN, "DAT_WTS", &colnum_out_datWts, &status);  
  fits_modify_vector_len(outfptr,colnum_out_datWts,newNchan*npol,&status);
  if (status) {fits_report_error(stderr, status); exit(1);}
  
  fits_get_colnum(outfptr, CASEINSEN, "DAT_SCL", &colnum_out_datScl, &status);  
  fits_modify_vector_len(outfptr,colnum_out_datScl,newNchan*npol,&status);
  if (status) {fits_report_error(stderr, status); exit(1);}
  
  fits_get_colnum(outfptr, CASEINSEN, "DAT_OFFS", &colnum_out_datOffs, &status);  
  fits_modify_vector_len(outfptr,colnum_out_datOffs,newNchan*npol,&status);
  
  
  if (status) {fits_report_error(stderr, status); exit(1);}
  printf("Reading input files\n");

  for (i=0;i<nIn;i++)
    {
      printf("Opening %s\n",inname[i]);
      fits_open_file(&infptr, inname[i], READONLY, &status);
      if (status) fits_report_error(stderr, status);
      fits_movnam_hdu(infptr, BINARY_TBL,(char *)"SUBINT",0,&status);
      fits_get_colnum(infptr, CASEINSEN, "DATA", &colnum_data_in, &status);  
      fits_get_colnum(infptr, CASEINSEN, "DAT_FREQ", &colnum_in_datFreq, &status);  
      fits_get_colnum(infptr, CASEINSEN, "DAT_WTS", &colnum_in_datWts, &status);  
      fits_get_colnum(infptr, CASEINSEN, "DAT_SCL", &colnum_in_datScl, &status);  
      fits_get_colnum(infptr, CASEINSEN, "DAT_OFFS", &colnum_in_datOffs, &status);  
      if (status) fits_report_error(stderr, status);
      
      for (k=0;k<nsubint;k++)
	{
	  printf("Reading subint %d\n",k);
	  writePos = 1+i*nchan;
	  
	  fits_read_col(infptr,TFLOAT,colnum_in_datFreq,k+1,1,nchan,&nullVal_f,dataArray,&initflag,&status);
	  if (i==0)
	    freqFirst = dataArray[0];
	  if (i==nIn-1)
	    freqLast = dataArray[nchan-1];
	  fits_write_col(outfptr,TFLOAT,colnum_out_datFreq,k+1,writePos,nchan,dataArray,&status);
	  //
	  fits_read_col(infptr,TFLOAT,colnum_in_datWts,k+1,1,nchan,&nullVal_f,dataArray,&initflag,&status);
	  fits_write_col(outfptr,TFLOAT,colnum_out_datWts,k+1,writePos,nchan,dataArray,&status);
	  if (status) {fits_report_error(stderr, status); exit(1);}

	  for (p=0;p<npol;p++)
	    {
	      //
	      // DAT_OFFS and DAT_SCL also change with polarisation
	      writePos = 1+p*newNchan+i*nchan;
	      fits_read_col(infptr,TFLOAT,colnum_in_datScl,k+1,1+p*nchan,nchan,&nullVal_f,dataArray,&initflag,&status);
	      fits_write_col(outfptr,TFLOAT,colnum_out_datScl,k+1,writePos,nchan,dataArray,&status);
	      //

	      fits_read_col(infptr,TFLOAT,colnum_in_datOffs,k+1,1+p*nchan,nchan,&nullVal_f,dataArray,&initflag,&status);
	      fits_write_col(outfptr,TFLOAT,colnum_out_datOffs,k+1,writePos,nchan,dataArray,&status);
	      fits_read_col(infptr,TBYTE,colnum_data_in,k+1,1+(p*nchan*nsblk)*bytespersample,nchan*bytespersample*nsblk,&nullVal,dataVals,&initflag,&status);
	      for (j=0;j<nsblk;j++)
		{
		  writePos = (long long)(((double)k*npol*newNchan*nsblk + (double)p*newNchan*nsblk + (double)j*newNchan + (double)i*nchan)*bytespersample);
		  if (writePos < 0 || writePos >= sizeWriteVals) {
		    printf("ERROR: writePos = %Ld (sizeWriteVals = %Ld)\n",writePos,sizeWriteVals);
		    exit(1);
		  }
		  memcpy(writeVals+writePos,dataVals+(int)(j*nchan*bytespersample),nchan*bytespersample);
		}
	    }
	}
      fits_close_file(infptr, &status);	      
    }
  printf("Writing the data\n");
  for (i=0;i<nsubint;i++)
    fits_write_col(outfptr,TBYTE,colnum_data_out,i+1,1,npol*newNchan*nsblk*bytespersample,writeVals+(long long)((double)i*(double)newNchan*(double)nsblk*(double)npol*(double)bytespersample),&status);

  // Update more header parameters
  fits_movabs_hdu(outfptr, 1, NULL, &status);
  printf("Frequency range = %g and %g\n",freqLast,freqFirst);
  newObsFreq = fabs((freqLast + freqFirst)/2.0);
  fits_update_key(outfptr, TFLOAT, (char *)"OBSFREQ", &newObsFreq, NULL, &status );

  fits_movnam_hdu(outfptr, BINARY_TBL,(char *)"SUBINT",0,&status);
  fits_update_key(outfptr, TFLOAT, (char *)"REFFREQ", &newObsFreq, NULL, &status );

  
    // Deallocate memory
    free(dataVals);
    free(dataArray);
    free(writeVals);
    printf("Closing the output file\n");
    // Now close the file
    fits_close_file(outfptr,  &status);
    if (status) fits_report_error(stderr, status);
    
}

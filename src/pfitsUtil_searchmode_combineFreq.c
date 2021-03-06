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
#include <unistd.h>

#define VERSION 1.0
#define MAX_SUBBANDS 30

void usage(void);

void usage()
{
  fprintf(stdout, "pfitsUtil_searchmode_combineFreq [opts] inputs\n");
  fprintf(stdout, "  input       input files, at least 2 required\n");
  fprintf(stdout, "  -f          force combination\n");
  fprintf(stdout, "  -h          display usage\n");
  fprintf(stdout, "  -o output   output filename\n");
}

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
  double nullVal_d = 0;
  float *dataArray;
  double *dataArrayFP64;
  int subint_in;
  int subint_out;
  int initflag=0;
  float bytespersample;
  double newObsFreq;
  long long sizeWriteVals = 0;
  long long writePos;
  char cval[16];
  int nsubints[MAX_SUBBANDS];
  char tdim[16];
  int data_colnum;
  double freqFirst = 0;
  double freqLast = 0;
  double obsBw = 0;
  double newObsBw = 0;
  int force_merge = 0;

  int arg;
  while ((arg = getopt(argc, argv, "fho:")) != -1)
  {
    switch (arg)
    {
      case 'f':
        force_merge = 1;
        break;

      case 'h':
        usage ();
        return 0;

      case 'o':
        strcpy(outname, optarg);
        break;

      default:
        fprintf(stderr, "ERROR: Unrecognised command line option\n");
        usage();
        return EXIT_FAILURE;
    }
  }

  int num_input_files = argc - optind;
  if (num_input_files < 2)
  {
    fprintf(stderr, "ERROR: 2 command line arguments are required\n");
    usage();
    exit(EXIT_FAILURE);
  }

  if (strlen(outname) == 0)
  {
    fprintf(stderr, "ERROR: output filename must be specified\n");
    usage();
    return EXIT_FAILURE;
  }

  for (i=0; i<num_input_files; i++)
  {
    strcpy(inname[nIn],argv[optind + i]);
    nIn++;
  }
  sprintf(usename,"!%s",outname);
  
  printf("pfitsUtil_searchmode_combineFreq version %d\n",(int)(VERSION));
  
  // in case the nsubint are different - use the smallest one
  for (i=0; i<nIn; i++)
  {
    if ( !fits_open_file(&infptr, inname[i], READONLY, &status) )
    {
      fits_movabs_hdu(infptr, 1, NULL, &status);
      fits_read_key(infptr, TDOUBLE, (char *) "OBSBW", &obsBw, NULL, &status);
      newObsBw += obsBw;

      fits_movnam_hdu(infptr,BINARY_TBL,(char *)"SUBINT",0,&status);
      fits_read_key(infptr,TINT,"NAXIS2",&nsubints[i],NULL,&status);
      fits_close_file(infptr, &status);
    }
  }

  int nsubint = 1e9;
  int ifile = 0;
  for (i=0; i<nIn; i++)
  {
    if ((force_merge == 0) && (nsubints[i] != nsubints[0]))
    {
      printf ("Error: Nsubints different in files %d vs %d\n", nsubints[0], nsubints[i]);
      return (EXIT_FAILURE);
    }
    if (nsubints[i] < nsubint)
    {
      nsubint = nsubints[i];
      ifile = i;
    }
  }

  printf("Output file is %s\n",outname);

  /* Open the input file */
  if ( !fits_open_file(&infptr, inname[ifile], READONLY, &status) )
  {
    // Copy the first file to the output file
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
    fits_close_file(infptr, &status);
  }

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
  dataArrayFP64 = (double *)malloc(sizeof(double)*nchan); 
  printf("Size of writeVals = %lld\n",sizeWriteVals);

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
  fits_modify_vector_len(outfptr,colnum_out_datFreq,newNchan,&status);
  if (status) {fits_report_error(stderr, status); exit(1);}
  
  fits_get_colnum(outfptr, CASEINSEN, "DAT_WTS", &colnum_out_datWts, &status);  
  fits_modify_vector_len(outfptr,colnum_out_datWts,newNchan,&status);
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
      
            printf("Reading %d subints\n", nsubint);
      for (k=0;k<nsubint;k++)
        {
          writePos = 1+i*nchan;
          
          fits_read_col(infptr,TDOUBLE,colnum_in_datFreq,k+1,1,nchan,&nullVal_d,dataArrayFP64,&initflag,&status);
          if (i==0)
            freqFirst = dataArrayFP64[0];
          if (i==nIn-1)
            freqLast = dataArrayFP64[nchan-1];
          fits_write_col(outfptr,TDOUBLE,colnum_out_datFreq,k+1,writePos,nchan,dataArrayFP64,&status);
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
  printf("Frequency range = %lf to %lf\n", freqFirst, freqLast);
  newObsFreq = fabs((freqLast + freqFirst)/2.0);
  fits_update_key(outfptr, TDOUBLE, (char *)"OBSFREQ", &newObsFreq, NULL, &status );
  fits_update_key(outfptr, TDOUBLE, (char *)"OBSBW", &newObsBw, NULL, &status );

  fits_movnam_hdu(outfptr, BINARY_TBL,(char *)"SUBINT",0,&status);
  fits_update_key(outfptr, TDOUBLE, (char *)"REFFREQ", &newObsFreq, NULL, &status );

  
    // Deallocate memory
    free(dataVals);
    free(dataArray);
    free(dataArrayFP64);
    free(writeVals);
    printf("Closing the output file\n");
    // Now close the file
    fits_close_file(outfptr,  &status);
    if (status) fits_report_error(stderr, status);
 
  return 0;   
}

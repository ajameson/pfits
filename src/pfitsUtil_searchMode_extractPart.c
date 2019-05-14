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
  
  int nchan,nsblk,npol,nbit;
  long naxes[4];
  int naxis=3;
  int newNchan;
  char *dataVals;
  char *writeVals;
  char nullVal = 0;
  float nullVal_f = 0;
  float *dataArray;
  float *freqChan;
  int subint;
  int subint_in;
  int subint_out;
  int initflag=0;
  float bytespersample=1;
  unsigned long writePos;
  char cval[16];
  int nsubint=0;
  char tdim[16];
  int data_colnum;
  FILE *fout;
  
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
  
  fout = fopen(outname,"w");

  nchan = 512;
  nsblk = 1024;
  freqChan = (float *)malloc(sizeof(float)*nchan);
  dataVals = (char *)malloc(sizeof(char)*nchan*nsblk);
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

      subint = 1;
      fits_read_col(infptr,TFLOAT,colnum_in_datFreq,subint,1,nchan,&nullVal_f,freqChan,&initflag,&status);

      p = 0; // Only 1 pol
      //
      // DAT_OFFS and DAT_SCL also change with polarisation
      //      writePos = 1+p*newNchan+i*nchan;
      
      //      fits_read_col(infptr,TFLOAT,colnum_in_datScl,k+1,1+p*nchan,nchan,&nullVal_f,dataArray,&initflag,&status);
      //      fits_read_col(infptr,TFLOAT,colnum_in_datOffs,k+1,1+p*nchan,nchan,&nullVal_f,dataArray,&initflag,&status);
      
      fits_read_col(infptr,TBYTE,colnum_data_in,subint,1+(p*nchan*nsblk)*bytespersample,nchan*bytespersample*nsblk,&nullVal,dataVals,&initflag,&status);
      for (k=0;k<nchan;k++)
	{	  
	  for (j=0;j<nsblk;j++)
	    {
	      fprintf(fout,"%d %g %d\n",j,freqChan[k],(int)dataVals[j*nchan+k]);
	    }
	  fprintf(fout,"\n");
	}
      fits_close_file(infptr, &status);	      
    }
  fclose(fout);
  printf("Writing the data\n");
  for (i=0;i<nsubint;i++)
    fits_write_col(outfptr,TBYTE,colnum_data_out,i+1,1,npol*newNchan*nsblk*bytespersample,writeVals+(int)(i*newNchan*nsblk*npol*bytespersample),&status);

    // Deallocate memory
  free(dataVals);
    free(freqChan);
    //    free(writeVals);
    printf("Closing the output file\n");
    // Now close the file
    
}

// Routines to add some recent data files into a pre-prepared template file
//
// HAVE TO ROTATE THE PROFILES FIRST TO ACCOUNT FOR DIFFERENT PREDICTOR BEFORE COMBINING
// SHOULD USE A SINGLE PREDICTOR FOR ALL OF ANDREW'S FILES?
// gcc -lm -o pfitsUtil_foldmode_useTemplate pfitsUtil_foldmode_useTemplate.c -lcfitsio


#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <fitsio.h>

void setHeaderInformation(fitsfile *outfptr,char *fname,double *stt_offs_ref);
void copyData(fitsfile *outfptr,char *fname,double stt_offs_ref);
void rotateData(int16_t *dataVals,int16_t *rotDataVals,int nbin,int nchan,double dat_offs,double dat_offs_ref,double period);
void getPredictor(char *fname);

int main(int argc,char *argv[])
{
  char outname[1024];
  char templateFile[1024];
  char fname[26][1024];
  fitsfile *infptr,*outfptr;
  int status=0;
  char usename[1024];
  int i;
  int hdu=1;
  int nFile=0;
  double stt_offs_ref=0;
  
  for (i=1;i<argc;i++)
    {
      if (strcmp(argv[i],"-o")==0)
	strcpy(outname,argv[++i]);
      else if (strcmp(argv[i],"-t")==0)
	strcpy(templateFile,argv[++i]);
      else
	strcpy(fname[nFile++],argv[i]);
    }
  sprintf(usename,"!%s",outname);

  printf("Creating output file\n");
  
  // Copy the first file to the output file
  /* Open the input file */
  if ( !fits_open_file(&infptr, templateFile, READONLY, &status) )
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
      fits_close_file(infptr, &status);
    }
  
  /* if error occured, print out error message */
  if (status) fits_report_error(stderr, status);
  printf("Setting header information\n");

  // Get final predictor
  getPredictor(fname[0]);

  // Set up header information from the first input file
  setHeaderInformation(outfptr,fname[0],&stt_offs_ref);
  
  for (i=0;i<nFile;i++)
    copyData(outfptr,fname[i],stt_offs_ref);
  
  fits_close_file(outfptr,&status);    
}


void getPredictor(char *fname)
{
  fitsfile *infptr;
  int status=0;
  int colnum;
  char nval[1024]="UNKNOWN";
  char **line;
  int anynul=0;
  long nrows;
  int i;
  FILE *fout;
  char str[4096];
  long double mjd0,mjd1;
  double freqLow,freqHigh;
  int ntimecoeff;
  int nfreqcoeff;
  int seglen;

  freqLow = 700;
  freqHigh = 4100;

  ntimecoeff = 12;
  nfreqcoeff = 2;
  seglen     = 48000; // Same defaults as T2psrcoeff 
  
  line = (char **)malloc(sizeof(char *));
  line[0] = (char *)malloc(sizeof(char)*1024);     
       
  if ( !fits_open_file(&infptr, fname, READONLY, &status) )
    {
      fits_movnam_hdu(infptr,BINARY_TBL,(char *)"PSRPARAM",0,&status);
      fits_get_colnum(infptr,CASEINSEN,"PARAM",&colnum,&status);        
      fits_get_num_rows(infptr,&nrows,&status);
      fout = fopen("pfitsUtil_ephem.par","w");
      for (i=1;i<=nrows;i++)
	{
	  fits_read_col_str(infptr,colnum,i,1,1,nval,line,&anynul,&status);
	  fprintf(fout,"%s\n",line[0]);
	}
      fclose(fout);

      // Create a predictor covering the entire band
      mjd0 = 58277.186874884299;
      mjd1 = 58277.22854155096666;
      sprintf(str,"tempo2 -pred \"PKS %.15Lf %.15Lf %lf %lf %d %d %d\" -f pfitsUtil_ephem.par",mjd0,mjd1,freqLow,freqHigh,ntimecoeff,nfreqcoeff,seglen);
      system(str);
    }

  free(line[0]);
  free(line);

}

void copyData(fitsfile *outfptr,char *fname,double stt_offs_ref)
{
  fitsfile *infptr;
  int status=0;
  int colnum_data_in,colnum_in_datFreq,colnum_in_datWts,colnum_in_datScl;
  int colnum_in_datOffs;
  int colnum_data_out;
  int i,j,p;
  long writePos;
  
  int npol = 4;
  int nchan = 128;
  int nbin = 1024;
  int subint_in = 1;
  int subint_out = 1;
  int nchanTot = nchan*26;
  int initflag=0;
  double obsFreq;
  long freqOff;
  
  double f0 = 704;
  double f1 = 4032;
  double absChanBW = (f1-f0)/nchanTot;
  double stt_offs;
  double period;
  
  int16_t *dataVals,*rotDataVals;
  int16_t nullVal=0;
  
  fits_movnam_hdu(outfptr,BINARY_TBL,(char *)"SUBINT",0,&status);
  fits_get_colnum(outfptr, CASEINSEN, "DATA", &colnum_data_out, &status);  

  fits_open_file(&infptr, fname, READONLY, &status);
  if (status) fits_report_error(stderr, status);

  // What is the frequency range for this file
  fits_read_key(infptr,TDOUBLE,"OBSFREQ",&obsFreq,NULL,&status);  
  fits_read_key(infptr,TDOUBLE,"STT_OFFS",&stt_offs,NULL,&status);  

  printf("Obs freq = %g\n",obsFreq);
  printf("Offset = %g versus %g, difference is %g\n",stt_offs,stt_offs_ref,stt_offs-stt_offs_ref);
  
  fits_movnam_hdu(infptr,BINARY_TBL,(char *)"SUBINT",0,&status);
  fits_get_colnum(infptr, CASEINSEN, "DATA", &colnum_data_in, &status);  
  fits_get_colnum(infptr, CASEINSEN, "DAT_FREQ", &colnum_in_datFreq, &status);  
  fits_get_colnum(infptr, CASEINSEN, "DAT_WTS", &colnum_in_datWts, &status);  
  fits_get_colnum(infptr, CASEINSEN, "DAT_SCL", &colnum_in_datScl, &status);  
  fits_get_colnum(infptr, CASEINSEN, "DAT_OFFS", &colnum_in_datOffs, &status);  
  if (status) fits_report_error(stderr, status);

  dataVals = (int16_t *)malloc(sizeof(int16_t)*nchan*nbin);
  rotDataVals = (int16_t *)malloc(sizeof(int16_t)*nchan*nbin);

  printf("PERIOD IS INCORRECTLY SET\n");
  period = 0.089403949593904; // MUST GET THIS FROM THE PREDICTOR

  if (obsFreq >= f0 && obsFreq <= f1)
    {
      for (p=0;p<npol;p++)
	{
	  fits_read_col(infptr,TSHORT,colnum_data_in,subint_in,1+p*nchan*nbin,nbin*nchan,&nullVal,dataVals,&initflag,&status);
	  rotateData(dataVals,rotDataVals,nbin,nchan,stt_offs,stt_offs_ref,period);      
	  freqOff = nchanTot - (obsFreq-128/2.0-f0)/absChanBW;      
	  writePos = 1+p*nchanTot*nbin + freqOff*nbin; 
	  printf("writePos = %d, obsFreq = %g, f0 = %g, freqOff = %d\n",writePos,obsFreq,f0,freqOff);
	  fits_write_col(outfptr,TSHORT,colnum_data_out,subint_out,writePos,nchan*nbin,rotDataVals,&status);
	}
    }
	
  free(dataVals);
  free(rotDataVals);
}

void rotateData(int16_t *dataVals,int16_t *rotDataVals,int nbin,int nchan,double dat_offs,double dat_offs_ref,double period)
{
  double offsetBinV;
  int offsetBinI;
  int i,j,readPos;
  int newBin;
  
  offsetBinV = nbin*(dat_offs-dat_offs_ref)/period;
  offsetBinI = (int)(offsetBinV + 0.5); // Should do this by a rotation
  //  offsetBinI = 0;

  offsetBinI = 0;
  printf("offsets are %g %d\n",offsetBinV,offsetBinI);
  for (i=0;i<nchan;i++)
    {     
      for (j=0;j<nbin;j++)
	{
	  newBin = j+offsetBinI;
	  while (newBin < 0) newBin+=nbin;
	  while (newBin >= nbin) newBin-=nbin;
	  //	  printf("newbin = %d %d\n",i,newBin);
	  
	  rotDataVals[i*nbin+newBin] = dataVals[i*nbin+j];
	}
    }
}

void setHeaderInformation(fitsfile *outfptr,char *fname,double *stt_offs_ref)
{
  int status=0;
  fitsfile *infptr;
  char str[1024];
  double dval;
  int newHDUtype;
  
  if ( !fits_open_file(&infptr, fname, READONLY, &status) )
    {
      fits_movabs_hdu(infptr, 1, NULL, &status);
      fits_movabs_hdu(outfptr, 1, NULL, &status);
      fits_read_key(infptr,TSTRING,"OBSERVER",str,NULL,&status);
      fits_update_key(outfptr, TSTRING, (char *)"OBSERVER", str, NULL, &status );  
      fits_read_key(infptr,TSTRING,"PROJID",str,NULL,&status);
      fits_update_key(outfptr, TSTRING, (char *)"PROJID", str, NULL, &status );  
      fits_read_key(infptr,TSTRING,"DATE-OBS",str,NULL,&status);
      fits_update_key(outfptr, TSTRING, (char *)"DATE-OBS", str, NULL, &status );  
      fits_read_key(infptr,TSTRING,"SRC_NAME",str,NULL,&status);
      fits_update_key(outfptr, TSTRING, (char *)"SRC_NAME", str, NULL, &status );  
      fits_read_key(infptr,TSTRING,"RA",str,NULL,&status);
      fits_update_key(outfptr, TSTRING, (char *)"RA", str, NULL, &status );  
      fits_read_key(infptr,TSTRING,"DEC",str,NULL,&status);
      fits_update_key(outfptr, TSTRING, (char *)"DEC", str, NULL, &status );  
      fits_read_key(infptr,TSTRING,"OBS_MODE",str,NULL,&status);
      fits_update_key(outfptr, TSTRING, (char *)"OBS_MODE", str, NULL, &status );  
      fits_read_key(infptr,TSTRING,"FD_POLN",str,NULL,&status);
      fits_update_key(outfptr, TSTRING, (char *)"FD_POLN", str, NULL, &status );  
      fits_read_key(infptr,TDOUBLE,(char *)"STT_IMJD",&dval,NULL,&status);
      fits_update_key(outfptr, TDOUBLE, (char *)"STT_IMJD", &dval, NULL, &status );  
      fits_read_key(infptr,TDOUBLE,(char *)"STT_SMJD",&dval,NULL,&status);
      fits_update_key(outfptr, TDOUBLE, (char *)"STT_SMJD", &dval, NULL, &status );  
      fits_read_key(infptr,TDOUBLE,(char *)"STT_OFFS",&dval,NULL,&status);
      fits_update_key(outfptr, TDOUBLE, (char *)"STT_OFFS", &dval, NULL, &status );  
      *stt_offs_ref = dval; // Set a reference offset
      
      // Copy the PSRPARAM table in its entirety
      printf("COPYING PSRPARAM TABLE\n");
      fits_movnam_hdu(infptr,BINARY_TBL,(char *)"PSRPARAM",0,&status);
      fits_movnam_hdu(outfptr,BINARY_TBL,(char *)"PSRPARAM",0,&status);
      fits_delete_hdu(outfptr, &newHDUtype, &status);
      fits_copy_hdu(infptr, outfptr, 0, &status);

      // Copy the HISTORY table in its entirety
      fits_movnam_hdu(infptr,BINARY_TBL,(char *)"HISTORY",0,&status);
      fits_movnam_hdu(outfptr,BINARY_TBL,(char *)"HISTORY",0,&status);
      fits_delete_hdu(outfptr, &newHDUtype, &status);
      fits_copy_hdu(infptr, outfptr, 0, &status);

      //
      fits_movnam_hdu(outfptr,BINARY_TBL,(char *)"POLYCO",0,&status);
      fits_delete_hdu(outfptr, &newHDUtype, &status);
      fits_movnam_hdu(outfptr,BINARY_TBL,(char *)"SPECKURT",0,&status);
      fits_delete_hdu(outfptr, &newHDUtype, &status);

      //
      fits_movnam_hdu(infptr,BINARY_TBL,(char *)"T2PREDICT",0,&status);
      fits_movnam_hdu(outfptr,BINARY_TBL,(char *)"T2PREDICT",0,&status);
      fits_delete_hdu(outfptr, &newHDUtype, &status);
      fits_copy_hdu(infptr, outfptr, 0, &status);
      
      // Copy information from SUBINT header
      fits_movnam_hdu(infptr,BINARY_TBL,(char *)"SUBINT",0,&status);
      fits_movnam_hdu(outfptr,BINARY_TBL,(char *)"SUBINT",0,&status);
      fits_read_key(infptr,TDOUBLE,(char *)"DM",&dval,NULL,&status);
      fits_update_key(outfptr, TDOUBLE, (char *)"DM", &dval, NULL, &status );  
      //      fits_read_key(infptr,TDOUBLE,(char *)"TBIN",&dval,NULL,&status);
      //      fits_update_key(outfptr, TDOUBLE, (char *)"TBIN", &dval, NULL, &status );  
      dval = 0;
      fits_update_key(outfptr, TDOUBLE, (char *)"TBIN", &dval, NULL, &status );  
      fits_close_file(infptr,&status);
    }


  if (status) fits_report_error(stderr, status);
  

}

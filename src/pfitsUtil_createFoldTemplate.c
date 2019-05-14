// Routines to produce an empty PSRFITS fold mode file
// containing the bands required for the Parkes UWL receiver system
//
// G. Hobbs
//
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <fitsio.h>

void writePrimaryHeaderParameters(fitsfile *fp,int npol,int nbin,int nsub,
				  int nsubBands,int chanPerSubband);
void createSubint(fitsfile *fp,int sub,int npol,int nbin,int nsub,
		  int nsubBands,int chanPerSubband);


int main(int argc,char *argv[])
{
  int  nsubBands = 26;
  int chanPerSubband = 128;
  int  npol = 4;
  int  nbin = 1024;
  int  nsub = 1;
  char fname[1024];
  fitsfile *fp;
  int status=0;
  int i;
  
  sprintf(fname,"!template_fold.rf(psrheader.fits)");
  fits_create_file(&fp,fname,&status); fits_report_error(stdout,status);
  // Write default header parameters
  writePrimaryHeaderParameters(fp,npol,nbin,nsub,nsubBands,chanPerSubband);

  for (i=0;i<nsub;i++)
    createSubint(fp,i,npol,nbin,nsub,nsubBands,chanPerSubband);
  
  fits_close_file(fp,&status); fits_report_error(stdout,status);
}

void writePrimaryHeaderParameters(fitsfile *fp,int npol,int nbin,int nsub,
				  int nsubBands,int chanPerSubband)
{
  double tVal;
  int iVal;
  int status=0;
  
  fits_movabs_hdu(fp, 1, NULL, &status );
  fits_write_date(fp, &status);
  fits_update_key(fp, TSTRING, (char *)"OBSERVER", (char *)"UNSET", NULL, &status );
  fits_update_key(fp, TSTRING, (char *)"PROJID", (char *)"UNSET", NULL, &status );
  fits_update_key(fp, TSTRING, (char *)"TELESCOP", (char *)"PARKES", NULL, &status );
  tVal = -4554231.6; fits_update_key(fp, TDOUBLE, (char *)"ANT_X", &tVal,NULL,&status);
  tVal = 2816759.1; fits_update_key(fp, TDOUBLE, (char *)"ANT_Y", &tVal,NULL,&status);
  tVal = -3454036.1; fits_update_key(fp, TDOUBLE, (char *)"ANT_Z", &tVal,NULL,&status);
  fits_update_key(fp, TSTRING, (char *)"FRONTEND", (char *)"UWL", NULL, &status );
  fits_update_key(fp, TSTRING, (char *)"BACKEND", (char *)"Medusa", NULL, &status );
  fits_update_key(fp, TSTRING, (char *)"OBS_MODE", (char *)"CAL", NULL, &status );
  tVal = 1; fits_update_key(fp, TDOUBLE, (char *)"CAL_FREQ", &tVal,NULL,&status); // DEFAULT VALUE
  fits_update_key(fp, TSTRING, (char *)"CAL_MODE", (char *)"SYNC", NULL, &status );
  tVal = 0.5*(4032.0+704); fits_update_key(fp, TDOUBLE, (char *)"OBSFREQ", &tVal,NULL,&status);
  tVal = 4032.0-704.; fits_update_key(fp, TDOUBLE, (char *)"OBSBW", &tVal,NULL,&status);
  iVal = chanPerSubband*nsubBands; fits_update_key(fp, TINT, (char *)"OBSNCHAN", &iVal,NULL,&status);
  fits_update_key(fp, TSTRING, (char *)"SRC_NAME", (char *)"UNSET", NULL, &status );
  fits_update_key(fp, TSTRING, (char *)"COORD_MD", (char *)"J2000", NULL, &status );
  tVal = 2000.; fits_update_key(fp, TDOUBLE, (char *)"EQUINOX", &tVal,NULL,&status);
  tVal = 53000.; fits_update_key(fp, TDOUBLE, (char *)"STT_IMJD", &tVal,NULL,&status); // DEFAULT VAL
  tVal = 0.; fits_update_key(fp, TDOUBLE, (char *)"STT_SMJD", &tVal,NULL,&status); // DEFAULT VAL
  tVal = 0.; fits_update_key(fp, TDOUBLE, (char *)"STT_OFFS", &tVal,NULL,&status); // DEFAULT VAL
}

void createSubint(fitsfile *fp,int sub,int npol,int nbin,int nsub,
		  int nsubBands,int chanPerSubband)
{
  int status=0;
  int ncol;
  int subint = sub+1;
  char tdim[64];
  long naxes[4];
  int naxis=3;
  int nchanTotal = nsubBands*chanPerSubband;
  double chan_bw = -128.0/chanPerSubband;
  float *dat_freq,*dat_wts,*dat_offs,*dat_scl,*dataVals;
  int colnum;

  double tsub = 1; // Default sub-int time
  double offsSub = 0.5; //Default offset time
  double fmax,fmin;
  int indxval = 0;
  int i;
  
  fits_movnam_hdu(fp,BINARY_TBL,(char *)"SUBINT",0,&status);
  if (status) {fits_report_error(stdout,status); exit(1);}
  fits_update_key(fp, TINT, (char *)"NAXIS2", &subint, NULL, &status );
  if (status) {fits_report_error(stdout,status); exit(1);}
  fits_get_num_cols(fp,&ncol,&status);

  if (sub==0) // First sub-integration
    {
      int dval_0 = 0;
      int dval_1 = 1;

      
      fits_update_key(fp, TSTRING, (char *)"INT_TYPE", (char *)"TIME", NULL, &status );
      fits_update_key(fp, TSTRING, (char *)"INT_UNIT", (char *)"SEC", NULL, &status );
      fits_update_key(fp, TSTRING, (char *)"SCALE", (char *)"FluxDen", NULL, &status );
      if (npol==4) fits_update_key(fp, TSTRING, (char *)"POL_TYPE", (char *)"AABBCRCI", NULL, &status );
      else if (npol==1) fits_update_key(fp, TSTRING, (char *)"POL_TYPE", (char *)"AA+BB", NULL, &status );
      else if (npol==2) fits_update_key(fp, TSTRING, (char *)"POL_TYPE", (char *)"AA,BB", NULL, &status );
      else {printf("ERROR: Not sure what to do with NPOL = %d\n",npol); exit(1);}
      
      fits_update_key(fp, TINT, (char *)"NPOL",&npol,NULL,&status);  
      fits_update_key(fp, TINT, (char *)"NBIN",&nbin,NULL,&status);  
      fits_update_key(fp, TINT, (char *)"NBIN_PRD",&nbin,NULL,&status);  
      fits_update_key(fp, TINT, (char *)"PHS_OFFS",&dval_0,NULL,&status);  
      fits_update_key(fp, TINT, (char *)"NBITS",&dval_1,NULL,&status);  
      fits_update_key(fp, TINT, (char *)"ZERO_OFF",&dval_0,NULL,&status);  
      fits_update_key(fp, TINT, (char *)"NSUBOFFS",&dval_0,NULL,&status);  
      fits_update_key(fp, TINT, (char *)"NCHAN",&nchanTotal,NULL,&status);  
      fits_update_key(fp, TDOUBLE, (char *)"CHAN_BW",&chan_bw,NULL,&status);  
      //      fits_update_key(fp, TDOUBLE, (char *)"TBIN",&tbin,NULL,&status);  
      //      fits_update_key(fp, TDOUBLE, (char *)"DM",&(control->dm),NULL,&status);  
      //      fits_update_key(fp, TINT, (char *)"RM",&dval_0,NULL,&status);  
      fits_update_key(fp, TINT, (char *)"NCHNOFFS",&dval_0,NULL,&status);  
      fits_update_key(fp, TINT, (char *)"NSBLK",&dval_1,NULL,&status);  
    }

  dat_freq = (float *)malloc(sizeof(float)*nchanTotal);
  dat_wts = (float *)malloc(sizeof(float)*nchanTotal);
  dat_offs = (float *)malloc(sizeof(float)*nchanTotal*npol);
  dat_scl = (float *)malloc(sizeof(float)*nchanTotal*npol);
  dataVals = (float *)malloc(sizeof(float)*nchanTotal*nbin*npol);

  fits_get_colnum(fp,CASEINSEN,"TSUBINT",&colnum,&status);
  fits_write_col(fp,TDOUBLE,colnum,subint,1,1,&tsub,&status);
  fits_get_colnum(fp,CASEINSEN,"OFFS_SUB",&colnum,&status);
  fits_write_col(fp,TDOUBLE,colnum,subint,1,1,&offsSub,&status);

  fmax = 4032; fmin = 704.0;
  
  for (i=0;i<nchanTotal;i++)
    {
      dat_freq[i] = fmax + chan_bw*i + chan_bw/2.0;
      dat_wts[i] = 1;
    }
  for (i=0;i<nchanTotal*npol;i++)
    {
      dat_offs[i] = 0;
      dat_scl[i] = 1;
    }
  for (i=0;i<nchanTotal*npol*nbin;i++)
    dataVals[i] = 0.;

  fits_get_colnum(fp,CASEINSEN,"INDEXVAL",&colnum,&status);
  fits_write_col(fp,TINT,colnum,subint,1,1,&indxval,&status);

  fits_get_colnum(fp, CASEINSEN, "DATA", &colnum, &status);  
  fits_modify_vector_len (fp, colnum, (nchanTotal*npol*nbin), &status); 
  if (status) {fits_report_error(stdout,status); exit(1);}
  naxes[0] = nbin;
  naxes[1] = nchanTotal;
  naxes[2] = npol;
  sprintf(tdim,"TDIM%d",colnum);
  fits_delete_key(fp, tdim, &status);
  fits_write_tdim(fp, colnum, naxis, naxes, &status);
  fits_write_col(fp,TFLOAT,colnum,subint,1,nchanTotal*npol*nbin,dataVals,&status);

  fits_get_colnum(fp, CASEINSEN, "DAT_FREQ", &colnum, &status);
  fits_modify_vector_len (fp, colnum, nchanTotal, &status); 
  fits_write_col(fp,TFLOAT,colnum,subint,1,nchanTotal,dat_freq,&status);
  fits_get_colnum(fp, CASEINSEN, "DAT_WTS", &colnum, &status);
  fits_modify_vector_len (fp, colnum, nchanTotal, &status); 
  fits_write_col(fp,TFLOAT,colnum,subint,1,nchanTotal,dat_wts,&status);
  fits_get_colnum(fp, CASEINSEN, "DAT_OFFS", &colnum, &status);
  fits_modify_vector_len (fp, colnum, nchanTotal*npol, &status); 
  fits_write_col(fp,TFLOAT,colnum,subint,1,nchanTotal*npol,dat_offs,&status);
  fits_get_colnum(fp, CASEINSEN, "DAT_SCL", &colnum, &status);
  fits_modify_vector_len (fp, colnum, nchanTotal*npol, &status); 
  fits_write_col(fp,TFLOAT,colnum,subint,1,nchanTotal*npol,dat_scl,&status);
  

  
  free(dat_freq);
  free(dat_wts);
  free(dat_offs);
  free(dat_scl);
  free(dataVals);
}

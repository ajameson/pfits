/* Software to simulate a PSRFITS search-mode file */
/*
 * Options include:
 * 
 * - white noise 
 * - pulse train 
 * - effect of dispersion 
 * - realistic pulse modulation 
 * - low level red noise (NOT DONE)
 * - impulsive RFI (NOT DONE)
 * - effect of pulse shape variation across the band (NOT DONE) --- use multiple components to define the pulse
 * - effect of scintillation (NOT DONE)
 * - effect of level setting (NOT DONE)
 * - signal in multiple telescope beams (NOT DONE)
 * - simulating binary pulsar (NOT DONE)
 * - include dead channels
 *
 * The user can select
 *
 * - sample time
 * - nsblk
 * - nchan
 * - rms of the white noise
 * - parameters of the pulse shape
 * - dispersion measure
 * - parameters of the red noise
 * - parameters of the pulse shape variation
 * - information of RFI
 * - smoothing time constant
 * - nbits
 */

#include <stdio.h>
#include <math.h>
#include <string.h>
#include "T2toolkit.h"
#include "fitsio.h"

void digitise(double *val,int nchan,int nbit,unsigned char *cvals,double offset,float *smmean);
void createHeader(fitsfile *fp,double freq_c,double bw,int nchan,char *raj,char *decj,
		  float scanLen,int stt_imjd,double stt_smjd);
void createSubintHeader(fitsfile *fp,double tbin,int nbits,int nchan,double chan_bw,int nsblk);
void addSubint(fitsfile *fp,int snum,int nchan,double freq_c,double bw,unsigned char *blkData,int nbit,int nsblk);

int main(int argc, char *argv[])
{
  int nchan = 64;
  double dm = 200;
  double tsamp = 1e-3;
  double tspan = 4.096*10;
  int nsblk = 4096;
  double white = 2;
  double offset = 128;
  double pulse_w = 0.01;
  double pulse_amp = 15;
  int i,j;
  char sig[nchan][nsblk];
  double val[nchan],t;
  double phase;
  double period=0.2;
  double phi0=0.5;
  double phi;
  long seed = TKsetSeed();
  int nbit = 8;
  unsigned char cvals[nchan*nsblk];
  int status=0;
  double freq_c = 1440;
  double bw     = 256;
  char raj[128] = "12:12:00.1";
  char decj[128] = "-05:12:00.1";
  int stt_imjd = 56000.0;
  double stt_smjd = 0.123456;
  double tDM,f0,f1;
  double smconst=0.1; // Smoothing time constant
  int smbin= smconst/tsamp; // Smoothing number of bins
  float run_meanval[nchan][smbin];
  float sm_mean[nchan];
  float prevVal,curVal;
  float *ptr;
  float *ptr_e;

  int s=0;
  int newHDUtype;
  long pn;
  float pamp[(int)(tspan/period+1)];
  int ndead=0;
  int dead[2] = {40,20};
  int k,isdead;
  float addVal=0;
  int nsub = ceil(tspan/tsamp/nsblk);
  float levelChange[nchan];
  double chanbw;


  char fname[128]="simulate.sf(psrheader.fits)";
  fitsfile *fp;

  // Open the file
  fits_create_file(&fp,fname,&status);
  fits_report_error(stderr,status);
  if (status)
    {
      printf("Error opening file >%s<\n",fname);
      exit(1);
    }

  
  // Write the header for the psrfits file
  createHeader(fp,freq_c,bw,nchan,raj,decj,tspan,stt_imjd,stt_smjd);

  // Delete unwanted tables
  fits_movnam_hdu(fp, BINARY_TBL, "BANDPASS", 0, &status);  fits_delete_hdu(fp, &newHDUtype, &status);   fits_report_error(stdout,status);
  fits_movnam_hdu(fp, BINARY_TBL, "COHDDISP", 0, &status);  fits_delete_hdu(fp, &newHDUtype, &status);   fits_report_error(stdout,status);
  fits_movnam_hdu(fp, BINARY_TBL, "PSRPARAM", 0, &status);  fits_delete_hdu(fp, &newHDUtype, &status);   fits_report_error(stdout,status);
  fits_movnam_hdu(fp, BINARY_TBL, "POLYCO", 0, &status);  fits_delete_hdu(fp, &newHDUtype, &status);   fits_report_error(stdout,status);
  fits_movnam_hdu(fp, BINARY_TBL, "T2PREDICT", 0, &status);  fits_delete_hdu(fp, &newHDUtype, &status);   fits_report_error(stdout,status);
  fits_movnam_hdu(fp, BINARY_TBL, "FLUX_CAL", 0, &status);  fits_delete_hdu(fp, &newHDUtype, &status);   fits_report_error(stdout,status);
  fits_movnam_hdu(fp, BINARY_TBL, "CAL_POLN", 0, &status);  fits_delete_hdu(fp, &newHDUtype, &status);   fits_report_error(stdout,status);
  fits_movnam_hdu(fp, BINARY_TBL, "FEEDPAR", 0, &status);  fits_delete_hdu(fp, &newHDUtype, &status);   fits_report_error(stdout,status);
  //  fits_movnam_hdu(fp, BINARY_TBL, "DIG_STAT", 0, &status);  fits_delete_hdu(fp, &newHDUtype, &status);   fits_report_error(stdout,status);
  //  fits_movnam_hdu(fp, BINARY_TBL, "DIG_CNTS", 0, &status);  fits_delete_hdu(fp, &newHDUtype, &status);   fits_report_error(stdout,status);

  // Write the data
  createSubintHeader(fp,tsamp,nbit,nchan,bw/(double)nchan,nsblk);

  // Create pulse amplitudes
  for (i=0;i<tspan/period;i++)
    pamp[i] = pow(10,0.2*TKgaussDev(&seed));

  // Set up running mean
  ptr = run_meanval[0];
  ptr_e = run_meanval[0]+smbin;
  for (i=0;i<nchan;i++)
    levelChange[i]=0.0;
  
  // Must sort out the minus sign!! 
  chanbw = -bw/(double)nchan;

  //  nsub =1 ;
  for (s=0;s<nsub;s++)
    {
      printf("Subintegration: %d\n",s+1);
      for (j=0;j<nsblk;j++)
	{
	  

	  t = s*nsblk*tsamp+j*tsamp;
	  //if (t > 1.3) addVal=10;

	  for (i=0;i<nchan;i++)
	    {
	      // Check for dead channels
	      isdead=0;
	      for (k=0;k<ndead;k++)
		{
		  if (i==dead[k]) {isdead=1; break;}
		}
	      if (isdead==1)
		{
		  val[i]=0;
		  //		  printf("DEAD\n");
		}
	      else
		{
		  // Add white noise
		  val[i] = TKgaussDev(&seed)*white;
		  
		  // Add pulse train
		  //		  printf("have %g %g %g\n",freq_c,bw,chanbw);
		  f0 = freq_c+bw/2.0;
		  f1 = f0+i*chanbw;
		  
		  tDM = -(4.148808e-3*dm*(pow(f0/1000.0,-2)-pow(f1/1000.0,-2)));
		  phi = phi0;
		  phase = fmod(t-tDM+phi*period,period);
		  pn = (int)((t-tDM+phi*period)/period+0.5);
		  //	      printf("pnumber = %d\n",pn);
		  //	  printf("phase = %g\n",phase);
		  //addVal=0;
		  val[i] += pamp[pn]*pulse_amp*exp(-pow((phase-period/2),2)/2.0/pulse_w/pulse_w)+addVal; //-levelChange[i];
		  //		  if (i==1)
		  //		    printf("Have %g %g %g %d\n",val[i],sm_mean[i],levelChange[i],smbin);
		}

	      //  float sm_mean[nchan];
	      //  float prevVal;
	      prevVal = *(ptr+i*smbin);
	      curVal = val[i]; //+levelChange[i];
	      //	      val[i] -= sm_mean[i];
	      sm_mean[i] = sm_mean[i]+(curVal-prevVal)/smbin;
	      //	      if (i==0) printf("mean = %g\n",sm_mean[i]);
	      //	      levelChange[i]+=sm_mean[i];
	      *(ptr+i*smbin) = val[i]; //+levelChange[i];
	    }
	  ptr++;
	  if (ptr == ptr_e) ptr=run_meanval[0];
	  digitise(val,nchan,nbit,cvals+j*nchan/8*nbit,offset,sm_mean);
	  //      printf("val = %d\n",cvals[j*nchan]);
	}
      // Now write out the subint
      addSubint(fp,s,nchan,freq_c,bw,cvals,nbit,nsblk);
      }
  printf("pnumber = %d\n",pn);

  // Close the file
  printf("Writing the file\n");
  fits_close_file(fp,&status);
}

void addSubint(fitsfile *fp,int snum,int nchan,double freq_c,double bw,unsigned char *blkData,int nbit,int nsblk)
{
  int status = 0;
  int colnum;
  int i;
  float datFreq[nchan];
  float datWts[nchan];
  float datOffs[nchan];
  float datScl[nchan];
  int ival;
  int samplesperbyte = 8/nbit;
  long naxes[4];
  char cval[128];
  double chanbw = bw/(double)nchan;
  
  printf("Starting add subint\n");
  for (i=0;i<nchan;i++)
    {
      datFreq[i] = freq_c-bw/2.0 + i*chanbw;
      datWts[i] = datOffs[i] = datScl[i] = 1.0;
    }
  //  sprintf(cval,"%dX",nchan*nsblk);
  //  fits_update_key(fp, TSTRING, "TFORM20", cval, NULL, &status);
  //  fits_report_error(stderr,status);

  fits_report_error(stderr,status);
  printf("About to insert %d\n",snum);
  fits_insert_rows(fp,snum,1,&status);
  if (status)
    {
      fits_report_error(stderr,status);
      exit(1);
    }
  printf("Have inserted some rows\n");
  fits_get_colnum(fp,CASEINSEN,"INDEXVAL",&colnum,&status); fits_report_error(stderr,status);
  ival = snum+1;  fits_write_col(fp,TINT,colnum,snum+1,1,1,&ival,&status);
  if (status!=0)
    {
      fits_report_error(stderr,status);
      exit(1);
    }
  printf("Have written the first column\n");
  fits_get_colnum(fp,CASEINSEN,"DAT_FREQ",&colnum,&status); fits_report_error(stderr,status);
  fits_modify_vector_len(fp,colnum,nchan,&status); fits_report_error(stderr,status);
  fits_write_col(fp,TFLOAT,colnum,snum+1,1,nchan,datFreq,&status);
  if (status!=0)
    {
      fits_report_error(stderr,status);
      exit(1);
    }


  printf("Have written the first column 2\n");
  fits_get_colnum(fp,CASEINSEN,"DAT_WTS",&colnum,&status); fits_report_error(stderr,status);
  fits_modify_vector_len(fp,colnum,nchan,&status); fits_report_error(stderr,status);
  fits_write_col(fp,TFLOAT,colnum,snum+1,1,nchan,datWts,&status);
  printf("Have written the first column 3\n");
  fits_get_colnum(fp,CASEINSEN,"DAT_OFFS",&colnum,&status); fits_report_error(stderr,status);
  fits_modify_vector_len(fp,colnum,nchan,&status); fits_report_error(stderr,status);
  fits_write_col(fp,TFLOAT,colnum,snum+1,1,nchan,datOffs,&status);
  printf("Have written the first column 4\n");
  fits_get_colnum(fp,CASEINSEN,"DAT_SCL",&colnum,&status); fits_report_error(stderr,status);
  fits_modify_vector_len(fp,colnum,nchan,&status); fits_report_error(stderr,status);
  fits_write_col(fp,TFLOAT,colnum,snum+1,1,nchan,datScl,&status);

  printf("About towrite the data\n");
  // Write the data
  naxes[0] = nchan;
  naxes[1] = 1; // Npol
  naxes[2] = nsblk;

  fits_get_colnum(fp,CASEINSEN,"DATA",&colnum,&status); fits_report_error(stderr,status);

  fits_modify_vector_len (fp, colnum, (nchan*nsblk), &status); fits_report_error(stderr,status);
  printf("Writing tdim\n");
  fits_delete_key(fp, "TDIM20", &status); // THIS SHOULD NOT BE HARDCODED
  fits_write_tdim(fp,colnum,3,naxes,&status);fits_report_error(stderr,status);
  printf("Done write\n");
  fits_write_col_byt(fp,colnum,snum+1,1,nsblk*nchan/samplesperbyte,blkData,&status);
  fits_report_error(stderr,status);
}


void createSubintHeader(fitsfile *fp,double tbin,int nbits,int nchan,double chan_bw,int nsblk)
{
  int status=0;
  char name[128];
  char str[128];
  char cval[128];

  strcpy(name,"SUBINT");
  fits_movnam_hdu(fp,BINARY_TBL,name,0,&status);fits_report_error(stderr,status);
  if (nbits==1)
    {
      sprintf(cval,"16X");//,nchan*nsblk);
      //      sprintf(cval,"%dX",nchan*nsblk);
      fits_update_key(fp, TSTRING, "TFORM20", cval, NULL, &status);
      fits_report_error(stderr,status);
    }
  printf("In the subint\n");
  fits_update_key_str(fp,"INT_TYPE","TIME","",&status);  fits_report_error(stderr,status);
  fits_update_key_str(fp,"INT_UNIT","SEC","",&status);  fits_report_error(stderr,status);
  fits_update_key_str(fp,"SCALE","FluxDen","",&status);  fits_report_error(stderr,status);
  fits_update_key_str(fp,"NPOL","1","",&status);  fits_report_error(stderr,status);
  fits_update_key_str(fp,"NBIN","1","",&status);  fits_report_error(stderr,status);
  fits_update_key_str(fp,"ZERO_OFF","0","",&status);  fits_report_error(stderr,status);
  fits_update_key_str(fp,"SIGNINT","0","",&status);  fits_report_error(stderr,status);
  fits_update_key_str(fp,"NSUBOFFS","0","",&status);  fits_report_error(stderr,status);

  sprintf(str,"%.3f",tbin); fits_update_key_str(fp,"TBIN",str,"",&status);  fits_report_error(stderr,status);
  sprintf(str,"%d",nchan); fits_update_key_str(fp,"NCHAN",str,"",&status);  fits_report_error(stderr,status);
  // NOTE: I HAVE A MINUS SIGN HERE ....
  sprintf(str,"%.3f",-chan_bw); fits_update_key_str(fp,"CHAN_BW",str,"",&status);  fits_report_error(stderr,status);
  sprintf(str,"%d",nsblk); fits_update_key_str(fp,"NSBLK",str,"",&status);  fits_report_error(stderr,status);
  sprintf(str,"%d",nbits); fits_update_key_str(fp,"NBITS",str,"",&status);  fits_report_error(stderr,status);
  printf("Complete creating header\n");
}

void digitise(double *val,int nchan,int nbit,unsigned char *cvals,double offset,float *smmean)
{
  int i,j;
  unsigned char tc;
  double bit_level = 0;

  if (nbit==8)
    {
      for (i=0;i<nchan;i++)
	cvals[i] = (char)(val[i]+offset-smmean[i]);
    }
  else if (nbit==1)
    {
      long n=0;
      for (i=0;i<nchan/(8/nbit);i++)
	{
	  tc=0;

	  for (j=0;j<8/nbit;j++)
	    {
	      //	      printf("Here with %g\n",val[i*(8/nbit)+j]);
	      //	      if (val[i*(8/nbit)+j] > bit_level)
	      	      //	      if (j==0 || j==3)
	      if (val[n]-smmean[i] > bit_level)
		//if (val[n] > bit_level)
		tc = tc | (1 << (7-j));
	      n++;
	    }
	  cvals[i] = tc;
	}
    }
  else
    {
      printf("Unable to process %d bit data\n",nbit);
    }
}


void createHeader(fitsfile *fp,double freq_c,double bw,int nchan,char *raj,char *decj,float scanLen,int stt_imjd,double stt_smjd)
{
  int status=0;
  char str[128];

  fits_movabs_hdu(fp,1,NULL,&status);


  fits_write_date(fp,&status);
  printf("Creating header\n");
  fits_update_key_str(fp,"OBSERVER","PFITS_SIM","",&status);  fits_report_error(stderr,status);
  printf("First line\n");
  fits_update_key_str(fp,"PROJID","P999","",&status);  fits_report_error(stderr,status);
  fits_update_key_str(fp,"TELESCOP","SIMULATE","",&status);  fits_report_error(stderr,status);
  fits_update_key_str(fp,"FRONTEND","SIMULATE","",&status);  fits_report_error(stderr,status);
  fits_update_key_str(fp,"IBEAM","1","",&status);  fits_report_error(stderr,status);
  fits_update_key_str(fp,"NRCVR","2","",&status);  fits_report_error(stderr,status);
  fits_update_key_str(fp,"BACKEND","SIMULATE","",&status);  fits_report_error(stderr,status);
  fits_update_key_str(fp,"OBS_MODE","SEARCH","",&status);  fits_report_error(stderr,status);
  fits_update_key_str(fp,"SRC_NAME","SIMULATE","",&status);  fits_report_error(stderr,status);
  fits_update_key_str(fp,"TRK_MODE","TRACK","",&status);  fits_report_error(stderr,status);
  fits_update_key_str(fp,"CAL_MODE","OFF","",&status);  fits_report_error(stderr,status);
  fits_update_key_str(fp,"CAL_FREQ","0","",&status);  fits_report_error(stderr,status);
  fits_update_key_str(fp,"CAL_DCYC","0","",&status);  fits_report_error(stderr,status);
  fits_update_key_str(fp,"CAL_PHS","0","",&status);  fits_report_error(stderr,status);
  fits_update_key_str(fp,"ANT_X","-4554231.5","",&status);  fits_report_error(stderr,status);
  fits_update_key_str(fp,"ANT_Y","2816759.1","",&status);  fits_report_error(stderr,status);
  fits_update_key_str(fp,"ANT_Z","-3454036.3","",&status);  fits_report_error(stderr,status);
  fits_update_key_str(fp,"COORD_MD","J2000","",&status);  fits_report_error(stderr,status);
  fits_update_key_str(fp,"EQUINOX","2000","",&status);  fits_report_error(stderr,status);
  fits_update_key_str(fp,"FD_HAND","1","",&status);  fits_report_error(stderr,status);
  fits_update_key_str(fp,"FD_SANG","-90","",&status);  fits_report_error(stderr,status);
  fits_update_key_str(fp,"FD_XYPH","0","",&status);  fits_report_error(stderr,status);
  fits_update_key_str(fp,"FD_POLN","LIN","",&status);  fits_report_error(stderr,status);
  fits_update_key_str(fp,"BMAJ","0.24","",&status);  fits_report_error(stderr,status);
  fits_update_key_str(fp,"BMIN","0.24","",&status);  fits_report_error(stderr,status);
  fits_update_key_str(fp,"BPA","0","",&status);  fits_report_error(stderr,status);
  fits_update_key_str(fp,"BE_PHASE","0","",&status);  fits_report_error(stderr,status);
  fits_update_key_str(fp,"BE_DCC","0","",&status);  fits_report_error(stderr,status);
  fits_update_key_str(fp,"TCYCLE","0","",&status);  fits_report_error(stderr,status);
  fits_update_key_str(fp,"BE_DELAY","0","",&status);  fits_report_error(stderr,status);
  fits_update_key_str(fp,"STT_OFFS","0","",&status);  fits_report_error(stderr,status);
  fits_update_key_str(fp,"RA",raj,"",&status);  fits_report_error(stderr,status);
  fits_update_key_str(fp,"STT_CRD1",raj,"",&status);  fits_report_error(stderr,status);
  fits_update_key_str(fp,"STP_CRD1",raj,"",&status);  fits_report_error(stderr,status);
  fits_update_key_str(fp,"DEC",decj,"",&status);  fits_report_error(stderr,status);
  fits_update_key_str(fp,"STT_CRD2",decj,"",&status);  fits_report_error(stderr,status);
  fits_update_key_str(fp,"STP_CRD2",decj,"",&status);  fits_report_error(stderr,status);

  fits_update_key_flt(fp,"OBSFREQ",freq_c,5,"",&status);  fits_report_error(stderr,status);
  fits_update_key_flt(fp,"OBSBW",bw,5,"",&status);  fits_report_error(stderr,status);
  fits_update_key_flt(fp,"OBSNCHAN",nchan,5,"",&status);  fits_report_error(stderr,status);
  fits_update_key_flt(fp,"SCANLEN",scanLen,6,"",&status);  fits_report_error(stderr,status);
  sprintf(str,"%d",stt_imjd); fits_update_key_str(fp,"STT_IMJD",str,"",&status);  fits_report_error(stderr,status);
  sprintf(str,"%.9f",stt_smjd); fits_update_key_str(fp,"STT_SMJD",str,"",&status);  fits_report_error(stderr,status);
}

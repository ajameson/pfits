// 58 seconds - dm 0
// Just reading - 3 seconds
// All in one function -- 2:01

// code for pfits

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <math.h>
#include "pfits.h"

void setMask(int maskN,int maskP,unsigned char **mask);
unsigned int extractBit(unsigned char byte,unsigned int pos);
//int extractBit(char byte,int pos);

void closeFitsFile(fitsfile *fp)
{
  int status=0;
  fits_close_file(fp,&status);
  fits_report_error(stderr,status);
  if (status)
    {
      printf("Error closing file\n");
      exit(1);
    }

}

fitsfile * openFitsFile(char *fname)
{
  fitsfile *fp;
  int status=0;
  fits_open_file(&fp,fname,READONLY,&status);
  fits_report_error(stderr,status);
  if (status)
    {
      printf("Error opening file >%s<\n",fname);
      exit(1);
    }
  return fp;
}

dSet * initialiseDset()
{
  dSet *data;
  if (!(data = (dSet *)malloc(sizeof(dSet))))
    {
      printf("Unable to allocate enough memory to create a data set\n");
      exit(1);
    }
  data->pheaderSet = 0;
  data->subintTable = 0;
  data->psrparamTable = 0;
  data->phead.dm = -1; // Set to impossible value as default
  data->phead.period = -1; // Set to impossible value as default
  return data;
}

void readPhead(dSet *data)
{
  printf("Reading header");
}

void freeDset(dSet *data)
{
  int i;

  if (data->pheaderSet == 1)
    {
      for (i=0;i<data->phead.nhead;i++)
    {
      free(data->phead.keyname[i]);
      free(data->phead.val[i]);
      free(data->phead.comment[i]);
    }
      free(data->phead.keyname);
      free(data->phead.val);
      free(data->phead.comment);
    }

  free(data);
}

void loadPrimaryHeader(fitsfile *fp,dSet *data)
{
  int status=0;
  int nkey=-1;
  int morekeys=-1;
  int i;
  char keyname[128],val[128],comment[128];

  fits_get_hdrspace(fp,&nkey,&morekeys,&status);
  data->pheaderSet = 1;

  data->phead.nhead = nkey;
  // Allocate memory
  data->phead.keyname = (char **)malloc(sizeof(char *)*nkey);
  data->phead.val = (char **)malloc(sizeof(char *)*nkey);
  data->phead.comment = (char **)malloc(sizeof(char *)*nkey);
  for (i=0;i<nkey;i++)
    {
      data->phead.keyname[i] = (char *)malloc(sizeof(char)*128);
      data->phead.val[i] = (char *)malloc(sizeof(char)*128);
      data->phead.comment[i] = (char *)malloc(sizeof(char)*128);
    }
  data->pheaderSet=1;

  // Complete allocating memory

  for (i=1;i<=nkey;i++)
    {
      fits_read_keyn(fp,i+1,data->phead.keyname[i-1],data->phead.val[i-1],data->phead.comment[i-1],&status);
      if (strcmp(data->phead.keyname[i-1],"OBSFREQ")==0)
    sscanf(data->phead.val[i-1],"%f",&(data->phead.freq));
      else if (strcmp(data->phead.keyname[i-1],"STT_IMJD")==0)
    sscanf(data->phead.val[i-1],"%d",&(data->phead.imjd));
      else if (strcmp(data->phead.keyname[i-1],"STT_SMJD")==0)
    sscanf(data->phead.val[i-1],"%f",&(data->phead.smjd));
      else if (strcmp(data->phead.keyname[i-1],"STT_OFFS")==0)
    sscanf(data->phead.val[i-1],"%f",&(data->phead.stt_offs));
      else if (strcmp(data->phead.keyname[i-1],"OBSBW")==0)
    sscanf(data->phead.val[i-1],"%f",&(data->phead.bw));
    }
  // Read specific parameters
  fits_read_key(fp,TSTRING,"OBS_MODE",data->phead.obsMode,NULL,&status);
  fits_read_key(fp,TSTRING,"SRC_NAME",data->phead.source,NULL,&status);
  if (status)
    {
      fits_report_error(stderr,status);
      exit(1);
    }

  // Now load information from the subintegration table
  fits_movnam_hdu(fp,BINARY_TBL,"SUBINT",1,&status);
  if (status)
    {
      printf("No subintegration table\n");
      data->subintTable=0;
      status=0;
    }
  else
    {
      data->subintTable=1;
      fits_read_key(fp,TINT,"NAXIS2",&(data->phead.nsub),NULL,&status);
      if (status)
    {
      printf("Reading naxis2\n");
      fits_report_error(stderr,status);
      exit(1);
    }     
      fits_read_key(fp,TINT,"NCHAN",&(data->phead.nchan),NULL,&status);
      if (status)
    {
      printf("Reading nchan\n");
      fits_report_error(stderr,status);
      exit(1);
    }
      
      fits_read_key(fp,TFLOAT,"ZERO_OFF",&(data->phead.zeroOff),NULL,&status);
      if (status)
    {
      printf("Reading zero_off\n");
      fits_report_error(stderr,status);
      data->phead.zeroOff = 0;
      status=0;
    }

      fits_read_key(fp,TINT,"NBITS",&(data->phead.nbits),NULL,&status);
      if (status)
    {
      printf("Reading nbits\n");
      fits_report_error(stderr,status);
      exit(1);
    }

      fits_read_key(fp,TINT,"NPOL",&(data->phead.npol),NULL,&status);
      if (status)
    {
      printf("Reading npol\n");
      fits_report_error(stderr,status);
      exit(1);
    }
      
      fits_read_key(fp,TINT,"NSBLK",&(data->phead.nsblk),NULL,&status);
      if (status)
    {
      printf("Reading nsblk\n");
      fits_report_error(stderr,status);
      exit(1);
    }

      fits_read_key(fp,TINT,"NBIN",&(data->phead.nbin),NULL,&status);
      if (status)
    {
      printf("Reading nbin\n");
      fits_report_error(stderr,status);
      exit(1);
    }

      //      printf("nbin = %d (%d)\n",data->phead.nbin,status);
      fits_read_key(fp,TFLOAT,"CHAN_BW",&(data->phead.chanbw),NULL,&status);
      if (data->phead.chanbw < 0 && data->phead.bw > 0)
    data->phead.bw*=-1;
      
      fits_read_key(fp,TFLOAT,"TBIN",&(data->phead.tsamp),NULL,&status);
      
    }
  fits_movnam_hdu(fp,BINARY_TBL,"PSRPARAM",1,&status);
  if (status)
    {
#ifdef PFITS_DEBUG
      printf("No PSRPARM table\n");
#endif
      data->psrparamTable=0;
      status=0;
    }
  else
    {
      int len,i,colnum;
      char **line,str1[1024],str2[1024];
      data->psrparamTable=1;
      char nval[128]="UNKNOWN";
      int anynul=0;
      float tt;
      fits_read_key(fp,TINT,"NAXIS2",&len,NULL,&status);

      fits_get_colnum(fp,CASEINSEN,"PARAM",&colnum,&status);
      if (status) {
    printf("Unable to find data in the psrparam table in FITS file\n");
    exit(1);
      }

      line = (char **)malloc(sizeof(char *));
      line[0] = (char *)malloc(sizeof(char)*1024); 

      for (i=0;i<len;i++)
    {
      fits_read_col_str(fp,colnum,i+1,1,1,nval,line,&anynul,&status);
      if (sscanf(line[0],"%s %s",str1,str2)==2)
        {
          if (strcasecmp(str1,"DM")==0)
        sscanf(str2,"%f",&(data->phead.dm));
          if (strcasecmp(str1,"F0")==0)
        {
          sscanf(str2,"%f",&tt);
          data->phead.period = 1.0/tt;
        }
        }
      //      printf("Read: %s\n",line[0]);
    }
      //      printf("Lenght = %d\n",len);
  free(line[0]);
  free(line);

    }
    
}

void displayHeaderInfo(dSet *data)
{
  if (data->pheaderSet == 0)
    printf("No header information loaded\n");
  else
    {
      printf("Source:                       %s\n",data->phead.source);
      printf("Observation mode:             %s\n",data->phead.obsMode);
      printf("Observation frequency:        %f\n",data->phead.freq);
      printf("Integer start time (d):       %d\n",data->phead.imjd);
      printf("Fractional start time (s):    %f\n",data->phead.smjd);
      printf("Start time offset (s):        %f\n",data->phead.stt_offs);
      printf("Observation bandwidth (MHz):  %f\n",data->phead.bw);
      printf("Number of channels:           %d\n",data->phead.nchan);
      printf("Number of subintegrations:    %d\n",data->phead.nsub);
      printf("Number of bits:               %d\n",data->phead.nbits);
      printf("Number of polarisations:      %d\n",data->phead.npol);
      printf("Number of samples per subint: %d\n",data->phead.nsblk);
      printf("Number of bins:               %d\n",data->phead.nbin);
      printf("Channel bandwidth:            %f\n",data->phead.chanbw);
      printf("Sample time:                  %f\n",data->phead.tsamp);
      printf("Zero offset:                  %f\n",data->phead.zeroOff);
      printf("Dispersion measure:           %f\n",data->phead.dm);
      printf("Period:                       %f\n",data->phead.period);
      printf("\n");
      printf("Observation time (s):         %g\n",data->phead.nsub*data->phead.nsblk*data->phead.tsamp);
    }
}

int extractFoldData(fitsfile *fp,dSet *data,float dm,float *fx,float *fy,float *freq_y,float *time_y,float *bpass, int sub0)
{
  int n=0;
  int status=0;
  int colnum;
  int i,j,k,l;
  int initflag=0;
  float nval=0;
  float ty[data->phead.nbin];
  float **offs;    // [data->phead.nsub];
  float **dat_scl; // [data->phead.nsub];
  double f0,chanbw,tdelay;
  int bn,cdelay;
  double bintime;
  int addDelay = 0; //500;
//  float bpass[data->phead.nchan*2];
  float bpass_offs[2];
  float bpass_scl[2];
  float meanVal,rmsVal;


  // get mean/RMS of off-pulse for scaling. no longer need OFFS/DAT_SCL
  printf("sub0 = %d\n",sub0);
  if (dm < 0)
    dm = data->phead.dm;
  
  // Need to remove a baseline from each polarisation channel before summing

  // Get first frequency channel
  // Central frequency
  f0 = data->phead.freq; //-data->phead.chanbw*data->phead.nchan/2;
  //  chanbw = data->phead.bw/data->phead.nchan;
  chanbw = data->phead.chanbw;

  bintime = (double)data->phead.period/(double)data->phead.nbin;

  fits_movnam_hdu(fp,BINARY_TBL,"BANDPASS",1,&status);
  if (status) {
    printf("Unable to move to bandpass table in FITS file\n");
    exit(1);
  }
  fits_get_colnum(fp,CASEINSEN,"DAT_OFFS",&colnum,&status);  
  fits_read_col_flt(fp,colnum,1,1,2,nval,bpass_offs,&initflag,&status);
  fits_get_colnum(fp,CASEINSEN,"DAT_SCL",&colnum,&status);  
  fits_read_col_flt(fp,colnum,1,1,2,nval,bpass_scl,&initflag,&status);
  // Now read the bandpass
  fits_get_colnum(fp,CASEINSEN,"DATA",&colnum,&status);  
  // NOTE: Starting at element 2 as element 1 is junk
  fits_read_col_flt(fp,colnum,1,2,data->phead.nchan,nval,bpass,&initflag,&status);
  fits_read_col_flt(fp,colnum,1,3+data->phead.nchan,data->phead.nchan,nval,bpass+data->phead.nchan,&initflag,&status);
  for (i=0;i<data->phead.nchan*2;i++)
    {
      if (i<data->phead.nchan)
    bpass[i] = bpass[i]*bpass_scl[0] + bpass_offs[0];
      else
    bpass[i] = bpass[i]*bpass_scl[1] + bpass_offs[1];
      //      printf("bpass: %d %g %g %g %g %g \n",i,bpass[i],bpass_scl[0],bpass_scl[1],bpass_offs[0],bpass_offs[1]);
    }
  //  exit(1);

  fits_movnam_hdu(fp,BINARY_TBL,"SUBINT",1,&status);
  if (status) {
    printf("Unable to move to subint table in FITS file\n");
    exit(1);
  }
  // REMOVED: No longer need dat_scl or offs. OFFS/DAT_SCL
//  offs = (float **)malloc(sizeof(float *)*data->phead.nsub);
//  dat_scl = (float **)malloc(sizeof(float *)*data->phead.nsub);
//  for (i=0;i<data->phead.nsub;i++)
//    {
//      offs[i] = (float *)malloc(sizeof(float)*data->phead.nchan*data->phead.npol);
//      dat_scl[i] = (float *)malloc(sizeof(float)*data->phead.nchan*data->phead.npol);
//    }
//  fits_get_colnum(fp,CASEINSEN,"DAT_OFFS",&colnum,&status);
//  if (status) {
//    printf("Unable to find DAT_OFFS in the subint table in FITS file\n");
//    exit(1);
//  }
//  for (i=0;i<data->phead.nsub;i++)
//    {
//      fits_read_col_flt(fp,colnum,i+1,1,data->phead.nchan*data->phead.npol,nval,offs[i],&initflag,&status);
//      //      printf("offs = %g\n",offs[i][5]);
//      //      offs[i] =0;
//    }
//
//  fits_get_colnum(fp,CASEINSEN,"DAT_SCL",&colnum,&status);
//  if (status) {
//    printf("Unable to find DAT_SCL in the subint table in FITS file\n");
//    exit(1);
//  }
//  for (i=0;i<data->phead.nsub;i++)
//    {
//      fits_read_col_flt(fp,colnum,i+1,1,data->phead.nchan*data->phead.npol,nval,dat_scl[i],&initflag,&status);
//      //      printf("dat_scl = %g\n",dat_scl[i][5]);
//      //            dat_scl[i]=1.0;
//    }

  fits_get_colnum(fp,CASEINSEN,"DATA",&colnum,&status);
  if (status) {
    printf("Unable to find data in the subint table in FITS file\n");
    exit(1);
  }
  if (sub0==0)
    {
      for (i=0;i<data->phead.nbin;i++)
    {
      fx[i] = i;
      fy[i] = 0;
    }
      //
      printf("Loading %d subintegrations\n",data->phead.nsub);
      printf("Number of frequency channels = %d\n",data->phead.nchan);
      printf("Number of polarisations = %d\n",data->phead.npol);
      for (i=0;i<data->phead.nchan*data->phead.nbin;i++)
    freq_y[i] = 0;
    }
  for (i=sub0;i<data->phead.nsub;i++) // *data->phead.nbin;i++)
    {
      for (j=0;j<data->phead.nbin;j++)
    time_y[i*data->phead.nbin+j] = 0;
    } 
  for (l=sub0;l<data->phead.nsub;l++)
    {
      for (j=0;j<data->phead.npol && j < 2;j++) // Do not add cross terms!
    {
      for (i=0;i<data->phead.nchan;i++)
        {
          // Must calculate the frequency of this channel
          // ... calculate the delay caused by the DM
          // ... dedisperse the subintegration
          
          //          tdelay = 4.15e-3*dm*(pow(f0/1000.0,-2)-pow((f0+chanbw*i)/1000.0,-2));
          tdelay = 4.15e-3*dm*(pow(f0/1000.0,-2)-pow((f0+(chanbw*i-chanbw*data->phead.nchan/2.0))/1000.0,-2));
          cdelay = nint(-tdelay/bintime);
          //          if (l==0 && j==0)
          //        printf("Have %g %g %g %d %d\n",dm,f0,f0+chanbw*i,i,cdelay);
          fits_read_col_flt(fp,colnum,l+1,j*(data->phead.nchan*data->phead.nbin)+i*data->phead.nbin+1,data->phead.nbin,nval,ty,&initflag,&status);
          meanVal=0;
          //for (k=0;k<data->phead.nbin;k++)
        //{
         // ty[k] = ((ty[k]+offs[l][j*data->phead.nchan+i])*dat_scl[l][j*data->phead.nchan+i]); //+offs[l][j*data->phead.nchan+i]);

          /*          if (j==0)
            ty[k] -= bpass[i];
          else if (j==1)
            ty[k] -= bpass[data->phead.nchan+i];
          else
            {
              ty[k] -= sqrt(bpass[data->phead.nchan+i]*bpass[i]);
              } */
          //meanVal+=ty[k];
    //    }

          getbaseline(ty,data->phead.nbin,0.35,&meanVal,&rmsVal);
          if (rmsVal == 0.0 ) rmsVal=1.0;
          //          printf("Val = %g\n",ty[10]);
          for (k=0;k<data->phead.nbin;k++)
        {
          //          if (i==10 && l==10)
          //            printf("Orig value = %g\n",ty[k]);

          //ty[k] = ((ty[k]+offs[l][j*data->phead.nchan+i])*dat_scl[l][j*data->phead.nchan+i]); //+offs[l][j*data->phead.nchan+i]);
          ty[k] -= meanVal;
          ty[k] /= rmsVal;
          // Subtract bandpass
          /*          if (j==0)
            ty[k] -= bpass[i];
          else if (j==1)
            ty[k] -= bpass[data->phead.nchan+i];
          else
            {
              ty[k] -= sqrt(bpass[data->phead.nchan+i]*bpass[i]);
              } */
          //          if (l==10)
              
          //            printf("New value = %g\n",ty[k]);
          bn = k-cdelay + addDelay;
          //          bn = nint(fmod(k-tdelay/bintime,data->phead.nbin));
          while (bn >= data->phead.nbin)
            bn -= data->phead.nbin;
          while (bn < 0)
            bn += data->phead.nbin;
          freq_y[i*data->phead.nbin+k]+=(ty[k]); ///(float)(data->phead.npol*data->phead.nsub));
          time_y[l*data->phead.nbin+bn]+=(ty[k]);///(float)(data->phead.npol*data->phead.nchan));
          //          printf("timey = %g\n",time_y[l*data->phead.nbin+bn]);
          fy[bn]+=(ty[k]/(float)(data->phead.nchan*data->phead.npol*data->phead.nsub));
        }
        }
    }
    }
   
  // REMOVED: OFFS/DAT_SCL table no longer needed.
  //for (i=0;i<data->phead.nsub;i++)
  //  {
  //    free(offs[i]);
  //    free(dat_scl[i]);
  //  }
  //free(offs);
  //free(dat_scl);*/
#ifdef PFITS_DEBUG
  printf("Status = %d\n",status);
#endif
}

int extractPolData(fitsfile *fp,dSet *data,int pol,float *arr,float t1,float t2)
{
  int nsblk;
  long s1,s2,sub_1,samp_1,sub_2,samp_2,s,r1,r2;
  int nchan = data->phead.nchan;
  int nbits = data->phead.nbits;
  int npol  = data->phead.npol;
  int samplesperbyte = 8/data->phead.nbits;
  float chanVal[nchan];
  long ipos=0;
  int status=0;
  int colnum;
  int initflag=0;
  unsigned char *cval;
  unsigned char nval = '0';
  int i,c,count=0,sa;
  unsigned char chVals[nchan];
  int smoothSamp=1;
  int sm;
  int t=0;
  float tempVal;
  long pos=0,p,cc=0;
#ifdef PFITS_DEBUG
  printf("Have t1/t2 %g %g %g\n",t1,t2,(double)data->phead.tsamp);
#endif
  s1 = (long)(t1/data->phead.tsamp);
  s2 = (long)(t2/data->phead.tsamp);

#ifdef PFITS_DEBUG
  printf("samples per byte = %d\n",samplesperbyte);
  printf("Searching for %g %g\n",(double)s1,(double)s2);
#endif

  nsblk = data->phead.nsblk;
  if (s1 < 0) s1=0;
  findPosSample(data,s1,&sub_1,&samp_1);
  findPosSample(data,s2,&sub_2,&samp_2);
  if (sub_1 < 0) sub_1=0;

#ifdef PFITS_DEBUG
  printf("Have samples (%d,%d) (%d,%d)\n",sub_1,samp_1,sub_2,samp_2);
#endif
  fits_movnam_hdu(fp,BINARY_TBL,"SUBINT",1,&status);
  if (status) {
    printf("Unable to move to subint table in FITS file\n");
    exit(1);
  }
  fits_get_colnum(fp,CASEINSEN,"DATA",&colnum,&status);
  if (status) {
    printf("Unable to find data in the subint table in FITS file\n");
    exit(1);
  }
  cval = (unsigned char *)malloc(sizeof(unsigned char)*nsblk*nchan/samplesperbyte*npol);
  t=0;
#ifdef PFITS_DEBUG
  printf("Have sub1/2 = %g %g\n",(double)sub_1,(double)sub_2);
#endif
  for (s=sub_1;s<=sub_2;s++)
  {
      //            printf("Loading sub: %d\n",s);
      if (s==sub_1)
    r1 = samp_1;
      else
    r1 = 0;
      if (s==sub_2)
    r2 = samp_2;
      else
    r2 = nsblk;
#ifdef PFITS_DEBUG
    printf("Reading from %d to %d for subint %d\n",r1,r2,s);
#endif

      // Read nchan data points
      for (sa=r1;sa<r2;sa++)
    {
      fits_read_col_byt(fp,colnum,s+1,sa*nchan/samplesperbyte+1,nchan/samplesperbyte,nval,cval,&initflag,&status);
      if (status)
        {
          fits_report_error(stderr,status);
          exit(1);
        }
      bytesToValues(samplesperbyte, nchan, cval, chVals);
      tempVal = 0.0;
      //      val[count] = 0.0;
      for (i=0;i<nchan;i++)
        {
          //          printf("Here with %g\n",(float)chVals[i]);
                    p = cc+(s2-s1)*i;
          //          printf("pos-p = %d\n",p);
          //                    p=pos;
          if (samplesperbyte==8)
        {
          if (chVals[i]==0) arr[p] = 0.5;
          else arr[p] = -0.5;
        }
          else if (samplesperbyte==4)
        {
          if (chVals[i]==0) arr[p]=-2.5;
          else if (chVals[i]==1) arr[p] = -0.5;
          else if (chVals[i]==2) arr[p] = 0.5;
          else arr[p] = 2.5;
        }
          else
        {
          arr[p] = chVals[i];
          //          printf("Have got %g\n",arr[p]);
        }
          //          if (chVals[i]==0) arr[pos] = 0.5;
          //          else arr[pos] =- 0.5;
          //          pos++;
          pos++;

        }
      cc++;
      //      pos++;
    }
    }
#ifdef PFITS_DEBUG
    printf("Complete %d %d %d\n",pos,s2-s1,nchan);
#endif
  free(cval);
  return cc;
  //return pos;
}

int extractDataZeroDM(fitsfile *fp,dSet *data,long s1,long s2,float *mean,float *min,float *max,long maxVal,int *nsmooth)
{
  long sub_1,samp_1;
  long sub_2,samp_2;
  int  s,r1,r2;
  int nsblk;
  int nchan = data->phead.nchan;
  int nbits = data->phead.nbits;
  int npol  = data->phead.npol;
  int samplesperbyte = 8/data->phead.nbits;
  float chanVal[nchan];
  long ipos=0;
  int status=0;
  int colnum;
  int initflag=0;
  unsigned char *cval;
  unsigned char nval = '0';
  int i,c,count=0,sa;
  unsigned char chVals[nchan];
  int smoothSamp=1;
  int sm;
  int t=0;
  float tempVal;

#ifdef PFITS_DEBUG
  printf("Extract zeroDM smaplesperbyte = %d\n",samplesperbyte);
#endif
  // Do we require smoothing?
  if (s2 - s1 > maxVal)
  {
    smoothSamp = ceil((double)(s2-s1)/(double)maxVal);
  }

  nsblk = data->phead.nsblk;

  findPosSample(data,s1,&sub_1,&samp_1);
  findPosSample(data,s2,&sub_2,&samp_2);
  //  printf("Here with %d %d %d %d\n",sub_1,samp_1,sub_2,samp_2);
  // Go to the subint table

  fits_movnam_hdu(fp,BINARY_TBL,"SUBINT",1,&status);
  if (status) {
    printf("Unable to move to subint table in FITS file\n");
    exit(1);
  }
  fits_get_colnum(fp,CASEINSEN,"DATA",&colnum,&status);
  if (status) {
    printf("Unable to find data in the subint table in FITS file\n");
    exit(1);
  }
  cval = (unsigned char *)malloc(sizeof(unsigned char)*nsblk*nchan/samplesperbyte*npol);
  t=0;
  for (s=sub_1;s<=sub_2;s++)
  {
    if (s==sub_1)
    r1 = samp_1;
      else
    r1 = 0;
      if (s==sub_2)
    r2 = samp_2;
      else
    r2 = nsblk;

    // Read nchan data points
    for (sa=r1;sa<r2;sa++)
    {
      fits_read_col_byt(fp,colnum,s+1,sa*nchan/samplesperbyte+1,nchan/samplesperbyte,nval,cval,&initflag,&status);
      if (status)
      {
        fits_report_error(stderr,status);
        exit(1);
      }
      bytesToValues(samplesperbyte, nchan, cval, chVals);
      tempVal = 0.0;
      //      val[count] = 0.0;
      for (i=0;i<nchan;i++)
      {
        if (samplesperbyte==8)
        {
          if (chVals[i]==0) tempVal += 0.5;
          else tempVal -= 0.5;
        }
        else
        {
          tempVal+=chVals[i];
        }
      }
      tempVal/=(float)nchan;
      if (t==0)
      {
        mean[count] = min[count] = max[count] = tempVal;
      }
      else
      {
        if (min[count] > tempVal)
          min[count] = tempVal;
        if (max[count] < tempVal)
          max[count] = tempVal;
        mean[count] += tempVal;
      }
      t++;
      if (t==smoothSamp)
      {
        mean[count]/=(double)smoothSamp;
        count++;
        t=0;
      }
    }
  }
  //  printf("Complete %d %d\n",t,count);
  free(cval);
  *nsmooth = smoothSamp;
  if (count > maxVal)
  {
    printf("In pfits.c -- this should not have happened\n");
    printf("count = %d, maxVal = %d\n",(int)count,(int)maxVal);
    printf("s2-s1 = %d\n",(int)(s2-s1));
    printf("Exit ...\n");
    exit(1);
  }
  //  printf("Returning\n");
  return count;
}

int extractData(fitsfile *fp,dSet *data,long s1,long s2,float *mean,float *min,float *max,long maxVal,int *nsmooth,float dm)
{
  long sub_1,samp_1;
  long sub_2,samp_2;
  int  s,r1,r2;
  int nsblk;
  int nchan = data->phead.nchan;
  int nbits = data->phead.nbits;
  int npol  = data->phead.npol;
  int samplesperbyte = 8/data->phead.nbits;
  float chanVal[nchan];
  long ipos=0;
  int status=0;
  int colnum;
  int initflag=0;
  unsigned char *cval;
  unsigned char nval = '0';
  int i,c,count=0,sa;
  unsigned char chVals[nchan];
  int smoothSamp=1;
  int sm;
  unsigned int t=0;
  unsigned int j;
  float tempVal;
  int nsampDM;
  float tDM;
  float bw = data->phead.bw; // MHz
  float f0 = data->phead.freq+fabs(data->phead.bw)/2.0; // MHz // Highest frequency
  float chanbw = data->phead.bw/data->phead.nchan;
  unsigned char *ring_s;
  unsigned char *ring_e;
  unsigned char *ring_pos;
  unsigned char *ring_last;
  unsigned int bitCount;
  int cdelay,k,l;
  float tdelay;
  unsigned int pos;
  unsigned int **cde_byte;
  unsigned char **cde_bit;
  unsigned int nchanbyte;

  nchanbyte = nchan/samplesperbyte;

  // Do we require smoothing?
  printf("At here: smoothSamp = %d\n",smoothSamp);
  if (s2 - s1 > maxVal)
    {
      smoothSamp = ceil((double)(s2-s1)/(double)maxVal);
      printf("aaa: smoothSamp = %d\n",smoothSamp);
    }
  printf("2 At here: smoothSamp = %d\n",smoothSamp);
  nsblk = data->phead.nsblk;

  findPosSample(data,s1,&sub_1,&samp_1);
  findPosSample(data,s2,&sub_2,&samp_2);
  //  printf("Here with %d %d %d %d\n",sub_1,samp_1,sub_2,samp_2);
  // Go to the subint table

  fits_movnam_hdu(fp,BINARY_TBL,"SUBINT",1,&status);
  if (status) {
    printf("Unable to move to subint table in FITS file\n");
    exit(1);
  }
  fits_get_colnum(fp,CASEINSEN,"DATA",&colnum,&status);
  if (status) {
    printf("Unable to find data in the subint table in FITS file\n");
    exit(1);
  }
  // Number of extra samples required to process 
  // to de-disperse across the band
  tDM = fabs(4.15e-3*dm*(pow((f0)/1000.0,-2)-pow((f0-fabs(bw))/1000.0,-2)));
  nsampDM = ceil(tDM/data->phead.tsamp);

  printf("Here with tDM = %g %d %d\n",tDM,nsampDM,smoothSamp);

  cval = (unsigned char *)malloc(sizeof(unsigned char)*nchanbyte*npol*nsampDM);

  printf("Trying to allocate memory\n");
  cde_byte = (unsigned int **)malloc(sizeof(unsigned int *)*nsampDM);
  cde_bit  = (unsigned char **)malloc(sizeof(unsigned char *)*nsampDM);

  for (i=0;i<nsampDM;i++)
    {
      cde_byte[i] = (unsigned int *)malloc(sizeof(unsigned int)*nchan*npol*nsampDM);
      cde_bit[i] = (unsigned char *)malloc(sizeof(unsigned char)*nchan*npol*nsampDM);
    }
  printf("Allocated memory\n");

  // Now turn on the bits that correspond to the DM
  j=0;
  for (i=0;i<nchan;i++)
    {
      tdelay = 4.15e-3*dm*(pow(f0/1000.0,-2)-pow((f0-fabs(chanbw)*i)/1000.0,-2));
      cdelay = nint(-tdelay/data->phead.tsamp);
      //      printf("Have %g %d\n",tdelay,cdelay);
      for (j=0;j<nsampDM;j++)
    {
      k = (int)((float)i/(float)samplesperbyte)+cdelay*nchanbyte+j*nchan/samplesperbyte;
      if (k>(nsampDM)*nchanbyte)
        k-=nsampDM*nchanbyte;
      if (k>(nsampDM)*nchanbyte)
        {
          printf("ERROR\n");
          exit(1);
        }
      cde_byte[j][i] = k;
      cde_bit[j][i]  = i%8;
    }
    }

 printf("Ready\n");

  ring_s   = cval;
  ring_pos = cval;
  ring_e   = &cval[nchanbyte*npol*(nsampDM)];
     printf("Got here 2\n");
  // Should check if I'm running off the end of a subint
  s = sub_1;
  sa = samp_1;
  // Use a ring buffer to store the data
  // Fill it up
  count = 0;
  t=0;
  pos = 0;
  printf("Got here 3\n");
  for (i=0;i<s2-s1-nsampDM;i++)
    {
      //            printf("i = %d, s = %d, sa = %d\n",i,s,sa);
      fits_read_col_byt(fp,colnum,s+1,(sa)*nchanbyte+1,nchanbyte,nval,ring_pos,&initflag,&status);
      sa++;
      
      if (i >= nsampDM-1)
    {
      bitCount = 0;
      //      for (k=0;k<35;k++)
        {
          for (j=0;j<nchan;j++)
        {
          //          bitCount += extractBit(cval[cde_byte[t][j]],cde_bit[t][j]);
          // ONLY 1-BIT DATA
          bitCount += extractBit(cval[cde_byte[t][j]],j%8);
          
          //          bitCount+=cde_byte[t][j]+cde_bit[t][j];
        }
        }
      if (pos==0)
        {
          min[count]=max[count]=mean[count] = bitCount; 
        }
      else
        {
          if (min[count] > bitCount) min[count] = bitCount;
          else if (max[count] < bitCount) max[count] = bitCount;
          mean[count] += bitCount;
        }
      t++;
      pos++;
      if (pos==smoothSamp)
        {
          mean[count]/=(double)smoothSamp;
          count++;
          pos=0;
        } 
    }

      ring_pos += nchanbyte*npol;
      if (ring_pos == ring_e)
    {
      //      printf("Got to end of ring buffer %d %d\n",i,nsampDM);
      ring_pos = ring_s;
      t=0;
    } 
      if (sa == nsblk)
    {
      //      printf("Running off the end of a block\n");
      s++;
      sa = 0;
    } 
    }
  
  tempVal=0;
  //  printf("Got here\n");
  //  exit(1);
  if (status)
    {
      fits_report_error(stderr,status);
      exit(1);
    }


  //  printf("Complete %d %d\n",t,count);
  free(cval);
  for (i=0;i<nsampDM;i++)
    {
      free(cde_byte[i]);
      free(cde_bit[i]);
    }
  free(cde_byte);
  free(cde_bit);

  *nsmooth = smoothSamp;
  if (count > maxVal)
    {
      printf("In pfits.c -- this should not have happened\n");
      printf("count = %d, maxVal = %d\n",(int)count,(int)maxVal);
      printf("s2-s1 = %d\n",(int)(s2-s1));
      printf("Exit ...\n");
      exit(1);
    }
  //  printf("Returning\n");
  return count;
}


int writeData(fitsfile *fp,FILE *fout,dSet *data,long s1,long s2,float dm)
{
  long sub_1,samp_1;
  long sub_2,samp_2;
  int  s,r1,r2;
  int nsblk;
  int nchan = data->phead.nchan;
  int nbits = data->phead.nbits;
  int npol  = data->phead.npol;
  int samplesperbyte = 8/data->phead.nbits;
  float chanVal[nchan];
  long ipos=0;
  int status=0;
  int colnum;
  int initflag=0;
  unsigned char *cval;
  unsigned char nval = '0';
  int i,c,count=0,sa;
  unsigned char chVals[nchan];
  int smoothSamp=1;
  int sm;
  unsigned int t=0;
  unsigned int j;
  float tempVal;
  int nsampDM;
  float tDM;
  float bw = data->phead.bw; // MHz
  float f0 = data->phead.freq+fabs(data->phead.bw)/2.0; // MHz // Highest frequency
  float chanbw = data->phead.bw/data->phead.nchan;
  unsigned char *ring_s;
  unsigned char *ring_e;
  unsigned char *ring_pos;
  unsigned char *ring_last;
  unsigned int bitCount;
  int cdelay,k,l;
  float tdelay;
  unsigned int pos;
  unsigned int **cde_byte;
  unsigned char **cde_bit;
  unsigned int nchanbyte;

  nchanbyte = nchan/samplesperbyte;

  // Do we require smoothing?
  printf("At here: smoothSamp = %d\n",smoothSamp);
  nsblk = data->phead.nsblk;

  findPosSample(data,s1,&sub_1,&samp_1);
  findPosSample(data,s2,&sub_2,&samp_2);
  //  printf("Here with %d %d %d %d\n",sub_1,samp_1,sub_2,samp_2);
  // Go to the subint table

  fits_movnam_hdu(fp,BINARY_TBL,"SUBINT",1,&status);
  if (status) {
    printf("Unable to move to subint table in FITS file\n");
    exit(1);
  }
  fits_get_colnum(fp,CASEINSEN,"DATA",&colnum,&status);
  if (status) {
    printf("Unable to find data in the subint table in FITS file\n");
    exit(1);
  }
  // Number of extra samples required to process 
  // to de-disperse across the band
  tDM = fabs(4.15e-3*dm*(pow((f0)/1000.0,-2)-pow((f0-fabs(bw))/1000.0,-2)));
  nsampDM = ceil(tDM/data->phead.tsamp);


  cval = (unsigned char *)malloc(sizeof(unsigned char)*nchanbyte*npol*nsampDM);

  printf("Trying to allocate memory\n");
  cde_byte = (unsigned int **)malloc(sizeof(unsigned int *)*nsampDM);
  cde_bit  = (unsigned char **)malloc(sizeof(unsigned char *)*nsampDM);

  for (i=0;i<nsampDM;i++)
    {
      cde_byte[i] = (unsigned int *)malloc(sizeof(unsigned int)*nchan*npol*nsampDM);
      cde_bit[i] = (unsigned char *)malloc(sizeof(unsigned char)*nchan*npol*nsampDM);
    }
  printf("Allocated memory\n");

  // Now turn on the bits that correspond to the DM
  j=0;
  for (i=0;i<nchan;i++)
    {
      tdelay = 4.15e-3*dm*(pow(f0/1000.0,-2)-pow((f0-fabs(chanbw)*i)/1000.0,-2));
      cdelay = nint(-tdelay/data->phead.tsamp);
      //      printf("Have %g %d\n",tdelay,cdelay);
      for (j=0;j<nsampDM;j++)
    {
      k = (int)((float)i/(float)samplesperbyte)+cdelay*nchanbyte+j*nchan/samplesperbyte;
      if (k>(nsampDM)*nchanbyte)
        k-=nsampDM*nchanbyte;
      if (k>(nsampDM)*nchanbyte)
        {
          printf("ERROR\n");
          exit(1);
        }
      cde_byte[j][i] = k;
      cde_bit[j][i]  = i%8;
    }
    }

 printf("Ready\n");

  ring_s   = cval;
  ring_pos = cval;
  ring_e   = &cval[nchanbyte*npol*(nsampDM)];
     printf("Got here 2\n");
  // Should check if I'm running off the end of a subint
  s = sub_1;
  sa = samp_1;
  // Use a ring buffer to store the data
  // Fill it up
  count = 0;
  t=0;
  pos = 0;
  printf("Got here 3\n");
  fprintf(fout,"# %s %d %g %g\n",data->phead.source,data->phead.imjd,data->phead.freq,data->phead.tsamp);
  for (i=0;i<s2-s1-nsampDM;i++)
    {
      //            printf("i = %d, s = %d, sa = %d\n",i,s,sa);
      fits_read_col_byt(fp,colnum,s+1,(sa)*nchanbyte+1,nchanbyte,nval,ring_pos,&initflag,&status);
      sa++;
      
      if (i >= nsampDM-1)
    {
      bitCount = 0;
      //      for (k=0;k<35;k++)
        {
          for (j=0;j<nchan;j++)
        {
          //          bitCount += extractBit(cval[cde_byte[t][j]],cde_bit[t][j]);
          // ONLY 1-BIT DATA
          bitCount += extractBit(cval[cde_byte[t][j]],j%8);
          
          //          bitCount+=cde_byte[t][j]+cde_bit[t][j];
        }
        }
        fprintf(fout,"%d\n",bitCount);
        t++;
    }
      
      ring_pos += nchanbyte*npol;
      if (ring_pos == ring_e)
    {
      //      printf("Got to end of ring buffer %d %d\n",i,nsampDM);
      ring_pos = ring_s;
      t=0;
    } 
      if (sa == nsblk)
     {
       //      printf("Running off the end of a block\n");
       s++;
       sa = 0;
     } 
     }
  
  tempVal=0;
  //  printf("Got here\n");
  //  exit(1);
  if (status)
    {
      fits_report_error(stderr,status);
      exit(1);
    }
  
  
  //  printf("Complete %d %d\n",t,count);
  free(cval);
   for (i=0;i<nsampDM;i++)
     {
       free(cde_byte[i]);
       free(cde_bit[i]);
     }
   free(cde_byte);
   free(cde_bit);

   return count;
}



// Return the folded, dedispersed profile for a grid of DM values 
// and fold-periods
//
/*int foldDM(fitsfile *fp,dSet *data,long s1,long s2,float ***fold,float *dm,double *period,int ndm,int nperiod,int nbin)
{
  long sub_1,samp_1;
  long sub_2,samp_2;
  int  s,r1,r2;
  int nsblk;
  int nchan = data->phead.nchan;
  int nbits = data->phead.nbits;
  int npol  = data->phead.npol;
  int samplesperbyte = 8/data->phead.nbits;
  float chanVal[nchan];
  long ipos=0;
  int status=0;
  int colnum;
  int initflag=0;
  unsigned char *cval;
  unsigned char nval = '0';
  int i,c,count=0,sa;
  unsigned char chVals[nchan];
  int smoothSamp=1;
  int sm;
  int t=0;
  unsigned int j;
  float tempVal;
  int nsampDM;
  float tDM;
  float bw = data->phead.bw; // MHz
  float f0 = data->phead.freq+fabs(data->phead.bw)/2.0; // MHz // Highest frequency
  float chanbw = data->phead.bw/data->phead.nchan;
  unsigned char *ring_s;
  unsigned char *ring_e;
  unsigned char *ring_pos;
  unsigned char *ring_last;
  unsigned int bitCount;
  int cdelay,k,l;
  float tdelay;
  int pos;
  unsigned int **cde_byte;
  unsigned char **cde_bit;
  unsigned int nchanbyte;

  nchanbyte = nchan/samplesperbyte;

  // Do we require smoothing?
  printf("At here: smoothSamp = %d\n",smoothSamp);
  if (s2 - s1 > maxVal)
    {
      smoothSamp = ceil((double)(s2-s1)/(double)maxVal);
      printf("aaa: smoothSamp = %d\n",smoothSamp);
    }
  printf("2 At here: smoothSamp = %d\n",smoothSamp);
  nsblk = data->phead.nsblk;

  findPosSample(data,s1,&sub_1,&samp_1);
  findPosSample(data,s2,&sub_2,&samp_2);
  printf("Here with %d %d %d %d\n",sub_1,samp_1,sub_2,samp_2);
  // Go to the subint table

  fits_movnam_hdu(fp,BINARY_TBL,"SUBINT",1,&status);
  if (status) {
    printf("Unable to move to subint table in FITS file\n");
    exit(1);
  }
  fits_get_colnum(fp,CASEINSEN,"DATA",&colnum,&status);
  if (status) {
    printf("Unable to find data in the subint table in FITS file\n");
    exit(1);
  }
  // Number of extra samples required to process 
  // to de-disperse across the band
  tDM = fabs(4.15e-3*dm*(pow((f0)/1000.0,-2)-pow((f0-fabs(bw))/1000.0,-2)));
  nsampDM = ceil(tDM/data->phead.tsamp);

  printf("Here with tDM = %g %d %d\n",tDM,nsampDM,smoothSamp);

  cval = (unsigned char *)malloc(sizeof(unsigned char)*nchanbyte*npol*nsampDM);

  printf("Trying to allocate memory\n");
  cde_byte = (unsigned int **)malloc(sizeof(unsigned int *)*nsampDM);
  cde_bit  = (unsigned char **)malloc(sizeof(unsigned char *)*nsampDM);

  for (i=0;i<nsampDM;i++)
    {
      cde_byte[i] = (unsigned int *)malloc(sizeof(unsigned int)*nchan*npol*nsampDM);
      cde_bit[i] = (unsigned char *)malloc(sizeof(unsigned char)*nchan*npol*nsampDM);
    }
  printf("Allocated memory\n");

  // Now turn on the bits that correspond to the DM
  j=0;
  for (i=0;i<nchan;i++)
    {
      tdelay = 4.15e-3*dm*(pow(f0/1000.0,-2)-pow((f0-fabs(chanbw)*i)/1000.0,-2));
      cdelay = nint(-tdelay/data->phead.tsamp);
      //      printf("Have %g %d\n",tdelay,cdelay);
      for (j=0;j<nsampDM;j++)
    {
      k = (int)((float)i/(float)samplesperbyte)+cdelay*nchanbyte+j*nchan/samplesperbyte;
      if (k>(nsampDM)*nchanbyte)
        k-=nsampDM*nchanbyte;
      if (k>(nsampDM)*nchanbyte)
        {
          printf("ERROR\n");
          exit(1);
        }
      cde_byte[j][i] = k;
      cde_bit[j][i]  = i%8;
    }
    }

 printf("Ready\n");
  ring_s   = cval;
  ring_pos = cval;
  ring_e   = &cval[nchanbyte*npol*(nsampDM)];
     printf("Got here 2\n");
  // Should check if I'm running off the end of a subint
  s = sub_1;
  sa = samp_1;
  // Use a ring buffer to store the data
  // Fill it up
  count = 0;
  t=0;
  pos = 0;
  printf("Got here 3\n");
  for (i=0;i<s2-s1-nsampDM;i++)
    {
      //            printf("i = %d, s = %d, sa = %d\n",i,s,sa);
      fits_read_col_byt(fp,colnum,s+1,(sa)*nchanbyte+1,nchanbyte,nval,ring_pos,&initflag,&status);
      sa++;
      
      if (i >= nsampDM-1)
    {
      bitCount = 0;
      for (j=0;j<nchan;j++)
        bitCount += extractBit(cval[cde_byte[t][j]],cde_bit[t][j]);
      if (pos==0)
        {
          min[count]=max[count]=mean[count] = bitCount; 
        }
      else
        {
          if (min[count] > bitCount) min[count] = bitCount;
          else if (max[count] < bitCount) max[count] = bitCount;
          mean[count] += bitCount;
      
        }

      t++;
      pos++;
      if (pos==smoothSamp)
        {
          mean[count]/=(double)smoothSamp;
          count++;
          pos=0;
        } 
    }

      ring_pos += nchanbyte*npol;
      if (ring_pos == ring_e)
    {
      //      printf("Got to end of ring buffer %d %d\n",i,nsampDM);
      ring_pos = ring_s;
      t=0;
    } 
      if (sa == nsblk)
    {
      //      printf("Running off the end of a block\n");
      s++;
      sa = 0;
    } 
    }
  
  tempVal=0;
  //  printf("Got here\n");
  //  exit(1);
  if (status)
    {
      fits_report_error(stderr,status);
      exit(1);
    }


  //  printf("Complete %d %d\n",t,count);
  free(cval);
  for (i=0;i<nsampDM;i++)
    {
      free(cde_byte[i]);
      free(cde_bit[i]);
    }
  free(cde_byte);
  free(cde_bit);

  *nsmooth = smoothSamp;
  if (count > maxVal)
    {
      printf("In pfits.c -- this should not have happened\n");
      printf("count = %d, maxVal = %d\n",count,maxVal);
      printf("s2-s1 = %d\n",s2-s1);
      printf("Exit ...\n");
      exit(1);
    }
  //  printf("Returning\n");
  return count;
} */



void findPosSample(dSet *data,long s,long *sub,long *samp)
{
  int nsblk;

  nsblk = data->phead.nsblk;
#ifdef PFITS_DEBUG
  printf("Searching for sample %d\n",(int)s);
#endif
  
  *sub  = floor((double)s/(double)nsblk);
  *samp = s-*sub*nsblk;
}

// Return the number of channels
int readBandpass(fitsfile *fp,float *bpass)
{
  int status=0;
  int samplesperbyte;
  long nrows;
  int colnum;
  int j,pol,subint;
  int initflag=0;
  float val=1;
  float nval = 0;
  long npts;
  int start,end;
  int i,anynul;
  float dat_offs,dat_scl;
  float freq,bw,chanbw;
  long nchan;
  int naxis;
  int typecode;
  long repeat,width;

  // Go to bandpass table
  fits_movnam_hdu(fp,BINARY_TBL,"BANDPASS",1,&status);
  if (status) {
    printf("Unable to move to bandpass table in FITS file\n");
    exit(1);
  }

  fits_get_colnum(fp,CASEINSEN,"DAT_OFFS",&colnum,&status);
  fits_read_col(fp,TFLOAT,colnum,1,1,1,&nval,&dat_offs,&anynul,&status);
  if (status){
    printf("Error in bandpass header 1\n");
    fits_report_error(stderr,status);
    exit(1);
  }
  fits_get_colnum(fp,CASEINSEN,"DAT_SCL",&colnum,&status);
  fits_read_col(fp,TFLOAT,colnum,1,1,1,&nval,&dat_scl,&anynul,&status);
  if (status){
    printf("Error in bandpass header\n");
    fits_report_error(stderr,status);
    exit(1);
    }

  fits_get_colnum(fp,CASEINSEN,"DATA",&colnum,&status);
  if (status){
    fits_report_error(stderr,status);
    exit(1);
  }
  fits_get_coltype(fp,colnum,&typecode,&repeat,&width,&status);
  nchan = (int)repeat-1; // Don't read the first channel
  printf("nc = %d, na = %d\n",(int)nchan,(int)naxis);
  for (i=1;i<(int)nchan;i++)
    {
      fits_read_col(fp,TFLOAT,colnum,1,i+1,1,&nval,&val,&anynul,&status);
      bpass[i-1] = val*dat_scl + dat_offs;
      if (status){
    fits_report_error(stderr,status);
    exit(1);
      }
    }
  fits_close_file(fp,&status);  
  return nchan;
}


void bytesToValues(int samplesperbyte, int nchan, unsigned char *val2, 
           unsigned char *values)
{
  int i,j;
  int pos=0;
  int index = 0;  

  for (i = 0; i < nchan/samplesperbyte; i++) 
    {
      void (*pointer)(int, unsigned char[], int*);
      switch (samplesperbyte)
    {
    case 1:
      pointer = eightBits;
      break;
      
    case 2:
      pointer = fourBits;
      break;
      
    case 4:
      pointer = twoBits;
      break;
      
    case 8:
      pointer = oneBit;
      break;
    }
      //      printf("i = %d; Here with index = %d, val2 = %d\n",i,index,val2[i]);
      pointer(val2[i], values, &index);
    }
  //  for (i=0;i<nchan;i++)
  //    printf("vals = %g\n",(float)values[i]);
  //  printf("Finsihed\n");
}

void eightBits(int eight_bit_number, unsigned char *results, int *index)
{
    // 135 subtracted from the original number to make the range from -8 to something
    // 0.5 is then added to remove bias

    //results[(*index)++] = eight_bit_number - 135.5;
  results[(*index)++] = eight_bit_number; // -31.5; // -31.5;
}

void fourBits(int eight_bit_number, unsigned char *results, int *index)
{
    // anding the least significant 4 bits with 0000 1111 will give the (signed) 4-bit number
    // shifting right 4 bits will produce the other (first) 4-bit numbers
    // 0.5 is added to each number to compensate for the bias, making the range -7.5 -> +7.5

    float tempResults[2];
    int i;
    for (i = 0; i < 2; i++) {
        int andedNumber = eight_bit_number & 15;
        tempResults[i] = andedNumber; // - 7.5;
        eight_bit_number = eight_bit_number >> 4;
    }

    for (i = 1; i >= 0; i--)
        results[(*index)++] = tempResults[i];
}

void twoBits(int eight_bit_number, unsigned char *results, int *index)
{
    // anding the least significant 2 bits with 0000 0011 will give the (signed 2-bit number
    // shifting right 2 bits will produce the next (previous) 2-bit number
    // the numbers are adjusted to compensate for bias and to give a rms of ~0.9
    //
    // 1 -> 2.0
    // 0 -> +0.5
    // -1 -> -0.5
    // -2 -> -2.0

    float tempResults[4];
    int i;

    for (i = 0; i < 4; i++) {
        int andedNumber = eight_bit_number & 3;
        switch (andedNumber) {
            case 0:
          tempResults[i] = 0; //-2.5;
                break;
            case 1:
          tempResults[i] = 1; // -0.5;
                break;
            case 2:
          tempResults[i] = 2; //0.5;
                break;
            case 3:
          tempResults[i] = 3; // 2.5;
                break;
        }
        eight_bit_number = eight_bit_number >> 2;
    }

    for (i = 3; i >= 0; i--)
      {
        results[(*index)++] = tempResults[i];
    //    printf("ret: %d %g\n",(*index),tempResults[i]);
      }
}

void oneBit(int eight_bit_number, unsigned char *results, int *index)
{
    // anding the least significant bit with 0000 0001 will give the (signed 1-bit number
    // shifting right 1 bit will produce the next (previous) 1-bit number
    //
    // 0.5 is added to each number to compensate for bias
    //

    float tempResults[8];
    int i;

    for (i = 0; i < 8; i++) {
        int andedNumber = eight_bit_number & 1;
    // tempResults[i] = andedNumber ? 0.5 : -0.5;
    tempResults[i] = andedNumber ? 0 : 1;
        eight_bit_number = eight_bit_number >> 1;
    }

    for (i = 7; i >= 0; i--)
      {
    results[(*index)++] = tempResults[i];
    //    printf("This pt: %d %g\n",*index,tempResults[i]);
      }
}

int andCount(unsigned char *cval,unsigned char *mask,int n)
{
  unsigned char andRes;
  unsigned int i,nbit=0;
  unsigned char c=0;
  for (i=0;i<n;i++)
    {
      //      printf("Comparing %d %d\n",(int)cval[i],(int)mask[i]);

      andRes = cval[i] & mask[i];
      for (c=0;andRes;c++)
          andRes &= andRes-1;
      nbit+=c;
    }
  // Now count the bits
  return nbit;
}

int nint(float val)
{
  if (val >= 0)
    return (int)(val+0.5);
  return -(int)(-val+0.5);
}

void setMask(int maskN,int maskP,unsigned char **mask)
{
  int mpos1,mpos2;
  unsigned char temp;
  // Must set the correct bit to '1' in the mask
  mpos1 = (int)(maskP/8.0);
  mpos2 = (int)(maskP-mpos1*8);
  if (mpos2==0) temp = 1;
  else if (mpos2==1) temp = 2;
  else if (mpos2==2) temp = 4;
  else if (mpos2==3) temp = 8;
  else if (mpos2==4) temp = 16;
  else if (mpos2==5) temp = 32;
  else if (mpos2==6) temp = 64;
  else if (mpos2==7) temp = 128;
  mask[maskN][mpos1] = temp ^ mask[maskN][mpos1];
}


unsigned int extractBit(unsigned char byte,unsigned int pos)
{
  return (byte >> pos) & 0x01;
}

void getbaseline(float* prof, int nbin, float fwindow,float* mean, float *rms){
    int i;
    int window=fwindow*nbin;
    float sum=0;
    float minsum=0;
    for (i=0;i<window;i++){
        sum+=prof[i];
    }
    minsum=sum;
    float totalsum=0;
    int offbin=0;
    for (i=0;i<nbin;i++){
        if (sum < minsum){
            offbin=i;
            minsum=sum;
        }
        sum-=prof[i];
        sum+=prof[(window+i)%nbin];
    }
    (*mean)=minsum/(float)window;
    sum=0;
    for (i=0;i<window;i++){
        float v = prof[(i+offbin)%nbin]-(*mean);
        sum+=v*v;
    }
    (*rms)=sqrt(sum/(float)window);


    return;

}

//
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "pfits.h"
#include <cpgplot.h>

#define MAX_RET 4096

void plotFrequency(fitsfile *fp,dSet *d,float mint,float maxt);
void plotBpass(float *arr,int nchan,int npts);
void draw_histogram(float *x,int count,int nbin,float minx,float maxx,int normalise,float maxy,int colour,int log,float offset,int histOutput);

void help()
{
  printf("-dm <dm>      set the dispersion measure\n");
  printf("-f filename   PSRFITS file name\n");
  printf("-g gr         set graphics device\n");
  printf("-t1 <start>   start time in seconds from beginning of observation\n");
  printf("-t2 <end>     end time in seconds from beginning of observation\n");
  printf("\n\n\n");
}

int main(int argc,char *argv[])
{
  int i;
  char fname[128];
  dSet *data; 
  float freq,bw,chanbw;
  int nchan,npol;
  float mx,my,mx2,my2;
  float binw;
  float mean[MAX_RET];
  float min[MAX_RET];
  float max[MAX_RET];

  float omean[MAX_RET]; // Original values
  float omin[MAX_RET];
  float omax[MAX_RET];
  long onp,osamp1,osamp2;

  float fx[MAX_RET];
  float dm=0;
  int np;
  long samp1 = -1;
  long samp2 = -1;
  float minx,maxx,miny=-1,maxy=260;
  int nsmooth,onsmooth;
  int findXrange=1;
  int findYrange=1;
  int plotSubint=-1;
  float ts=-1,te=-1;
  char key;
  char grDev[128]="1/xs";
  int interactive=1;
  int stopit=0;
  fitsfile *fp;
  float tx[2],ty[2];

  help();
  


  data = initialiseDset();

  for (i=0;i<argc;i++)
    {
      if (strcmp(argv[i],"-f")==0)
	strcpy(fname,argv[++i]);
      else if (strcmp(argv[i],"-g")==0)
	{
	  strcpy(grDev,argv[++i]);
	  interactive=0;
	}
      else if (strcmp(argv[i],"-t1")==0)
	sscanf(argv[++i],"%f",&ts);
      else if (strcmp(argv[i],"-t2")==0)
	sscanf(argv[++i],"%f",&te);
      else if (strcmp(argv[i],"-dm")==0)
	sscanf(argv[++i],"%f",&dm);
      else if (strcmp(argv[i],"-stop")==0)
	stopit=1;
    }
  
  fp   = openFitsFile(fname); 
  loadPrimaryHeader(fp,data);
  displayHeaderInfo(data);
  if (ts!=-1)
    samp1 = ts/data->phead.tsamp;
  if (te!=-1)
    samp2 = te/data->phead.tsamp;

  if (samp1==-1) samp1=0;
  if (samp2==-1) samp2=data->phead.nsub*data->phead.nsblk; 

  if (dm==0)
    np = extractDataZeroDM(fp,data,samp1,samp2,mean,min,max,MAX_RET,&nsmooth);
  else
    np = extractData(fp,data,samp1,samp2,mean,min,max,MAX_RET,&nsmooth,dm);

  if (stopit==1)
    exit(1);
  cpgbeg(0,grDev,1,1);
  cpgask(0);
  //    exit(1);
  printf("np = %d\n",np);
  for (i=0;i<np;i++)
    {
      omean[i] = mean[i];
      omin[i]  = min[i];
      omax[i]  = max[i];
    }
  onsmooth = nsmooth;
  onp = np;
  osamp1 = samp1;
  osamp2 = samp2;


  printf("\n\n");
  printf("f       show a frequency-time plot for the zoom region\n");
  printf("l       ascii listing of dispersed data in zoom region\n");
  printf("q       quit\n");
  printf("r       reload data within the zoom region\n");
  printf("s       show subintegration boundaries\n");
  printf("u       unzoom to specified region given the last reload\n");
  printf("U       show the entire observation\n");
  printf("w       create .wav file for zoom region\n");
  printf("z       zoom into a specified region\n");
  printf("\n\n");
  do {
    for (i=0;i<np;i++)
      {
	fx[i] = samp1*data->phead.tsamp + (i+0.5)*data->phead.tsamp*nsmooth ;
	if (findXrange==1)
	  {
	    if (i==0)
	      minx = maxx = fx[i];
	    else if (minx > fx[i]) minx = fx[i];
	    else if (maxx < fx[i]) maxx = fx[i];
	  }
	if (findYrange==1)
	  {
	    if (i==0)
	      {
		miny = min[i];
		maxy = max[i];
	      }
	    if (miny > min[i]) miny = min[i];
	    if (maxy < max[i]) maxy = max[i];
	  }
	
      }
    findXrange=0;
    findYrange=0;
    cpgenv(minx,maxx,miny,maxy,0,1);
    cpglab("Time (s)","",fname);
    printf("nsmooth = %d\n",nsmooth);
    if (nsmooth > 1)
      {
	cpgsci(2); cpgbin(np,fx,min,1);
	cpgsci(3); cpgbin(np,fx,max,1);
      }
    cpgsci(1); cpgbin(np,fx,mean,1);
    //    cpgpt(np,fx,mean,4);
    if (plotSubint==1)
      {
	for (i=0;i<data->phead.nsub;i++)
	  {
	    tx[0] = tx[1] = i*data->phead.nsblk*data->phead.tsamp;
	    ty[0] = miny; ty[1] = maxy;
	    cpgsci(7); cpgsls(4); cpgline(2,tx,ty); cpgsci(1); cpgsls(1);
	  }
      }
    cpgcurs(&mx,&my,&key);
    if (key=='z')
      {
	cpgband(2,0,mx,my,&mx2,&my2,&key);
	if (mx > mx2) {maxx = mx; minx = mx2;}
	else {maxx = mx2; minx = mx;}
	
	if (my > my2) {maxy = my; miny = my2;}
	else {maxy = my2; miny = my;}	   
	findXrange=0;
	findYrange=0;
      }
    else if (key=='l')
      {
	char fname[128];
	FILE *fout;
	long s1,s2;

	printf("Enter filename: ");
	scanf("%s",fname);
	if (!(fout = fopen(fname,"w")))
	  printf("Unable to open file >%s<\n",fname);
	else
	  {
	    s1 = (int)(minx/data->phead.tsamp);
	    s2 = (int)(maxx/data->phead.tsamp);
	    if (s1 < 0) s1=0;
	    if (s2 >= data->phead.nsblk*data->phead.nsub)
	      s2 = data->phead.nsblk*data->phead.nsub-1;
	    
	    writeData(fp,fout,data,s1,s2,dm);	    
	    fclose(fout);
	  }
      }
    else if (key=='h') // Plot histogram of the values
      {
	float *vals;
	int nbin;
	float hist_minx;
	float hist_maxx;
	int count;
	double sx,sx2,mean,sdev;
	float fx[2],fy[2];

	if (data->phead.nbits==1) {hist_minx = -1; hist_maxx = 1; nbin = 4;}
	if (data->phead.nbits==2) {hist_minx = -4; hist_maxx = 4; nbin = 10;}
	if (data->phead.nbits==4) {hist_minx = 0; hist_maxx = 16; nbin = 16;}
	if (data->phead.nbits==8) {hist_minx = 0; hist_maxx = 255; nbin = 255;}
	if (!(vals = (float *)malloc(sizeof(float)*((maxx-minx)/data->phead.tsamp)*(data->phead.nchan))))
	  {
	    printf("Sorry: unable to allocate enough memory - please choose a smaller region\n");
	  }
	else
	  {
	    count = extractPolData(fp,data,1,vals,minx,maxx);
	    cpgend();
	    cpgbeg(0,"4/xs",1,1);
	    sx=0;
	    sx2=0;
	    for (i=0;i<count;i++)
	      {
		sx+=vals[i];
		sx2+=vals[i]*vals[i];
	      }
	    mean = sx/(double)count;
	    sdev = sqrt(sx2/(double)(count)-pow(sx/(double)(count),2));
	    printf("Mean = %g\n",sx/(double)count);
	    printf("Sdev = %g\n",sqrt(sx2/(double)(count)-pow(sx/(double)(count),2)));
	    //	  printf("have: %g\n",vals[i]);
	    draw_histogram(vals,count,nbin,hist_minx,hist_maxx,1,-1,-1,0,0.5,0);
	    cpglab("Level","","");
	    cpgsci(2);
	    fx[0]=fx[1] = mean;
	    fy[0] =0; fy[1] = 1;
	    cpgline(2,fx,fy);
	    fx[0]=fx[1] = mean+sdev;
	    cpgsls(4);
	    cpgline(2,fx,fy);
	    fx[0]=fx[1] = mean-sdev;
	    cpgline(2,fx,fy);
	    cpgsci(1); cpgsls(1);
	    cpgend();
	    cpgbeg(0,grDev,1,1);
	    free(vals);
	  }
      }
    else if (key=='u')
      {
	findXrange=1;
	findYrange=1;
      }
    else if (key=='U')
      {
	findXrange=1;
	findYrange=1;
	np = onp;
	samp1 = osamp1;
	samp2 = osamp2;
	nsmooth = onsmooth;
	for (i=0;i<np;i++)
	  {
	    mean[i] = omean[i];
	    min[i] = omin[i];
	    max[i] = omax[i];
	  }
      }
    else if (key=='f')
      {
	printf("plotting frequency\n");
	cpgend();
	printf("Here minx = %g, maxx = %g\n",minx,maxx);
	plotFrequency(fp,data,minx,maxx);
	cpgbeg(0,grDev,1,1);
	cpgask(0);
      }
    else if (key=='r')
      {
	findXrange=1;
	findYrange=1;
	samp1 = (int)(minx/data->phead.tsamp);
	samp2 = (int)(maxx/data->phead.tsamp);
	printf("New sample range (a): %d %d\n",samp1,samp2);
	if (samp1 < 0) samp1=0;
	if (samp2 >= data->phead.nsblk*data->phead.nsub)
	  samp2 = data->phead.nsblk*data->phead.nsub-1;
	printf("New sample range: %d %d\n",samp1,samp2);
	if (dm==0)
	  np = extractDataZeroDM(fp,data,samp1,samp2,mean,min,max,MAX_RET,&nsmooth);
	else
	  np = extractData(fp,data,samp1,samp2,mean,min,max,MAX_RET,&nsmooth,dm);
	printf("Complete extraction; np = %d, nsmooth = %d\n",np,nsmooth);
      }
    else if (key=='s')
      plotSubint*=-1;
  } while (key!='q');
  cpgend();
  freeDset(data);
  closeFitsFile(fp);
}

void plotFrequency(fitsfile *fp,dSet *d,float mint,float maxt)
{
  float tr[6];
  float min,max;
  int i,j,n=0,i0,i1,p;
  int npts = (int)(((maxt-mint)/d->phead.tsamp)+0.5);
  float *arr;
  float my,mx,mx2,my2;
  char key;
  float minx,maxx,miny,maxy;
  float ominx,omaxx,ominy,omaxy;

  printf("Allocated arr size: %d, npts = %d, %g %g\n",npts*d->phead.nchan,npts,maxt,mint);
  if (!(arr = (float *)malloc(sizeof(float)*npts*d->phead.nchan)))
    {
      printf("plotFreq: unable to allocate memory\n");
      exit(1);
    }
  if (d->phead.npol==1)
    cpgbeg(0,"2/xs",1,1);
  else  
    cpgbeg(0,"2/xs",2,2);
  printf("mint/maxt = %g %g\n",mint,maxt);
  npts = extractPolData(fp,d,0,arr,mint,maxt);

  //  for (i=0;i<npts;i++)
  //    printf("Res: %f\n",arr[i]);

  tr[0] = mint;
  tr[1] = d->phead.tsamp;
  tr[2] = 0;
  tr[3] = d->phead.freq-d->phead.bw/2;
  tr[4] = 0;
  tr[5] =  (float)d->phead.bw/d->phead.nchan;


  /*  tr[0] = 0;
  tr[1] = 1;
  tr[2] = 0;
  tr[3] = 0;
  tr[4] = 0;
  tr[5] = 1;*/


  //  for (i=0;i<1024;i++)
  //    printf("arr = %g\n",arr[i]);

  //  for (p=0;p<d->phead.npol;p++)
    {
      n=0;
      for (i=0;i<npts*d->phead.nchan;i++)
	{
	  if (i==0)
	    {
	      min = max = arr[n];
	    }
	  else 
	    {
	      if (min > arr[n]) min = arr[n];
	      if (max < arr[n]) max = arr[n];
	    }
	  n++;
	}
      printf("min = %g, max = %g, n = %d\n",min,max,n);
      printf("nchan = %d, npts = %d\n",d->phead.nchan,npts);
      //cpgenv(mint,maxt,d->phead.freq-d->phead.bw/2,d->phead.freq+d->phead.bw/2,0,1);
      //      cpgenv(d->phead.freq-d->phead.bw/2,d->phead.freq+d->phead.bw/2,mint,maxt,0,1);
      //      cpgenv(0,3000,0,3000,0,1);
      //cpgenv(0,d->phead.nchan,0,npts,0,1);
      printf("doing the image\n");
      //      cpggray(arr,npts,d->phead.nchan,1,npts,1,d->phead.nchan,min,max,tr);
      cpgask(0);
      ominx = minx = mint;
      omaxx = maxx = maxt;
      ominy = miny = d->phead.freq-d->phead.bw/2;
      omaxy = maxy = d->phead.freq+d->phead.bw/2;
      do {
	cpgenv(minx,maxx,miny,maxy,0,1);
	cpglab("Time (s)","Frequency (MHz)","");
	cpggray(arr,npts,d->phead.nchan,1,npts,1,d->phead.nchan,min,max,tr);
	
	cpgcurs(&mx,&my,&key);

	if (key=='z')
	  {
	    cpgband(2,0,mx,my,&mx2,&my2,&key);
	    if (mx > mx2) {maxx = mx; minx = mx2;}
	    else {maxx = mx2; minx = mx;}
	    
	    if (my > my2) {maxy = my; miny = my2;}
	    else {maxy = my2; miny = my;}	   
	  }
	else if (key=='u')
	  {
	    minx = ominx;
	    maxx = omaxx;
	    miny = ominy;
	    maxy = omaxy;
	  }
	else if (key=='b')
	  {
	    cpgend();
	    cpgbeg(0,"3/xs",1,1);
	    plotBpass(arr,d->phead.nchan,npts);
	    cpgend();
	    cpgbeg(0,"2/xs",1,1);
	    cpgask(0);
	  }
	
      } while (key!='q');
      printf("Done the image\n");
	       }
  cpgend();
  free(arr);
}

void plotBpass(float *arr,int nchan,int npts)
{
  float bpass[nchan],x[nchan];
  int i,j;
  double sum;
  float minx=0;
  float maxx = nchan;
  float maxy=0;
  float miny=0;
  float my,mx,mx2,my2;
  char key;

  float ominx,omaxx,ominy,omaxy;
  cpgask(0);
  for (i=0;i<nchan;i++)
    {
      x[i] = i;
      sum=0.0;
      for (j=0;j<npts;j++)
	sum+=arr[i*npts+j];
      bpass[i] = sum/(double)npts;
      printf("npts for channel %d = %d\n",i,npts);
      if (i==0){ miny = maxy = bpass[i];}
      if (bpass[i] > maxy) maxy = bpass[i];
      if (bpass[i] < miny) miny = bpass[i];
    }
  ominx = minx;
  omaxx = maxx;
  ominy = miny;
  omaxy = maxy;
  do {
    cpgenv(minx,maxx,miny,maxy,0,1);
    cpglab("Channel number","Mean value","");
    cpgbin(nchan,x,bpass,1);
	cpgcurs(&mx,&my,&key);

	if (key=='z')
	  {
	    cpgband(2,0,mx,my,&mx2,&my2,&key);
	    if (mx > mx2) {maxx = mx; minx = mx2;}
	    else {maxx = mx2; minx = mx;}
	    
	    if (my > my2) {maxy = my; miny = my2;}
	    else {maxy = my2; miny = my;}	   
	  }
	else if (key=='u')
	  {
	    minx = ominx;
	    maxx = omaxx;
	    miny = ominy;
	    maxy = omaxy;
	  }
  }while (key!='q');
}

void draw_histogram(float *x,int count,int nbin,float minx,float maxx,int normalise,float maxy,int colour,int log,float offset,int histOutput)
{
  int i,j;
  float binvalX[nbin],binvalY[nbin];
  float binsize;
  float area;
  float highest;
  float scale;
  char dummy[10];
  printf("Offset = %f\n",offset);
  for (i=0;i<nbin;i++)
    {
      binvalX[i]=i*(maxx-minx)/nbin+minx+(maxx-minx)/nbin*offset;
      binvalY[i]=0.0;
    }

  for (i=0;i<count;i++)
    {
      /*      if ((int)((x[i]-minx)/((maxx-minx)/nbin))>=0 
        && (int)((x[i]-minx)/((maxx-minx)/nbin)) <100)
	binvalY[(int)((x[i]-minx)/((maxx-minx)/nbin))]++; */
      for (j=0;j<nbin-1;j++)
	{
	  if (binvalX[j] <= x[i] && binvalX[j+1]>x[i])
	    {
	      binvalY[j]++;
	    }
	}
    }
  if (normalise>0)
    {
      area=0.0;
      for (i=0;i<nbin;i++)
      area+=binvalY[i];
      for (i=0;i<nbin;i++)
      binvalY[i]=binvalY[i]/area*normalise;
    }
  /* Normalise so that highest point lies at 1 */
  highest=0.0;
  for (i=0;i<nbin;i++)
    {
      if (highest<binvalY[i])
	highest=binvalY[i];
    }
  if (maxy<0) maxy = (float)highest+0.1*highest;
  printf("Maxy = %f %d\n",maxy,normalise);
  if (normalise<0)
    {
      if (highest!=0.0)
	{
	  for (i=0;i<nbin;i++)
	    binvalY[i]/=highest;
	}
    }
  if (colour==-1)
    {
      printf("In here with %f %f\n",maxx,maxy);
      if (log==1)
      cpgenv(minx, maxx, 0, maxy, 0, 10); 
      else
      cpgenv(minx, maxx, 0, maxy, 0, 0); 
      cpgsls(1);
      cpgsci(1);
    }
  else if (colour==-2)
    {
      cpgswin(minx,maxx,0,maxy);
      cpgsci(1);
      cpgsls(1);
      cpgslw(1);

    }
  else
    {
      cpgsci(colour);
      /*      cpgsci(1); */
      cpgsls(1);
      cpgslw(1);
      if (colour==2) cpgslw(5);
      if (colour==3) cpgsls(3);
      if (colour==4) cpgsls(5);

    }
  if (histOutput==1)
    {
      for (i=0;i<nbin;i++)
	printf("%d %f %f\n",i,binvalX[i],binvalY[i]);
    }
  cpgbin(nbin,binvalX,binvalY,0);
  cpgsls(1);
  cpgsci(1);
  cpgslw(1);
}

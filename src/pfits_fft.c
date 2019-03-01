/* software to obtain an interactive Fourier transform of a data set */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "pfits.h"
#include "fftw3.h"
#include <cpgplot.h>

void calculateFFT(float *val,int np,double *outSpec);


int main(int argc,char *argv[])
{
  int i;
  char fname[128];
  dSet *data; 
  fitsfile *fp;
  float ts=-1,te=-1;
  float dm=0;
  float *mean,*max,*min;
  int nsmooth,np;
  int samp1=-1,samp2=-1;
  float *outSpecR,*outSpecI,*outSpecP,*x;
  double *outSpec;
  int ptype=1;
  float minx,maxx,miny,maxy;
  float ominx,omaxx,ominy,omaxy;
  float mx,my,mx2,my2;
  char key='x';
  float overlayF;
  int   nHarm=0;
  double meanv=0;
  int nm=0;
  float fx[2],fy[2];
  int scaleY=0;

  data = initialiseDset();

  for (i=0;i<argc;i++)
    {
      if (strcmp(argv[i],"-f")==0)
	strcpy(fname,argv[++i]);
      else if (strcmp(argv[i],"-t1")==0)
	sscanf(argv[++i],"%f",&ts);
      else if (strcmp(argv[i],"-t2")==0)
	sscanf(argv[++i],"%f",&te);
      else if (strcmp(argv[i],"-dm")==0)
	sscanf(argv[++i],"%f",&dm);
      else if (strcmp(argv[i],"-overlay")==0)
	{
	  sscanf(argv[++i],"%f",&overlayF);
	  sscanf(argv[++i],"%d",&nHarm);
	}
	       
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
  printf("Number of points = %d\n",samp2-samp1);

  if (!(min = (float *)malloc(sizeof(float)*(samp2-samp1))))
    {
      printf("Unable to allocate memory\n"); exit(1);
    }
  if (!(max = (float *)malloc(sizeof(float)*(samp2-samp1))))
    {
      printf("Unable to allocate memory\n"); exit(1);
    }
  if (!(mean = (float *)malloc(sizeof(float)*(samp2-samp1))))
      {
      printf("Unable to allocate memory\n"); exit(1);
    }

  printf("samp1 = %d, samp2 = %d\n",samp1,samp2);
  np = extractData(fp,data,samp1,samp2,mean,min,max,samp2-samp1,&nsmooth,dm);
  printf("nsmooth = %d\n",nsmooth);
  printf("np = %d\n",np);
  
  // Calculate the FFT
  printf("Calculating FFT\n");
  outSpecR = (float *)malloc(sizeof(float)*np);
  outSpecI = (float *)malloc(sizeof(float)*np);
  outSpecP = (float *)malloc(sizeof(float)*np);
  x = (float *)malloc(sizeof(float)*np);
  outSpec = (double *)malloc(sizeof(double)*np*2);
  printf("Allocated memory\n");
  calculateFFT(mean,np,outSpec);
  printf("Complete calc\n");
  for (i=0;i<np/2;i++)
    {
      outSpecR[i] = outSpec[2*i];
      outSpecI[i] = outSpec[2*i+1];
      outSpecP[i] = pow(outSpec[2*i],2)+pow(outSpec[2*i+1],2);
      x[i] = (i)/((samp2-samp1)*data->phead.tsamp);
      if (i==0)
	{
	  minx = maxx = x[i];
	  miny = maxy = outSpecP[i];
	}
      else
	{
	  if (i==1){miny = maxy = outSpecP[i];}
	  if (minx > x[i]) minx = x[i];
	  if (maxx < x[i]) maxx = x[i];
	  if (miny > outSpecP[i]) miny = outSpecP[i];
	  if (maxy < outSpecP[i]) maxy = outSpecP[i];
	}
      if (x[i]>110 && x[i] < 990)
	{
	  //	  printf("spec %g\n",outSpecP[i]);
	  meanv+=outSpecP[i];
	  nm++;
	}
    }
  printf("Number of points in mean count = %d\n",nm);
  meanv/=(double)nm;
  ominx = minx;
  omaxx = maxx;
  ominy = miny;
  omaxy = maxy;

  minx = 90;
  maxx = 1010;

  cpgbeg(0,"/xs",1,1);
  cpgask(0);
  do {
    if (scaleY==1)
      {
	int t=1;
	scaleY=0;
	for (i=0;i<np;i++)
	  {
	    if (x[i] > minx && x[i] < maxx)
	      {
		if (t==1)
		  {
		    t=0;
		    if (ptype==1)
		      miny = maxy = outSpecP[i];		    
		    else if (ptype==2)
		      miny = maxy = outSpecR[i];		    
		    else if (ptype==3)
		      miny = maxy = outSpecI[i];		    
		  }
		else
		  {
		    if (ptype==1)
		      {
			if (miny > outSpecP[i]) miny = outSpecP[i];
			if (maxy < outSpecP[i]) maxy = outSpecP[i];
		      }
		    else if (ptype==2)
		      {
			if (miny > outSpecR[i]) miny = outSpecR[i];
			if (maxy < outSpecR[i]) maxy = outSpecR[i];
		      }
		    else if (ptype==3)
		      {
			if (miny > outSpecI[i]) miny = outSpecI[i];
			if (maxy < outSpecI[i]) maxy = outSpecI[i];
		      }
		  }
	      }
	  }
      }
    cpgenv(minx,maxx,miny,maxy,0,1);
    cpglab("Frequency (Hz)","Power","");

    if (ptype==1)
      {
	cpgsci(7); cpgline(np/2,x,outSpecP); cpgsci(1);
      }
    else if (ptype==2)
      {
	cpgsci(7); cpgline(np/2,x,outSpecR); cpgsci(1);
      }
    else if (ptype==3)
      {
	cpgsci(7); cpgline(np/2,x,outSpecI); cpgsci(1);
      }

    fx[0] = minx;
    fx[1] = maxx;
    fy[0] = fy[1] = meanv;
    cpgsci(3); cpgline(2,fx,fy); cpgsci(1);
    fy[0] = fy[1] = meanv*2.996;
    cpgsci(3); cpgline(2,fx,fy); cpgsci(1);
    fy[0] = fy[1] = meanv*4.608;
    cpgsci(3); cpgline(2,fx,fy); cpgsci(1);
    fy[0] = fy[1] = meanv*6.908;
    cpgsci(3); cpgline(2,fx,fy); cpgsci(1);
    fy[0] = fy[1] = meanv*10.29; // Should be 1 point above here
    cpgsci(3); cpgline(2,fx,fy); cpgsci(1);
    if (nHarm>0)
      {
	float fx[2],fy[2];
	fy[0] = miny; fy[1] = maxy;
	for (i=0;i<nHarm;i++)
	  {
	    fx[0] = fx[1] = (i+1)*overlayF;
	    cpgsci(2); cpgsls(4); cpgline(2,fx,fy); cpgsls(1); cpgsci(1);
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
      }
    else if (key=='A')
      {
	printf("frequency = %g, period = %g, val = %g\n",mx,1.0/mx,my);
      }
    else if (key=='u')
      {
	minx = ominx;
	maxx = omaxx;
	miny = ominy;
	maxy = omaxy;
      }
    else if (key=='1')
      {ptype=1; scaleY=1;}
    else if (key=='2')
      {ptype=2; scaleY=1;}
    else if (key=='3')
      {ptype=3; scaleY=1;}
  } while (key!='q');
  cpgend();
  free(x);
  free(outSpecR);
  free(outSpecI);
  free(outSpec);
  free(outSpecP);
  free(min);
  free(max);
  free(mean);
}

void calculateFFT(float *val,int np,double *outSpec)
{
  double *pf;
  fftw_complex *output;
  fftw_plan transform_plan;
  int i;
  printf("Calculating FFT\n");
  pf = (double *)malloc(sizeof(double)*np);
  for (i=0;i<np;i++)
    pf[i] = (double)val[i];

  output = (fftw_complex*)outSpec;
  transform_plan = fftw_plan_dft_r2c_1d(np,pf,output,FFTW_ESTIMATE);
  fftw_execute(transform_plan);
  fftw_destroy_plan(transform_plan);

  free(pf);
}

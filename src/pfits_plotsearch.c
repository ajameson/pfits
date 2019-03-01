//
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "pfits.h"
#include <cpgplot.h>

#define MAX_RET 1024
//#define PFITS_DEBUG

void draw_histogram(float *x,int count,int nbin,float minx,float maxx,int normalise,float maxy,int colour,int log,float offset,int histOutput,int verbose, int plain);
void plotFrequency(float* arr,dSet *d,float mint,float maxt,int npts, int verbose, int* xres, int* yres);
void get_scale (int from, int to, float * width, float * height);
void set_resolution (int width_pixels, int height_pixels);

void help()
{
  printf("-dm <dm>      set the dispersion measure\n");
  printf("-f filename   PSRFITS file name\n");
  printf("-g gr         set graphics device\n");
  printf("-t1 <start>   start time in seconds from beginning of observation\n");
  printf("-t2 <end>     end time in seconds from beginning of observation\n");
  printf("-v            verbose output\n");
  printf("\n\n\n");
}

int main(int argc,char *argv[])
{
  char verbose=0;
  int i;
  char fname[128] = "";
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
  int npts;
  int np;
  long samp1 = 0;
  long samp2 = 10;
  float minx,maxx,miny=-1,maxy=260;
  int nsmooth,onsmooth;
  int findXrange=1;
  int findYrange=1;
  int plotSubint=-1;
  float ts=-10,te=0;
  char key;
  char grDev[128]="/png";
  int interactive=1;
  int stopit=0;
  fitsfile *fp;
  float tx[2],ty[2];

#ifdef DEBUG
  help();
#endif

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
    else if (strcmp(argv[i],"-v")==0)
      verbose=1;
  }

  if (strlen(fname) == 0)
  {
    fprintf (stderr, "ERROR: no filename provided\n");
    help();
    return EXIT_FAILURE;
  }

  fp   = openFitsFile(fname); 
  loadPrimaryHeader(fp,data);
#ifdef PFITS_DEBUG
  displayHeaderInfo(data);
#endif

  samp1 = ts/data->phead.tsamp;
  samp2 = te/data->phead.tsamp;

  if (samp1 < 0){
    samp1+=data->phead.nsub*data->phead.nsblk;
  }
  if (samp2 < 0){
    samp2+=data->phead.nsub*data->phead.nsblk;
  }

  if (samp2 == 0){
    te=samp2=data->phead.nsub*data->phead.nsblk;
  }

  if (dm==0)
    np = extractDataZeroDM(fp,data,samp1,samp2,mean,min,max,MAX_RET,&nsmooth);
  else
    np = extractData(fp,data,samp1,samp2,mean,min,max,MAX_RET,&nsmooth,dm);

  if (stopit==1)
    exit(1);

  int xres[2] = {300, 1024};
  int yres[2] = {225, 768};

#ifdef PFITS_DEBUG
  printf("np = %d\n",np);
#endif
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
  char devname[128];
  unsigned ires;
  for (ires=0; ires<2; ires++)
  {
    sprintf(devname, "pfits_timeseries.%dx%d.png/png", xres[ires], yres[ires]);
    if (cpgopen (devname) != 1)
      fprintf (stderr, "error opening plot device\n");
    set_resolution (xres[ires], yres[ires]);

    if (xres[ires] >= 640)
    {
      cpgenv(minx,maxx,miny,maxy,0,1);
      cpglab("Time (s)","","");
    }
    else
    {
      cpgswin(minx,maxx,miny,maxy);
      cpgsvp(0, 1, 0, 1);
      cpgbox("BCNST", 0.0, 0.0, "BCNST", 0.0, 0.0);
    }
    if (verbose)
      printf("nsmooth = %d\n",nsmooth);
    if (nsmooth > 1)
    {
      cpgsci(5); cpgbin(np,fx,min,1);
      cpgsci(7); cpgbin(np,fx,max,1);
    }
    cpgsci(1); cpgbin(np,fx,mean,1);
    if (plotSubint==1)
    {
      for (i=0;i<data->phead.nsub;i++)
      {
        tx[0] = tx[1] = i*data->phead.nsblk*data->phead.tsamp;
        ty[0] = miny; ty[1] = maxy;
        cpgsci(7); cpgsls(4); cpgline(2,tx,ty); cpgsci(1); cpgsls(1);
      }
    }
    cpgclos();
  }

  float *vals;
  if (verbose)
    printf("ALLOCATING %d floats\n",(int)((maxx-minx)/data->phead.tsamp)*(data->phead.nchan));
  if (!(vals = (float *)malloc(sizeof(float)*(int)(((maxx-minx+0.5)/data->phead.tsamp)*(data->phead.nchan)))))
  {
    printf("Sorry: unable to allocate enough memory - please choose a smaller region\n");
    exit(2);
  }
  npts = extractPolData(fp,data,0,vals,minx,maxx);

  char plot_histogram = 1;
  if (plot_histogram)
  {
    int nbin;
    float hist_minx;
    float hist_maxx;
    int count=npts;
    double sx,sx2,mean,sdev;
    float fx[2],fy[2];

    if (data->phead.nbits==1) {hist_minx = -1;hist_maxx = 1; nbin = 4;}
    if (data->phead.nbits==2) {hist_minx = -4; hist_maxx = 4; nbin = 10;}
    if (data->phead.nbits==4) {hist_minx = 0; hist_maxx = 16; nbin = 16;}
    if (data->phead.nbits==8) {hist_minx = 0; hist_maxx = 255; nbin = 255;}
    count *= data->phead.nchan;
    sx=0;
    sx2=0;
    for (i=0;i<count;i++)
    {
      sx+=vals[i];
      sx2+=vals[i]*vals[i];
    }
    mean = sx/(double)count;
    sdev = sqrt(sx2/(double)(count)-pow(sx/(double)(count),2));
    if (verbose)
    {
      printf("Mean = %g\n",sx/(double)count);
      printf("Sdev = %g\n",sqrt(sx2/(double)(count)-pow(sx/(double)(count),2)));
    }
    for (ires=0; ires<2; ires++)
    {
      sprintf(devname, "pfits_histogram.%dx%d.png/png", xres[ires], yres[ires]);
      if (cpgopen (devname) != 1)
        fprintf (stderr, "error opening plot device\n");
      set_resolution (xres[ires], yres[ires]);

      int plain = xres[ires] < 640;
      if (verbose)
      {
        fprintf (stdout, "draw_histogram count=%d nbin=%d minx=%f maxx=%f normalize=%d maxy=%f color=%d log=%d offset=%f histOutput=%d verbose=%d plain=%d\n", count, nbin, hist_minx, hist_maxx, 1, -1, -3, 0, 0.5, 0, verbose, plain);
      }
      draw_histogram(vals,count,nbin,hist_minx,hist_maxx,1,-1,-3,0,0.5,0,verbose, plain);
      cpgslw(4);
      if (xres[ires] >= 640)
        cpglab("Level","","");
      cpgsci(7);
      fx[0]=fx[1] = mean;
      fy[0] =0; fy[1] = 1;
      cpgline(2,fx,fy);
      fx[0]=fx[1] = mean+sdev;
      cpgsls(4);
      cpgline(2,fx,fy);
      fx[0]=fx[1] = mean-sdev;
      cpgline(2,fx,fy);
      cpgsci(1); cpgsls(1);
      cpgclos();
    }
  }

  {
    if (verbose)
    {
      printf("plotting frequency\n");
      printf("Here minx = %g, maxx = %g\n",minx,maxx);
    }
    plotFrequency(vals,data,minx,maxx,npts,verbose,xres,yres);
  }

  free(vals);
  freeDset(data);
  closeFitsFile(fp);
}

void plotFrequency(float* arr,dSet *d,float mint,float maxt,int npts,int verbose, int* xres, int* yres)
{
  float tr[6];
  float min,max;
  int i,j,n=0,i0,i1,p;
  float my,mx,mx2,my2;
  char key;
  float minx,maxx,miny,maxy;
  float ominx,omaxx,ominy,omaxy;

  if (verbose)
    printf("mint/maxt = %g %g\n",mint,maxt);

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
  
  char devname[128];
  int ires;
  for (ires=0; ires<2; ires++)
  {
    sprintf(devname, "pfits_freqtime.%dx%d.png/png", xres[ires], yres[ires]);
    if (d->phead.npol==1)
      cpgbeg(0,devname,1,1);
    else
      cpgbeg(0,devname,2,2);
  
    set_resolution (xres[ires], yres[ires]);

    for (p=0;p<d->phead.npol;p++)
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
      if (verbose)
      {
        printf("min = %g, max = %g, n = %d\n",min,max,n);
        printf("nchan = %d, npts = %d\n",d->phead.nchan,npts);
      }
      //cpgenv(mint,maxt,d->phead.freq-d->phead.bw/2,d->phead.freq+d->phead.bw/2,0,1);
      //      cpgenv(d->phead.freq-d->phead.bw/2,d->phead.freq+d->phead.bw/2,mint,maxt,0,1);
      //      cpgenv(0,3000,0,3000,0,1);
      //cpgenv(0,d->phead.nchan,0,npts,0,1);
      if (verbose)
        printf("doing the image\n");
      //      cpggray(arr,npts,d->phead.nchan,1,npts,1,d->phead.nchan,min,max,tr);
      cpgask(0);
      ominx = minx = mint;
      omaxx = maxx = maxt;
      ominy = miny = d->phead.freq-d->phead.bw/2;
      omaxy = maxy = d->phead.freq+d->phead.bw/2;

      if (xres[ires] >= 640)
      {
        cpgenv(minx,maxx,miny,maxy,0,1);
        cpglab("Time (s)","Frequency (MHz)","");
      }
      else
      {
        cpgpage();
        cpgsvp(0.0,1.0,0.0,1.0);
        cpgswin(minx,maxx,miny,maxy);
      }
      cpggray(arr,npts,d->phead.nchan,1,npts,1,d->phead.nchan,min,max,tr);
      if (verbose)
        printf("Done the image\n");
    }
  }
  cpgend();
}

void draw_histogram (float *x,int count,int nbin,float minx,float maxx,int normalise,float maxy,int colour,int log,float offset,int histOutput,int verbose, int plain)
{
  if (verbose)
    printf ("draw_histogram minx=%f maxx=%f\n", minx, maxx);

  int i,j;
  float binvalX[nbin+1],binvalY[nbin+1];
  unsigned long binvalY_int[nbin+1];
  float binsize;
  float area;
  float highest;
  float scale;
  char dummy[10];
  if (verbose)
  {
    printf("Offset = %f\n",offset);
    printf("nbin=%d\n", nbin);
  }
  // configure the X and Y values of the histogram
  for (i=0;i<nbin;i++)
  {
    binvalX[i+1]=i*(maxx-minx)/nbin+minx+(maxx-minx)/nbin*offset;
    binvalY[i+1]=0.0;
    binvalY_int[i+1]=0.0;
  }
  binvalY_int[0]=0;
  binvalY[0]=0;
  binvalX[0]=2*binvalX[1]-binvalX[2];
  nbin++; // we add this extra bin on the front to make pgplot behave

  const float bin_width = (maxx - minx) / (nbin -1);
  if (verbose)
    printf ("range=%f bin_width=%f nbin=%d\n", (maxx - minx), bin_width, nbin);

  // process each value adding it to the appropriate bins
  if (verbose)
    printf("count=%d\n", count);
  for (i=0;i<count;i++)
  {
    /*
    float o = x[i] - minx;
    float p = o / bin_width;
    float q = p + 0.0;
    float r = rintf(q);
    int   j = (int) r;
    if (i < 100)
      printf ("x=%f o=%f p=%f q=%f r=%f s=%d\n", x[i], o, p, q, r, j);
      */

    // determine the bin for x[i]
    j = (int) rintf(((x[i] - minx) / bin_width));

    if (j < 0)
      j = 0;
    if (j >= nbin)
      j = nbin-1;
    binvalY_int[j]++;
  }

  if (normalise>0)
  {
    area=0.0;
    for (i=0;i<nbin;i++)
      area += (float) binvalY_int[i];
    for (i=0;i<nbin;i++)
    {
      binvalY[i]=((float) binvalY_int[i])/area*normalise;
    }
  }

  /* Normalise so that highest point lies at 1 */
  highest=0.0;
  for (i=0;i<nbin;i++)
  {
    if (highest<binvalY[i])
      highest=binvalY[i];
  }
  if (maxy<0) maxy = (float)highest+0.1*highest;
  if (verbose)
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
    if (verbose)
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
  else if (colour==-3)
  {
    if (verbose)
      printf("In here with %f %f\n",maxx,maxy);
    if (plain)
    {
      cpgsvp(0.0,1.0,0.0,1.0);
      cpgswin(minx,maxx,0,maxy);
      if (log == 1)
        cpgbox("BCNST", 0.0, 0.0, "BCNSTL", 0.0, 0.0);
      else
        cpgbox("BCNST", 0.0, 0.0, "BCNST", 0.0, 0.0);
    }
    else
    {
      if (log==1)
        cpgenv(minx, maxx, 0, maxy, 0, 10); 
      else
        cpgenv(minx, maxx, 0, maxy, 0, 0); 
    }
    cpgsls(1);
    cpgsci(1);
    cpgslw(8);
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
  if (verbose && histOutput==1)
  {
    for (i=0;i<nbin;i++)
      printf("%d %f %f\n",i,binvalX[i],binvalY[i]);
  }
  cpgbin(nbin,binvalX,binvalY,0);
  cpgsls(1);
  cpgsci(1);
  cpgslw(1);
}

void get_scale (int from, int to, float * width, float * height)
{
  float j = 0;
  float fx, fy;
  cpgqvsz (from, &j, &fx, &j, &fy);

  float tx, ty;
  cpgqvsz (to, &j, &tx, &j, &ty);

  *width = tx / fx;
  *height = ty / fy;
}


void set_resolution (int width_pixels, int height_pixels)
{
  float width_scale, height_scale;
  width_pixels--;
  height_pixels--;

  get_scale (3, 1, &width_scale, &height_scale);

  float width_inches = width_pixels * width_scale;
  float aspect_ratio = height_pixels * height_scale / width_inches;

  cpgpap( width_inches, aspect_ratio );

  float x1, x2, y1, y2;
  cpgqvsz (1, &x1, &x2, &y1, &y2);
}


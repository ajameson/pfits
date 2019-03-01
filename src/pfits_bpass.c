// Routines to plot the bandpass stored in a psrfits file
//
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "pfits.h"
#include <cpgplot.h>

#define MAX_OVERLAY 1000

void help()
{
  printf("-f <filename>   File to process\n");
  printf("-noc1           Ignore the first channel\n");
  printf("-g <grdev>      Set the graphics device\n");
  exit(1);
}

void sortInt(int *zapChannels,int nzap);

int main(int argc,char *argv[])
{
  int i;
  char fname[128];
  dSet *data; 
  float freq,bw,chanbw;
  int nchan,npol;
  float bpass[4096];
  float fx[4096];
  float miny,maxy,minx,maxx;
  float ominy,omaxy,ominx,omaxx;
  float mx,my,mx2,my2;
  float binw;
  char key;
  char grDev[128]="/xs";
  int interactive=1;
  int noc1=0;
  int zapChannels[4096];
  int nzap=0;
  int overlay=-1;
  float overlayVal[MAX_OVERLAY];
  char overlayStr[MAX_OVERLAY][128];
  char overlayFile[128];
  int noverlay=0;
  fitsfile *fp;

  data = initialiseDset();
  

  for (i=0;i<argc;i++)
    {
      if (strcmp(argv[i],"-f")==0)
	strcpy(fname,argv[++i]);
      else if (strcmp(argv[i],"-noc1")==0)
	noc1=1;
      else if (strcmp(argv[i],"-g")==0)
	{
	  strcpy(grDev,argv[++i]);
	  interactive=0;
	}
      else if (strcmp(argv[i],"-h")==0)
	help();
      else if (strcmp(argv[i],"-overlay")==0)
	{
	  strcpy(overlayFile,argv[++i]);
	  overlay=1;
	}
    }
  if (overlay==1)
    {
      FILE *fin;
      char line[1024];
      noverlay=0;

      if (!(fin = fopen(overlayFile,"r")))
	  printf("Unable to open overlay file >%s<\n",overlayFile);
      else
	{
	  while (!feof(fin))
	    {
	      fgets(overlayStr[noverlay],1024,fin);
	      if (fscanf(fin,"%f",&overlayVal[noverlay])==1)
		{
		  if (overlayStr[noverlay][strlen(overlayStr[noverlay])-1] == '\n')
		    overlayStr[noverlay][strlen(overlayStr[noverlay])-1]='\0';
		  noverlay++;
		}
	    }
	  fclose(fin);
	}
    }
  fp   = openFitsFile(fname); 
  loadPrimaryHeader(fp,data);
  displayHeaderInfo(data);
  readBandpass(fp,bpass);
  nchan = data->phead.nchan;
  freq  = data->phead.freq;
  bw    = data->phead.bw;
  chanbw = data->phead.chanbw;

  for (i=0;i<nchan;i++)
    {
      fx[i] = freq-bw/2+(i+0.5)*chanbw;
      if (i==noc1)
	{
	  miny = maxy = bpass[i];
	  minx = maxx = fx[i];
	}
      else if (i!=0)
	{
	  if (bpass[i] > maxy) maxy = bpass[i];
	  if (bpass[i] < miny) miny = bpass[i];
	  if (fx[i] > maxx) maxx = fx[i];
	  if (fx[i] < minx) minx = fx[i];
	}
    }
  ominx = minx;
  omaxx = maxx;
  ominy = miny;
  omaxy = maxy;
  binw = fx[1]-fx[0];
  printf("Complete\n");

  cpgbeg(0,grDev,1,1);
  cpgask(0);
  do {
    cpgenv(minx,maxx,miny,maxy,0,1);
    cpglab("Frequency (MHz)","Amplitude (arbitrary)",fname);
    cpgbin(nchan-noc1,fx+noc1,bpass+noc1,0);
    if (overlay==1)
      {
	float tx[2],ty[2];
	cpgsls(4); cpgsci(2); cpgsch(0.8);
	for (i=0;i<noverlay;i++)
	  {
	    tx[0] = tx[1] = overlayVal[i];
	    ty[0] = miny;
	    ty[1] = maxy;
	    if (tx[1] > minx && tx[1] < maxx)
	      {
		cpgline(2,tx,ty);
		//		cpgtext(tx[1],ty[1]-0.05*(maxy-miny),overlayStr[i]);
		cpgptxt(tx[1]-0.004*(maxx-minx),ty[0]+0.05*(maxy-miny),90,0.0,overlayStr[i]);
	      }
	  }
	cpgsci(1); cpgsls(1); cpgsch(1);
      }
    if (interactive==1)
      {
	cpgcurs(&mx,&my,&key);
	if (key=='A')
	  {
	    int cc=-1;
	    int i;
	    for (i=0;i<nchan-1;i++)
	      {
		//		if ((bw > 0 && (mx > fx[i]-binw/2 && mx < fx[i]+binw/2)) ||
		//		    (bw < 0 && (mx > fx[i]+binw/2 && mx < fx[i]-binw/2)))
		if ((bw > 0 && (mx > fx[i] && mx < fx[i]+binw)) ||
		    (bw < 0 && (mx > fx[i] && mx < fx[i]+binw)))
		  {
		    cc = i;
		    break;
		  }
	      }
	    printf("mouse x = %g MHz, mouse y = %g, channel = %d, channel frequency = %g MHz\n",mx,my,cc,fx[cc]);
	  }
	else if (key=='X')
	  {
	    int cc=-1;
	    int i;
	    printf("Deleting %g %g %g\n",mx,fx[10],binw);
	    for (i=0;i<nchan-1;i++)
	      {
		//		if ((bw > 0 && (mx > fx[i]-binw/2 && mx < fx[i]+binw/2)) ||
		//		    (bw < 0 && (mx > fx[i]+binw/2 && mx < fx[i]-binw/2)))
		if ((bw > 0 && (mx > fx[i] && mx < fx[i]+binw)) ||
		    (bw < 0 && (mx > fx[i] && mx < fx[i]+binw)))
		  {
		    cc = i;
		    break;
		  }
	      }
	    printf("Want to delete = %d\n",cc);
	    if (cc != -1)
	      {
		bpass[cc] = 0;
		omaxy = bpass[noc1];
		zapChannels[nzap++] = cc;
		for (i=noc1;i<nchan;i++)
		  {
		    if (omaxy < bpass[i]) omaxy = bpass[i];
		  }
	      }
	  }
	else if (key=='z')
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
	else if (key=='l') // List the channels and frequencies to zap
	  {
	    int i;
	    sortInt(zapChannels,nzap);
	    printf("-------------------------------------------------------\n");
	    printf("Zap channels with first channel = 0\n\n");
	    for (i=0;i<nzap;i++)
	      printf("%d ",zapChannels[i]);
	    printf("\n\n");
	    printf("Zap channels with first channel = 1\n\n");
	    for (i=0;i<nzap;i++)
	      printf("%d ",zapChannels[i]+1);
	    printf("\n\n");
	    printf("Zap channels frequencies:\n\n");
	    for (i=0;i<nzap;i++)
	      printf("%g ",fx[zapChannels[i]]);
	    printf("\n\n");
	    printf("-------------------------------------------------------\n");
	  }
	else if (key=='%') // Enter percentage of the band edges to zap
	  {
	    float percent;
	    int i;

	    printf("Enter band edge percentage to zap ");
	    scanf("%f",&percent);
	    for (i=0;i<nchan;i++)
	      {
		if (i < nchan*percent/100.0 || i > nchan-(nchan*percent/100.0))
		  {
		    bpass[i] = 0;		    
		    zapChannels[nzap++] = i;

		  }
	      }
	    omaxy = bpass[noc1];
	    for (i=noc1;i<nchan;i++)
	      {
		if (omaxy < bpass[i]) omaxy = bpass[i];
	      }

	    // Unzoom
	    minx = ominx;
	    maxx = omaxx;
	    miny = ominy;
	    maxy = omaxy;
	  }

      }
  } while (key != 'q' && interactive==1);
  cpgend();
}

void sortInt(int *zapChannels,int nzap)
{
  int change,t1,t2;
  int i;

  do {
  change=0;
    for (i=0;i<nzap-1;i++)
      {
	if (zapChannels[i] > zapChannels[i+1])
	  {
	    t1 = zapChannels[i];
	    zapChannels[i] = zapChannels[i+1];
	    zapChannels[i+1] = t1;
	    change=1;
	  }
	
      }
  } while (change == 1);
}

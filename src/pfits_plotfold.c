#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "pfits.h"
#include <cpgplot.h>

int main(int argc,char *argv[])
{
  int i,j;
  char fname[128];
  dSet *data;
  fitsfile *fp;
  float *fx,*fy,*freq_y,*time_y,*bpass,*bx;
  float *fy_R,*freq_y_R,*time_y_R;
  float tx[4096],ty[4096];
  float minx,maxx,miny,maxy;
  float tr[6];
  float dm =-1; // Read from the header
  int png=-1;
  double tdelay;
  int cdelay;
  double chanbw;
  double f0;
  double bintime;
  int writeProfile=0;
  FILE *fout;
  int n=0;
  int wcache=0;
  int lcache=0;
  int sub0=0;
  char cacheFileWrite[128];
  char cacheFileLoad[128];
  double rotate=0.0;

  data = initialiseDset();

  for (i=0;i<argc;i++)
    {
      if (strcmp(argv[i],"-f")==0)
	strcpy(fname,argv[++i]);
      else if (strcmp(argv[i],"-png")==0)
	png=1;
      else if (strcmp(argv[i],"-wprofile")==0)
	writeProfile=1;
      else if (strcmp(argv[i],"-wcache")==0)
	{
	  strcpy(cacheFileWrite,argv[++i]);
	  wcache=1;
	}
      else if (strcmp(argv[i],"-lcache")==0)
	{
	  strcpy(cacheFileLoad,argv[++i]);
	  lcache=1;
	}
      else if (strcmp(argv[i],"-r")==0)
	{
	  rotate=atof(argv[++i]);
	}

    }
  printf("Opening file\n");
  fp   = openFitsFile(fname);
  printf("Loading header\n");
  loadPrimaryHeader(fp,data);
  fx = (float *)malloc(sizeof(float)*data->phead.nbin);
  fy = (float *)malloc(sizeof(float)*data->phead.nbin);
  fy_R = (float *)malloc(sizeof(float)*data->phead.nbin);
  bx = (float *)malloc(sizeof(float)*data->phead.nchan*2);
  bpass = (float *)malloc(sizeof(float)*data->phead.nchan*2);
  freq_y = (float *)malloc(sizeof(float)*data->phead.nbin*data->phead.nchan*2);
  time_y = (float *)malloc(sizeof(float)*data->phead.nbin*data->phead.nsub*2);
  freq_y_R = (float *)malloc(sizeof(float)*data->phead.nbin*data->phead.nchan*2);
  time_y_R = (float *)malloc(sizeof(float)*data->phead.nbin*data->phead.nsub*2);
  if (lcache==1)
    {
      FILE *fin;
      // Check and load the cache
      if (!(fin = fopen(cacheFileLoad,"r")))
	printf("Unable to load the cache file >%s<\n",cacheFileLoad);
      else
	{
	  char origName[128],temp[128];
	  int  origNsub,origNbin,origNchan;

	  /*	  fscanf(fin,"%s",origName);
	  if (strcmp(origName,fname)!=0)
	    {
	      printf("Warning: the file used to produce the cache (%s) has a different name to the current file (%s)\n",origName,fname);
	    }
	  */
	  fread(&origNsub,sizeof(int),1,fin);
	  printf("Number of subintegrations in cache = %d, new number of subintegrations = %d\n",origNsub,data->phead.nsub);
	  sub0 = origNsub;

	  fread(&origNbin,sizeof(int),1,fin);
	  if (origNbin != data->phead.nbin)
	    {
	      printf("ERROR: cached data has a different number of bins to the new data\n");
	      lcache=0;
	    }
	  fread(&origNchan,sizeof(int),1,fin);
	  if (origNchan != data->phead.nchan)
	    {
	      printf("ERROR: cached data has a different number of channels to the new data\n");
	      lcache=0;
	    }
	  if (lcache==1)
	    {
	      int i,j;
	      printf("Still loading the cache %d %d %d\n",data->phead.nchan,data->phead.nbin,data->phead.nsub);	      

	      for (i=0;i<origNchan;i++)
		{
		  for (j=0;j<origNbin;j++)
		    fread(&freq_y[i*data->phead.nbin+j],sizeof(float),1,fin);		
		}
	      printf("Loaded freq-phase\n");
	      for (i=0;i<origNsub;i++)
		{
		  for (j=0;j<origNbin;j++)
		    fread(&time_y[i*data->phead.nbin+j],sizeof(float),1,fin);
		}
	      printf("Loaded time-phase\n");
	      for (i=0;i<origNbin;i++)
		{
		  fread(&fy[i],sizeof(float),1,fin);
		  fx[i] = i;
		}
	      printf("Loaded profile\n");
 	    }

	  fclose(fin);
	}
    }

  //  if (lcache==0)
  extractFoldData(fp,data,dm,fx,fy,freq_y,time_y,bpass,sub0);

  for(i=0; i < data->phead.nchan; i++){
	  bx[i]=i%data->phead.nchan;
  }


  int R=rotate*data->phead.nbin;
  printf("Rotate by %lf phase, %d bins\n",rotate,R);
  memcpy(fy_R,fy+R ,(data->phead.nbin-R)*sizeof(float));
  memcpy(fy_R+(data->phead.nbin-R),fy,R*sizeof(float));
  for(i =0 ; i < data->phead.nchan; i++){
  memcpy(freq_y_R+i*data->phead.nbin,freq_y+i*data->phead.nbin+R ,(data->phead.nbin-R)*sizeof(float));
  memcpy(freq_y_R+i*data->phead.nbin+R,freq_y+i*data->phead.nbin ,R*sizeof(float));
  }
for(i =0 ; i < data->phead.nsub; i++){
  memcpy(time_y_R+i*data->phead.nbin,time_y+i*data->phead.nbin+R ,(data->phead.nbin-R)*sizeof(float));
  memcpy(time_y_R+i*data->phead.nbin+R,time_y+i*data->phead.nbin ,R*sizeof(float));
  }


  
  printf("Displaying header\n");
  displayHeaderInfo(data);
  
  miny = maxy = freq_y[0];
  for (i=0;i<data->phead.nbin*data->phead.nchan;i++)
    {
      if (miny > freq_y[i]) miny = freq_y[i];
      if (maxy < freq_y[i]) maxy = freq_y[i];
    }
  if (png==1)
    cpgbeg(0,"plot1.png/png",1,1);
  else
    cpgbeg(0,"1/xs",1,1);

  cpgenv(0,data->phead.nbin,0,data->phead.nchan,0,1);
  cpglab("Bin number","Frequency channel",fname);
  printf("Freq: min/max = %g %g\n",miny,maxy);

  tr[0] = 0;
  tr[1] = 1;
  tr[2] = 0;
  tr[3] = 0;
  tr[4] = 0;
  tr[5] = 1;



/* 
    float heat_l[] = {0.0, 0.2, 0.4, 0.6, 1.0};
    float heat_r[] = {0.0, 0.5, 1.0, 1.0, 1.0};
    float heat_g[] = {0.0, 0.0, 0.5, 1.0, 1.0};
    float heat_b[] = {0.0, 0.0, 0.0, 0.3, 1.0};
*/

    float heat_l[] = {0.0, 0.2, 0.4, 0.6, 1.0};
    float heat_r[] = {0.0, 0.5, 1.0, 1.0, 1.0};
    float heat_g[] = {0.0, 0.0, 0.5, 1.0, 1.0};
    float heat_b[] = {0.0, 0.0, 0.0, 0.3, 1.0};


    cpgctab (heat_l, heat_r, heat_g, heat_b, 5, 1.0,0.5);


  cpgimag(freq_y_R,data->phead.nbin,data->phead.nchan,1,data->phead.nbin,1,data->phead.nchan,miny,maxy,tr);
  // Overlay the DM line
  if (dm < 0)
    dm = data->phead.dm;
  // Want frequency of channel 0
  f0 = data->phead.freq-data->phead.chanbw*data->phead.nchan/2;
  chanbw = data->phead.chanbw;
  bintime = (double)data->phead.period/(double)data->phead.nbin;
  for (i=0;i<data->phead.nchan;i++)
    {
      tdelay = 4.15e-3*dm*(pow(f0/1000.0,-2)-pow((f0+chanbw*i)/1000.0,-2));
      cdelay = nint(-tdelay/bintime);
      while (cdelay >= data->phead.nbin)
      	cdelay -= data->phead.nbin;
      while (cdelay < 0)
      	cdelay += data->phead.nbin;
      //      printf("At this point: %d %d %g %g\n",i,cdelay,tdelay,bintime);
      if (n > 0 && fabs(cdelay-tx[n-1]) > 3)
	{
	  cpgsci(2); cpgline(n,tx,ty);
	  n=0;
	  tx[n] = (float)cdelay;
	  ty[n] = i;
	}
      else
	{
	  tx[n] = (float)cdelay;
	  ty[n] = i;
	  n++;
	}
      
    }

  cpgsci(2); cpgline(n,tx,ty);
  cpgend();
  

  tr[0] = 0;
  tr[1] = 1;
  tr[2] = 0;
  tr[3] = -0.5;
  tr[4] = 0;
  tr[5] = 1;

  miny = maxy = time_y[0];
  for (i=0;i<data->phead.nbin*data->phead.nsub;i++)
    {
      if (miny > time_y[i]) miny = time_y[i];
      if (maxy < time_y[i]) maxy = time_y[i];
    }
  printf("Time: min/max = %g %g\n",miny,maxy);
  if (png == 1)
    cpgbeg(0,"plot2.png/png",1,1);
  else
    cpgbeg(0,"2/xs",1,1);
  //  cpgbeg(0,"plot2.png/png",1,1);
  cpgenv(0,data->phead.nbin,0,data->phead.nsub,0,1);
  cpglab("Bin number","Subintegration",fname);

    cpgctab (heat_l, heat_r, heat_g, heat_b, 5, 1.0,0.5);

  cpgimag(time_y_R,data->phead.nbin,data->phead.nsub,1,data->phead.nbin,1,data->phead.nsub,miny,maxy,tr);
  cpgend();

  if (writeProfile==1)
    fout = fopen("profile.dat","w");


  // work out approx S/N
  
  float mean;
  float rms;
  getbaseline(fy,data->phead.nbin,0.35,&mean,&rms);

  int window;
  float sum;
  int bestwindow=1;
  int bestp0=0;
  float best=0;
  for (window=1; window < 0.5*data->phead.nbin; window*=2){
	  sum=0;
	  for (i=0;i<window;i++){
		  sum+=fy[i]-mean;
	  }
	  for (i=0;i<data->phead.nbin;i++){
		  if (sum/sqrt(window) > best){
			  bestp0=i;
			  bestwindow=window;
			  best=sum/sqrt(window);
		  }
		  sum-=fy[i]-mean;
		  sum+=fy[(window+i)%data->phead.nbin]-mean;
	  }
  }
  for (window=bestwindow/2; window < bestwindow*2; window++){
	  sum=0;
	  for (i=0;i<window;i++){
		  sum+=fy[i]-mean;
	  }
	  for (i=0;i<data->phead.nbin;i++){
		  if (sum/sqrt(window) > best){
			  bestp0=i;
			  bestwindow=window;
			  best=sum/sqrt(window);
		  }
		  sum-=fy[i]-mean;
		  sum+=fy[(window+i)%data->phead.nbin]-mean;
	  }
  }

  printf("bw: %d %d\n",bestwindow,bestp0);
  printf("S/N: %.1f Nsub: %d\n",best/rms,data->phead.nsub);

  for (i=0;i<data->phead.nbin;i++){
	  fy_R[i]/=rms;
	  fy_R[i]-=mean;
  }




  for (i=0;i<data->phead.nbin;i++)
    {
      if (i==0)
	{
	  minx = maxx = fx[i];
	  miny = maxy = fy_R[i];
	}
      if (minx > fx[i]) minx = fx[i];
      if (maxx < fx[i]) maxx = fx[i];
      if (miny > fy_R[i]) miny = fy_R[i];
      if (maxy < fy_R[i]) maxy = fy_R[i];
      if (writeProfile==1)
	fprintf(fout,"%g %g\n",fx[i],fy[i]);
    }
  if (writeProfile==1)
    fclose(fout);
  if (png==1)
    cpgbeg(0,"plot3.png/png",1,1);
  else
    cpgbeg(0,"3/xs",1,1);
  //
  cpgenv(minx,maxx,miny,maxy,0,0);
  cpglab("Bin number","Signal",fname);
  cpgline(data->phead.nbin,fx,fy_R);
  cpgend();

  for (i=0;i<data->phead.nchan*2;i++)
    {
      if (i==0)
	{
	  minx = maxx = bx[i];
	  miny = maxy = bpass[i];
	}
      if (minx > bx[i]) minx = bx[i];
      if (maxx < bx[i]) maxx = bx[i];
      if (miny > bpass[i]) miny = bpass[i];
      if (maxy < bpass[i]) maxy = bpass[i];
    }

  if (png==1)
    cpgbeg(0,"plot4.png/png",1,1);
  else
    cpgbeg(0,"4/xs",1,1);
  //
  cpgenv(minx,maxx,miny,maxy,0,0);
  cpglab("Frequency Channel","Amplitude",fname);
  cpgline(data->phead.nchan,bx,bpass);
  cpgsls(2);
  cpgline(data->phead.nchan,bx,bpass+data->phead.nchan);
  cpgend();



  // Write a cache file if requested
  if (wcache==1)
    {
      FILE *fout;
      if (!(fout = fopen(cacheFileWrite,"wb")))
	printf("Unable to open cache file: %s\n",cacheFileWrite);
      else
	{
	  // Write header information
	  fwrite(&(data->phead.nsub),sizeof(int),1,fout);
	  fwrite(&(data->phead.nbin),sizeof(int),1,fout);
	  fwrite(&(data->phead.nchan),sizeof(int),1,fout);

	  // Write the freq-phase plot
	  for (i=0;i<data->phead.nchan;i++)
	    {
	      for (j=0;j<data->phead.nbin;j++)
		{
		  fwrite(&freq_y[i*data->phead.nbin+j],sizeof(float),1,fout);
		  if (i==0 && j < 10) printf("Have: %g\n",freq_y[i*data->phead.nbin+j]);
		}
	    }
	  for (i=0;i<data->phead.nsub;i++)
	    {
	      for (j=0;j<data->phead.nbin;j++)
		fwrite(&time_y[i*data->phead.nbin+j],sizeof(float),1,fout);
	    }
	  for (i=0;i<data->phead.nbin;i++)
	    fwrite(&fy[i],sizeof(float),1,fout);
	  fclose(fout);
	}     
    }
  
  freeDset(data);
  closeFitsFile(fp);
  free(fx); free(fy); free(freq_y); free(time_y);
}

// Code to read in one or more FITS files and process them in some way to produce
// one or more output files
//
// gcc -lm -o fitsProcess fitsProcess.c
//
// Examples:
// To extract a specific frequency band from a PSRFITS fold-mode file:
// ./fitsProcess ofileID=file1 selectfreq=650,710 ifile=myfile.rf
//
// Step 1 combining or removing frequency channels
//
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include "fitsio.h"

#define MAX_FILENAME_LEN 1024
#define MAX_FILES 100

typedef struct fileStruct {
  char     filename[MAX_FILENAME_LEN];  // filename
  char     id[MAX_FILENAME_LEN];
  fitsfile *fp;        // cfitsio file pointer
  int      fileOpen;   // 0 = file not opened, 1 = file opened
  int      fileType;   // 0 = unknown, 1 = PSRFITS fold, 2 = PSRFITS search, 3 = SDFITS
  float    f0;         // frequency of first channel
  float    f1;         // frequency of last channel
  int      nchan;      // number of channels
  int      nbin;       // number of phase bins in a profile
  int      nsblk;      // number of samples in a block
  int      nbits;      // number of bits
  int      npol;       // number of polarisation products
  float    t0;         // start time of file
  float    nsub;       // number of subintegrations
  float    dt;         // sample time
} fileStruct;


// Function definitions
void createNewInFile(char *fname,fileStruct *infile,int *nInFiles);
void readHeaderInformation(fileStruct *infile, int nInFiles);
void readHeaderInformationPSRFITS_fold(fileStruct *infile,int nInFiles);
void createOutputFile(fileStruct *outfile, int nOutFiles, int outfileNum,fileStruct *infile,int nInFiles);

int main(int argc,char *argv[])
{
  fileStruct *infile;
  fileStruct *outfile;
  int nInFiles=0;
  int nOutFiles=0;
  int fileSet=0;
  int i;
  
  if (!(infile = (fileStruct *)malloc(sizeof(fileStruct)*MAX_FILES)))
    {
      printf("ERROR: Unable to allocate sufficient memory for files\n");
      exit(1);
    }
  if (!(outfile = (fileStruct *)malloc(sizeof(fileStruct)*MAX_FILES)))
    {
      printf("ERROR: Unable to allocate sufficient memory for files\n");
      exit(1);
    }
  // Read the inputs
  for (i=1;i<argc;i++)
    {
      if (strstr(argv[i],"=")!=NULL)
	{
	  // Have a command to set
	  if (strstr(argv[i],"ofileID")!=NULL) // Have an output file name
	    {
	      if (fileSet==1)
		nOutFiles++;
	      
	      fileSet=1;
	      strcpy(outfile[nOutFiles].id,argv[i]+8);
	      outfile[nOutFiles].fileOpen = 0;
	    }
	  else if (strstr(argv[i],"selectfreq")!=NULL)
	    {
	      char temp[1024];
	      char *tok;
	      strcpy(temp,argv[i]+11);
	      tok = strtok(temp,",");
	      sscanf(tok,"%f",&(outfile[nOutFiles].f0));
	      tok = strtok(NULL,",");
	      sscanf(tok,"%f",&(outfile[nOutFiles].f1));
	    }
	}
      else // Have an input filename
	{
	  createNewInFile(argv[i],infile,&nInFiles);
	}
    }
  if (fileSet==1) // Deal with the last output file
    nOutFiles++;
  
  printf("Have loaded %d input files\n",nInFiles);


  // Produce the outputs
  printf("Have %d output files\n",nOutFiles);
  for (i=0;i<nOutFiles;i++)
    {
      createOutputFile(outfile,nOutFiles,i,infile,nInFiles);
    }
  free(infile);
  free(outfile);
}

// Create the output files
void createOutputFile(fileStruct *outfile, int nOutFiles, int outfileNum,fileStruct *infile,int nInFiles)
{
  int useFile;
  int status=0;
  char extendName[1024];
  
  printf("processing %s: f0 = %g, f1 = %g\n",outfile[outfileNum].id,outfile[outfileNum].f0,
	 outfile[outfileNum].f1);

  // Must work out what input files I need and open them as needed
  useFile=0;
  if (infile[useFile].fileOpen == 0)
    {
      sprintf(extendName,"%s[SUBINT][#row > 5 && #row < 20]",infile[useFile].filename);
      fits_open_file(&(infile[useFile].fp),extendName,READONLY,&status);
      fits_report_error(stderr,status);
      infile[useFile].fileOpen = 1;
    }

  // Create the output file
  //  
  sprintf(extendName,"!%s","mynewfile.fits");
  //  sprintf(extendName,"a051015_132035.rf[SUBINT][#row > 3 && #row < 40]");
  printf("Extendname = >%s<\n",extendName);
  fits_create_file(&(outfile[outfileNum].fp),extendName,&status);
  fits_report_error(stderr,status);

  fits_copy_file(infile[useFile].fp,outfile[outfileNum].fp,1,1,1,&status);
  fits_report_error(stderr,status);
  fits_close_file(outfile[outfileNum].fp,&status);
  fits_report_error(stderr,status);

  
  // Should close any files now not needed
  fits_close_file(infile[useFile].fp,&status);
  fits_report_error(stderr,status);
	
}

// Reads in a file name from the command line and initialised the infile structure
void createNewInFile(char *fname,fileStruct *infile,int *nInFiles)
{
  strcpy(infile[*nInFiles].filename,fname);
  infile[*nInFiles].fileOpen = 0; // Will just read the header and then close the file
  infile[*nInFiles].fileType = 1; // PSRFITS fold (read for file extension)

  readHeaderInformation(infile,*nInFiles);
  (*nInFiles)++;
}

void readHeaderInformation(fileStruct *infile, int nInFiles)
{
  if (infile[nInFiles].fileType == 1) // PSRFITS fold mode
    readHeaderInformationPSRFITS_fold(infile,nInFiles);
}

void readHeaderInformationPSRFITS_fold(fileStruct *infile,int nInFiles)
{
  int status=0;
  int colnum=0;
  float nval=0;
  int initFlag=0;
  fitsfile *fp;
  int nsub = 0;
  int npol = 0;
  int nchan = 0;
  int nbin = 0;
  
  printf("Opening file: >%s<\n",infile[nInFiles].filename);
  fits_open_file(&fp,infile[nInFiles].filename,READONLY,&status);
  fits_report_error(stderr,status);

  // Read header information
  fits_movnam_hdu(fp,BINARY_TBL,"SUBINT",1,&status);
  // Number of subintegrations
  fits_read_key(fp,TINT,"NAXIS2",&nsub,NULL,&status);
  fits_read_key(fp,TINT,"NPOL",&npol,NULL,&status);
  fits_read_key(fp,TINT,"NCHAN",&nchan,NULL,&status);
  fits_read_key(fp,TINT,"NBIN",&nbin,NULL,&status);
  printf("NEED TO LOAD FREQUENCY VALUES AND TIME VALUES\n");
  
  printf("Loaded %s nsub = %d npol = %d nchan = %d nbin = %d\n",infile[nInFiles].filename,nsub,npol,nchan,nbin);

  infile[nInFiles].nchan = nchan;
  infile[nInFiles].npol  = npol;
  infile[nInFiles].nsub  = nsub;
  infile[nInFiles].nbin  = nbin;
    
  // Close the file
  fits_close_file(fp,&status);   fits_report_error(stderr,status);
}

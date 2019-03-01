// Present the header information about a PSRFITS file
//
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "pfits.h"

int main(int argc,char *argv[])
{
  char fname[128]="";
  dSet *data;
  int i;
  fitsfile *fp;
  
  for (i=0;i<argc;i++)
    {
      if (strcmp(argv[i],"-f")==0)
	strcpy(fname,argv[++i]);
    }

  data = initialiseDset();
  fp = openFitsFile(fname);
  loadPrimaryHeader(fp,data);
  printf("\n\n Specific header information \n\n");
  displayHeaderInfo(data);
  // Now display all info
  printf("\n\n All primary header information\n\n");
  for (i=0;i<data->phead.nhead;i++)
    printf("%-10.10s %-30.30s %s\n",data->phead.keyname[i],
	   data->phead.val[i],
	   data->phead.comment[i]);
	  

  closeFitsFile(fp);
  freeDset(data);
}

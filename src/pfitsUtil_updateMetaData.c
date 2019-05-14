#include <stdio.h>
#include <string.h>
#include <math.h>
#include <stdlib.h>
#include "fitsio.h"
#include <stdint.h>

#define VERSION 1.0

// gcc -lm -o pfitsUtil_updateMetaData pfitsUtil_updateMetaData.c -lcfitsio
//

int main(int argc,char *argv[])
{
  fitsfile *infptr;
  int status=0;
  char temp[1024];
  char hduStr[1024];
  char str[1024];
  char keyword[1024];
  char oldvalue[1024],comment[1024];
  char newval[1024];
  char card[FLEN_CARD],newcard[FLEN_CARD];
  int keytype;
  int i;
  
  if ( !fits_open_file(&infptr, argv[1], READWRITE, &status) )
    {
      if (status) fits_report_error(stderr, status);

      for (i=2;i<argc;i++)
	{
	  strcpy(str,argv[i]);
	  // Check for : symbol
	  if (strstr(str,":")!=NULL)
	    {
	      strcpy(keyword,strstr(str,":")+1);
	      printf("Have %s\n",keyword);
	      strcpy(strstr(str,":"),"\0");
	      printf("Have %s\n",str);
	      fits_movnam_hdu(infptr,BINARY_TBL,str,0,&status);
	    }
	  else
	    {
	      strcpy(keyword,str);
	      fits_movabs_hdu(infptr, 1, NULL, &status);
	    }
	  
	  strcpy(newval,argv[++i]);
	  printf("Have newval = [%s]\n",newval);

	  
	  if (fits_read_card(infptr,keyword,card,&status))
	      printf("Keyword: %s does not exist\n",keyword);
	  else
	    {
	      fits_parse_value(card,oldvalue,comment,&status);
	      printf("Before changing have: %s\n",card);
	      /* construct template for new keyword */
	      strcpy(newcard, keyword);     /* copy keyword name */
	      strcat(newcard, " = ");       /* '=' value delimiter */
	      if (strstr(newval," ")!=NULL) // Has a space in it
		{
		  sprintf(temp,"'%s'",newval);
		  strcpy(newval,temp);
		}
	      strcat(newcard, newval);     /* new value */
	      if (*comment) {
		strcat(newcard, " / ");  /* comment delimiter */
		strcat(newcard, comment);     /* append the comment */
	      }
	      printf("new card = [%s]\n",newcard);
	      //	      sprintf(newcard,"OBSERVER = 'George Hobbs'                               / Observer names[s]");
	      /* reformat the keyword string to conform to FITS rules */
	      fits_parse_template(newcard, card, &keytype, &status);
	      
	      /* overwrite the keyword with the new value */
	      fits_update_card(infptr, keyword, card, &status);
	      printf("Have changed to: >%s<\n",card);
	    }
	}

      fits_close_file(infptr, &status);
    }
  else
    printf("Unable to open file %s for reading and writing\n",argv[1]);
  

}

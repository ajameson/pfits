// 1) WHERE IS SUBINT:EPOCHS coming from?
// 2) HOW TO ENSURE TABLES IN THE FITS FILE ARE NOT DELETED WHEN OUTPUTTING FROM DSPSR/PSRCHIVE?

//
// Routine to copy a PSRFITS file, but use an up-to-date template
//
// gcc -lm -o pfitsUtil_copyToNew pfitsUtil_copyToNew.c -lcfitsio
// gcc -lm -o pfitsUtil_copyToNew_tycho pfitsUtil_copyToNew.c -L/pulsar/psr/software/20170525/src/util/anaconda2/lib cfitsio/libcfitsio.a -lcurl


#include <stdio.h>
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <fitsio.h>

//void setHeaderInformation(fitsfile *outfptr,char *fname,double *stt_offs_ref);
//void copyData(fitsfile *outfptr,char *fname,double stt_offs_ref);

void processExistingHDU(fitsfile *outfptr,fitsfile *infptr,int hduNum);
void updateObsdescr(fitsfile *outfptr,char *obsdescrFile);

int main(int argc,char *argv[])
{
  char outname[1024];
  char templateFile[1024];
  char fname[1024];
  char constantMD[1024]="UNSET";
  fitsfile *infptr,*outfptr;
  int status=0;
  int statusMovAbs=0;
  char usename[1024];
  int i;
  long numRows;
  int hdu=1;
  double stt_offs_ref=0;
  char keyname[128],val[128],comment[128];
  int updateOBSDESCR=0;
  char obsdescrFile[1024];
  int newHDUtype;
  
  for (i=1;i<argc;i++)
    {
      if (strcmp(argv[i],"-o")==0)
	strcpy(outname,argv[++i]);
      else if (strcmp(argv[i],"-t")==0)
	strcpy(templateFile,argv[++i]);
      else if (strcmp(argv[i],"-constantMD")==0)
	strcpy(constantMD,argv[++i]);
      else if (strcmp(argv[i],"-obsdescr")==0)
	{
	  updateOBSDESCR=1;
	  strcpy(obsdescrFile,argv[++i]);
	}
      else
	strcpy(fname,argv[i]);
    }
  sprintf(usename,"!%s(%s)",outname,templateFile);

  printf("Creating output file\n");
  fits_open_file(&infptr, fname, READONLY, &status);
  
  if ( !fits_create_file(&outfptr, usename, &status) )
    {
      // Go through each HDU in turn and get the suitable data and meta-data
      while( !fits_movabs_hdu(outfptr, hdu, NULL, &status) )
	{
	  printf("UP HERE\n");
	  // Get the name of this HDU
	  if (hdu==1) // PRIMARY
	    {
	      processExistingHDU(outfptr,infptr,hdu);
	    }
	  else 
	    {
	      fits_read_key(outfptr,TSTRING,"EXTNAME",keyname,comment,&status);
	      printf("keyname = %s\n",keyname);
	      // Does this HDU exist in the input file?
	      if (!fits_movnam_hdu(infptr,BINARY_TBL,keyname,0,&status))
		{
		  printf("HDU does exist\n");
		  processExistingHDU(outfptr,infptr,hdu);
		}
	      else
		{
		  printf("HDU does not exist\n");
		  status=0;

		}
	    }
	  printf("Processing HDU %d\n",hdu);
	  hdu++;
	}
      if (status == END_OF_FILE) status = 0;	 	      
    }

  // Now remove tables not of interest

  // Check if there's anything in the polyco table
  fits_movnam_hdu(outfptr,BINARY_TBL,"POLYCO",0,&status);
  fits_get_num_rows(outfptr,&numRows,&status);
  if (numRows == 0)    
    fits_delete_hdu(outfptr, &newHDUtype, &status); 

  // Check if there's anything in the PSRPARAM table
  fits_movnam_hdu(outfptr,BINARY_TBL,"PSRPARAM",0,&status);
  fits_get_num_rows(outfptr,&numRows,&status);
  if (numRows == 0)    
    fits_delete_hdu(outfptr, &newHDUtype, &status); 

  // Check if there's anything in the CAL_POLN table
  fits_movnam_hdu(outfptr,BINARY_TBL,"CAL_POLN",0,&status);
  fits_get_num_rows(outfptr,&numRows,&status);
  if (numRows == 0)    
    fits_delete_hdu(outfptr, &newHDUtype, &status); 

  // Check if there's anything in the FLUX_CAL table
  fits_movnam_hdu(outfptr,BINARY_TBL,"FLUX_CAL",0,&status);
  fits_get_num_rows(outfptr,&numRows,&status);
  if (numRows == 0)    
    fits_delete_hdu(outfptr, &newHDUtype, &status); 

  // Check if there's anything in the FEEDPAR table
  fits_movnam_hdu(outfptr,BINARY_TBL,"FEEDPAR",0,&status);
  fits_get_num_rows(outfptr,&numRows,&status);
  if (numRows == 0)    
    fits_delete_hdu(outfptr, &newHDUtype, &status); 

  // Check if there's anything in the SPECKURT table
  fits_movnam_hdu(outfptr,BINARY_TBL,"SPECKURT",0,&status);
  fits_get_num_rows(outfptr,&numRows,&status);
  if (numRows == 0)    
    fits_delete_hdu(outfptr, &newHDUtype, &status); 

  // Check if there's anything in the DIG_STAT table
  fits_movnam_hdu(outfptr,BINARY_TBL,"DIG_STAT",0,&status);
  fits_get_num_rows(outfptr,&numRows,&status);
  if (numRows == 0)    
    fits_delete_hdu(outfptr, &newHDUtype, &status); 

    // Check if there's anything in the COHDDISP table
  fits_movnam_hdu(outfptr,BINARY_TBL,"COHDDISP",0,&status);
  fits_get_num_rows(outfptr,&numRows,&status);
  if (numRows == 0)    
    fits_delete_hdu(outfptr, &newHDUtype, &status); 
  

  // Check if there's anything in the BANDPASS table
  fits_movnam_hdu(outfptr,BINARY_TBL,"BANDPASS",0,&status);
  fits_get_num_rows(outfptr,&numRows,&status);
  if (numRows == 0)    
    fits_delete_hdu(outfptr, &newHDUtype, &status); 

// Check if there's anything in the DIG_CNTS table
  fits_movnam_hdu(outfptr,BINARY_TBL,"DIG_CNTS",0,&status);
  fits_get_num_rows(outfptr,&numRows,&status);
  if (numRows == 0)    
    fits_delete_hdu(outfptr, &newHDUtype, &status); 

  
  // Now add in extra meta-data information that should remain constant
  // SHOULD RECORD FILE_VER SOMEWHERE
  if (strcmp(constantMD,"UNSET")==0)
    printf("No extra constant meta-data set\n");
  else
    {
      FILE *fin;
      char line[1024];
      char param[1024],val[1024];
      char keyword[1024];
      char oldvalue[1024],comment[1024];
      char newval[1024];
      char card[FLEN_CARD],newcard[FLEN_CARD];
      int keytype;

      fin = fopen(constantMD,"r");
      while (!feof(fin))
	{
	  if (fgets(line,1024,fin)!=NULL)
	    {
	      if (line[0]!='#')
		{
		  sscanf(line,"%s %s",param,val);
		  if (strcmp(param,"FILE_VER")==0)
		    printf("File version %s\n",val);
		  else
		    {
		      if (strstr(param,":")!=NULL)
			{
			  strcpy(keyword,strstr(param,":")+1);
			  strcpy(strstr(param,":"),"\0");
			  printf("Moving to %s\n",param);
			  fits_movnam_hdu(outfptr,BINARY_TBL,param,0,&status);
			}
		      else
			{
			  strcpy(keyword,param);
			  fits_movabs_hdu(outfptr, 1, NULL, &status);
			}
		      if (fits_read_card(outfptr,keyword,card,&status))
			printf("Keyword: >%s< does not exist\n",keyword);
		      else
			{
			  fits_parse_value(card,oldvalue,comment,&status);
			  printf("Before changing have: %s\n",card);
			  /* construct template for new keyword */
			  strcpy(newcard, keyword);     /* copy keyword name */
			  strcat(newcard, " = ");       /* '=' value delimiter */
			  strcat(newcard, val);     /* new value */
			  if (*comment) {
			    strcat(newcard, " / ");  /* comment delimiter */
			    strcat(newcard, comment);     /* append the comment */
			  }
			  /* reformat the keyword string to conform to FITS rules */
			  fits_parse_template(newcard, card, &keytype, &status);
			  
			  /* overwrite the keyword with the new value */
			  fits_update_card(outfptr, keyword, card, &status);
			  printf("Have changed to: >%s<\n",card);
			}
		    }
		  
		}
	    }
	}
      fclose(fin);
    }

  // Update observation description
  
  if (updateOBSDESCR==1)
    updateObsdescr(outfptr,obsdescrFile);

  // Add in meta-data that needs calculating (e.g., STT_LST)
  // STT_LST
  // CHAN_DM
  // BMAJ
  // BMIN
  // BPA
  // TBIN
  
  // Add in meta-data from other lots (e.g., digitiser statistics)
  // PNT_ID
  // SCANLEN
  // CAL_MODE
  
  
  
  
  /* if error occured, print out error message */
  if (status) {printf("AT END\n"); fits_report_error(stderr, status); exit(1);}
  printf("Finishing\n");
  fits_close_file(outfptr,&status);    
  fits_close_file(infptr,&status);
}


void updateObsdescr(fitsfile *outfptr,char *obsdescrFile)
{
  int status=0;
  FILE *fin;
  int colnum;
  int lineNum=1;
  char line[128];
  char *temp = &(line[0]);
  
  fits_movnam_hdu(outfptr,BINARY_TBL,"OBSDESCR",0,&status);
  if (status) {printf("updateObsdescr\n"); fits_report_error(stderr, status); exit(1);}

  if (!(fin = fopen(obsdescrFile,"r")))
    printf("Unable to open file %s\n",obsdescrFile);
  else
    {
      fits_get_colnum(outfptr, CASEINSEN, "DESCR", &colnum, &status);  
      if (status) {printf("updateObsdescr\n"); fits_report_error(stderr, status); exit(1);}      
      while (!feof(fin))
	{
	  if (fgets(line,128,fin)!=NULL)
	    {
	      if (line[strlen(line)-1]='\n') line[strlen(line)-1]='\0';
	      fits_write_col(outfptr,TSTRING,colnum,lineNum++,1,1,&temp,&status);
	      if (status) {printf("updateObsdescr B\n"); fits_report_error(stderr, status); exit(1);}
	    }
	}
    }
  if (status) {printf("updateObsdescr\n"); fits_report_error(stderr, status); exit(1);}      
}

void processExistingHDU(fitsfile *outfptr,fitsfile *infptr,int hduNum)
{
  int status=0;
  int i;
  char keyname[128],val[128],comment[128];
  char outComment[128],outValue[128];
  int keysexist=-1;
  int morekeys=-1;
  char str[1024];
  char cardIn[FLEN_CARD];
  char cardOut[FLEN_CARD];
  char newCard[FLEN_CARD];
  char card[FLEN_CARD];
  char colName[1024];
  char inValue[1024],origValue[1024];
  int keytype;
  int newCardVal;
  int colNumIn;
  int colNumOut;
  long nRows;
  int nCols;
  int sstat=0;
  int naxis;
  int maxdim=5;
  long naxes[maxdim];
  long newLen;
  int j;

  fits_get_hdrspace(outfptr,&keysexist,&morekeys,&status);  
  for (i=0;i<keysexist;i++)
    {
      fits_read_keyn(outfptr,i+1,keyname,val,comment,&status);
      if (strstr(keyname,"XTENSION")==NULL && strstr(keyname,"BITPIX")==NULL &&
	  strstr(keyname,"NAXIS")==NULL && strstr(keyname,"PCOUNT")==NULL &&
	  strstr(keyname,"GCOUNT")==NULL && strstr(keyname,"TFIELDS")==NULL &&
	  strstr(keyname,"TTYPE") == NULL && strstr(keyname,"TFORM")==NULL &&
	  strstr(keyname,"EXTNAME")==NULL && strstr(keyname,"TUNIT")==NULL &&
	  strstr(keyname,"EXTVER")==NULL && strstr(keyname,"COMMENT")==NULL &&
	  strstr(keyname,"HDRVER")==NULL && strstr(keyname,"REF_FREQ")==NULL) 
	{
	  fits_read_card(outfptr,keyname,cardOut,&status);
	  if (!fits_read_card(infptr,keyname,cardIn,&status))
	    {	      
	      fits_parse_value(cardIn,inValue,comment,&status);
	      if (status) {printf("Parse value\n"); fits_report_error(stderr, status); exit(1);}	      
	      fits_parse_value(cardIn,origValue,comment,&status);
	      fits_parse_value(cardOut,outValue,outComment,&status);
	      strcpy(newCard, keyname);
	      strcat(newCard," = ");
	      strcat(newCard, inValue);
	      if (*outComment)
		{
		  strcat(newCard, " / "); 
		  strcat(newCard, outComment);
		}
	      fits_parse_template(newCard, card, &keytype, &status);
	      fits_update_card(outfptr, keyname, card, &status); 
	    }
	  else
	    {
	      printf("Do not have %s in the input file\n",keyname);
	      status=0;
	    }
	}      
    }
  
  if (hduNum > 1)
    {
      // Now copy the data information
      // Get the number of rows from the input data file
      fits_get_num_rows(infptr,&nRows,&status); printf("nRows %d\n",(int)nRows);
      fits_get_num_cols(infptr,&nCols,&status); printf("nCols %d\n",(int)nCols);
      fits_insert_rows(outfptr,0,nRows,&status);
      if (status) {printf("inserting rows\n"); fits_report_error(stderr, status); exit(1);}
      printf("Got here\n");
      for (i=0;i<nCols;i++)
	{
	  fits_get_colname(infptr,0,"*",colName,&colNumIn,&sstat);
	  if (strcmp(colName,"REF_FREQ")==0)
	    {
	      printf("WARNING: IGNORING REF_FREQ\n");
	    }
	  else
	    {
	      fits_get_colnum(outfptr,0,colName,&colNumOut,&status);
	      if (status) {printf("error getting colName\n"); fits_report_error(stderr, status); exit(1);}
	      // Get the size
	      fits_read_tdim(infptr,colNumIn,maxdim,&naxis,naxes,&status);
	      if (status) {printf("have read tdim\n"); fits_report_error(stderr, status); exit(1);}
	      
	      // Must resize column
	      newLen = 1;
	      for (j=0;j<naxis;j++)
		newLen*=naxes[j];
	      
	      fits_modify_vector_len(outfptr,colNumOut,newLen,&status);
	      if (status) {printf("modify vector length error: \n"); fits_report_error(stderr, status); exit(1);}
	      if (strcmp(colName,"DATA")==0)
		{
		  char tdimStr[128];
		  sprintf(tdimStr,"TDIM%d",colNumOut);
		  fits_delete_key(outfptr, tdimStr, &status);
		  fits_write_tdim(outfptr,colNumOut,naxis,naxes,&status);
		}
	      
	      //      fits_write_tdim(outfptr,colNumOut,naxis,naxes,&status);
	      if (status) {printf("write tdim %d %d (%s)\n",naxis,naxes[0],colName); fits_report_error(stderr, status); exit(1);}
	      printf("Column %s (in = %d) (out = %d) size = %d\n",colName,colNumIn,colNumOut,naxis);
	      fits_copy_col(infptr,outfptr,colNumIn,colNumOut,0,&status);
	      
	      //fits_copy_col(infptr,outfptr,i+1,i+1,0,&status);
	    }
	}
      if (status) {printf("copying data\n"); fits_report_error(stderr, status); exit(1);}
    }


}

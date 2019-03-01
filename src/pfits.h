// Header information for pfits software

#include "fitsio.h"

typedef struct pheader {
  int nhead;
  //
  // All header information
  char **keyname;
  char **val;
  char **comment;
  //
  // Most useful header information
  int nchan; // Number of channels
  int nbits; // Number of bits
  int nsamp; // Number of samples
  int nsub;  // Number of subintegrations
  int nsblk; // Number of samples per subintegration
  int nbin;  // Number of bins
  int npol;  // Number of polarisations
  float tsamp; // Sample time
  float freq;  // Centre frequency for observation
  float bw;    // Observation bandwidth
  float zeroOff; // Zero offsets
  int   imjd;  // Integer start time (day)
  float smjd;  // Fraction start time (sec)
  float stt_offs; // Start time offset (seconds)
  float chanbw;   // Channel bandwidth (MHz)
  float dm;       // Pulsar's dispersion measure (cm^-3 pc)
  float period;   // Pulsar's period (s)   
  char obsMode[128]; // PSR, CAL, SEARCH
  char source[128];
} pheader;

typedef struct dSet {
  int pheaderSet;
  int subintTable;
  int psrparamTable;
  pheader phead;
} dSet;

fitsfile * openFitsFile(char *fname);
void closeFitsFile(fitsfile *fp);
//void initialiseDset(dSet *data);
dSet * initialiseDset();
void freeDset(dSet *data);
void readPhead(dSet *data);
void loadPrimaryHeader(fitsfile *fp,dSet *data);
int readBandpass(fitsfile *fp,float *bpass);
int extractPolData(fitsfile *fp,dSet *data,int pol,float *arr,float t1,float t2);
int extractDataZeroDM(fitsfile *fp,dSet *data,long s1,long s2,float *mean,float *min,float *max,long maxVal,int *nsmooth);
int extractData(fitsfile *fp,dSet *data,long s1,long s2,float *mean,float *min,float *max,long maxVal,int *nsmooth,float dm);
int extractFoldData(fitsfile *fp,dSet *data,float dm,float *fx,float *fy,float *freq_y,float *time_y,float* bpass,int sub0);
void findPosSample(dSet *data,long s,long *sub,long *samp);
void bytesToValues(int samplesperbyte, int nchan, unsigned char *val2, 
		   unsigned char *values);
void eightBits(int eight_bit_number, unsigned char *results, int *index);
void fourBits(int eight_bit_number, unsigned char *results, int *index);
void twoBits(int eight_bit_number, unsigned char *results, int *index);
void oneBit(int eight_bit_number, unsigned char *results, int *index);
int andCount(unsigned char *cval,unsigned char *mask,int n);
int foldDM(fitsfile *fp,dSet *data,long s1,long s2,float ***fold,float *dm,double *period,int ndm,int nperiod,int nbin);
int nint(float val);
int writeData(fitsfile *fp,FILE *fout,dSet *data,long s1,long s2,float dm);
int createWavFile(fitsfile *fp,FILE *fout,dSet *data,long s1,long s2,float dm);


void getbaseline(float* prof, int nbin, float fwindow,float* mean, float *rms);

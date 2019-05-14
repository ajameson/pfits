
noinst_LTLIBRARIES = libpfits.la

libpfits_la_headers = pfits.h

libpfits_la_SOURCES = pfits.c

bin_PROGRAMS = fitsProcess \
  pfits_fv \
  pfits_bpass \
  pfitsUtil_copyToNew \
  pfitsUtil_headerInformation \
  pfitsUtil_foldmode_combineFreq \
  pfitsUtil_foldmode_combineTime \
  pfitsUtil_foldmode_extractFreq \
  pfitsUtil_foldmode_extractTime \
  pfitsUtil_foldmode_useTemplate \
  pfitsUtil_searchmode_combineFreq \
  pfitsUtil_searchmode_combineTime \
  pfitsUtil_searchMode_extractPart \
  pfitsUtil_updateMetaData

LDADD = libpfits.la @CFITSIO_LIBS@

if HAVE_PGPLOT

bin_PROGRAMS += pfits_plotfold pfits_plotsearch pfits_plot pfits_bpass

AM_CFLAGS = @PGPLOT_CFLAGS@
LDADD += @PGPLOT_LIBS@

endif

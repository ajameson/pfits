
noinst_LTLIBRARIES = libpfits.la

libpfits_la_headers = pfits.h

libpfits_la_SOURCES = pfits.c

bin_PROGRAMS = pfits_fv pfits_bpass

LDADD = libpfits.la @CFITSIO_LIBS@

if HAVE_PGPLOT

bin_PROGRAMS += pfits_plotfold pfits_plotsearch pfits_plot pfits_bpass

AM_CFLAGS = @PGPLOT_CFLAGS@
LDADD += @PGPLOT_LIBS@

endif

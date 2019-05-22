#                                               -*- Autoconf -*-
# This files generates the configure file.
#
# If you change this file, you usually should run
# autoreconf -fvi
AC_PREREQ(2.61)

# This is where the version information is set
AC_INIT([PFITS], [0.1], [andrew.jameson@gmail.com])

AC_CONFIG_MACRO_DIR(config)

AC_CONFIG_SRCDIR([src/pfits.c])
AC_DISABLE_SHARED

# Automake initialization
AM_INIT_AUTOMAKE([-Wall -Werror foreign])
AM_CONFIG_HEADER([config.h])

# Checks for programs
m4_ifdef([AM_PROG_AR], [AM_PROG_AR])
AC_PROG_CC
AC_PROG_F77
AC_PROG_LIBTOOL

# Checks for libraries.
SWIN_LIB_CFITSIO
SWIN_LIB_PGPLOT

AC_CONFIG_FILES([
                 Makefile
                 src/Makefile
                ])
AC_OUTPUT


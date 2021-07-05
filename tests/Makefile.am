## Process this file with automake to produce Makefile.in

AM_CPPFLAGS = -I$(top_srcdir) -I$(top_srcdir)/src -I$(includedir)
LDADD = $(top_builddir)/src/.libs/libmop.a 

DEPS = $(top_builddir)/src/libmop.la

bin_PROGRAMS = \
	mop-test

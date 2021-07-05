## Process this file with automake to produce Makefile.in

AM_CPPFLAGS = -I$(top_srcdir) -I$(includedir) -DG_LOG_DOMAIN=\"MOP\"

bin_SCRIPTS=mop-config

BUILT_SOURCES= \
	mop-config

mop-config: mop-config.in

lib_LTLIBRARIES = libmop.la

libmop_la_LDFLAGS = -version-info $(LT_CURRENT):$(LT_REVISION):$(LT_AGE)\
	-release $(LT_RELEASE) -export-dynamic

libmop_la_SOURCES = \
	mop.c \
	interp.c \
	block.c \
	mopblock.h \
	mop-logging.c \
	mop-logging.h \
	mop.h

include_HEADERS = \
	mop.h

CLEANFILES = $(BUILT_SOURCES)

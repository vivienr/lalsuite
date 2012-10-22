# Makefile for SMEE

AM_CPPFLAGS = -I$(top_srcdir)/src -I$(top_srcdir)/src/lalapps
LDADD = $(top_builddir)/src/lalapps/liblalapps.la

if LALINFERENCE
smee_prog = LALSMEE
endif

bin_PROGRAMS = $(smee_prog)

LALSMEE_SOURCES = \
	LALSMEE.c \
	SMEELikelihood.c \
	SMEELikelihood.h \
	SMEEPrior.c \
	SMEEPrior.h \
	SMEEReadData.c \
	SMEEReadData.h 

## make sure lalapps.la is rebuilt when appropriate
BUILT_SOURCES =
include $(top_srcdir)/misc/vcsID.common

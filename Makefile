PROG=		bcftools

all: $(PROG)

# Adjust $(HTSDIR) to point to your top-level htslib directory
HTSDIR = ../htslib
include $(HTSDIR)/htslib.mk
HTSLIB = $(HTSDIR)/libhts.a

CC=			gcc
CFLAGS=		-g -Wall -Wc++-compat -O2
DFLAGS=
OBJS=		main.o bcfidx.o tabix.o \
			vcfstats.o vcfisec.o vcfmerge.o vcfquery.o vcffilter.o filter.o vcfsom.o \
            vcfnorm.o vcfgtcheck.o vcfsubset.o vcfannotate.o vcfroh.o \
            vcfcall.o mcall.o vcmp.o \
            ccall.o em.o prob1.o kmin.o # the original samtools calling
INCLUDES=	-I. -I$(HTSDIR)

prefix      = /usr/local
exec_prefix = $(prefix)
bindir      = $(exec_prefix)/bin
mandir      = $(prefix)/share/man
man1dir     = $(mandir)/man1

INSTALL = install -p
INSTALL_PROGRAM = $(INSTALL)
INSTALL_DATA    = $(INSTALL) -m 644


all:$(PROG)

# See htslib/Makefile
PACKAGE_VERSION  = 0.0.1
ifneq "$(wildcard .git)" ""
PACKAGE_VERSION := $(shell git describe --always --dirty)
version.h: $(if $(wildcard version.h),$(if $(findstring "$(PACKAGE_VERSION)",$(shell cat version.h)),,force))
endif
version.h:
	echo '#define BCFTOOLS_VERSION "$(PACKAGE_VERSION)"' > $@


.SUFFIXES:.c .o
.PHONY:all install lib test force

force:

.c.o:
		$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) $< -o $@

test: $(PROG)
		./test/test.pl

main.o: version.h bcftools.h
vcfcall.o: vcfcall.c call.h mcall.c prob1.h $(HTSDIR)/htslib/kfunc.h $(HTSDIR)/htslib/vcf.h
mcall.o ccall.o: call.h vcmp.h
vcffilter.o: filter.h
vcfsubset.o: filter.h
vcfnorm.o: rbuf.h
vcffilter.o: rbuf.h
vcfroh.o: rbuf.h

bcftools: $(HTSLIB) $(OBJS)
		$(CC) $(CFLAGS) -o $@ $(OBJS) $(HTSLIB) -lpthread -lz -lm


install: $(PROG)
		mkdir -p $(DESTDIR)$(bindir) $(DESTDIR)$(man1dir)
		$(INSTALL_PROGRAM) $(PROG) plot-vcfstats $(DESTDIR)$(bindir)
		$(INSTALL_DATA) bcftools.1 $(DESTDIR)$(man1dir)


cleanlocal:
		rm -fr gmon.out *.o a.out *.dSYM *~ $(PROG) version.h

clean:cleanlocal clean-htslib

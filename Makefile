PROG=		bcftools
TEST_PROG=  test/test-rbuf


all: $(PROG) $(TEST_PROG)

# Adjust $(HTSDIR) to point to your top-level htslib directory
HTSDIR = ../htslib
include $(HTSDIR)/htslib.mk
HTSLIB = $(HTSDIR)/libhts.a

CC=			gcc
CFLAGS=		-g -Wall -Wc++-compat -O2
DFLAGS=
OBJS=		main.o bcfidx.o tabix.o \
			vcfstats.o vcfisec.o vcfmerge.o vcfquery.o vcffilter.o filter.o vcfsom.o \
            vcfnorm.o vcfgtcheck.o vcfview.o vcfannotate.o vcfroh.o \
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


all:$(PROG) plugins

# See htslib/Makefile
PACKAGE_VERSION  = 0.0.1
ifneq "$(wildcard .git)" ""
PACKAGE_VERSION := $(shell git describe --always --dirty)
version.h: $(if $(wildcard version.h),$(if $(findstring "$(PACKAGE_VERSION)",$(shell cat version.h)),,force))
endif
version.h:
	echo '#define BCFTOOLS_VERSION "$(PACKAGE_VERSION)"' > $@


.SUFFIXES:.c .o
.PHONY:all install lib test force plugins

force:

.c.o:
		$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) $< -o $@

test: $(PROG) plugins test/test-rbuf
		./test/test.pl

PLUGINC = $(foreach dir, plugins, $(wildcard $(dir)/*.c))
PLUGINS = $(PLUGINC:.c=.so)

plugins: $(PLUGINS)

%.so: %.c
	$(CC) $(CFLAGS) $(INCLUDES) -fPIC -shared -o $@ $<

main.o: version.h $(HTSDIR)/version.h bcftools.h
vcfcall.o: vcfcall.c call.h mcall.c prob1.h $(HTSDIR)/htslib/kfunc.h $(HTSDIR)/htslib/vcf.h bcftools.h
mcall.o ccall.o: call.h vcmp.h bcftools.h
vcffilter.o: bcftools.h filter.h
vcfsubset.o: bcftools.h filter.h
vcfnorm.o: bcftools.h rbuf.h
vcffilter.o: bcftools.h rbuf.h
vcfroh.o: bcftools.h rbuf.h
vcfannotate.o: bcftools.h vcmp.h $(HTSDIR)/htslib/kseq.h
test/test-rbuf.o: rbuf.h test/test-rbuf.c

test/test-rbuf: test/test-rbuf.o
		$(CC) $(CFLAGS) -o $@ -lm -ldl $<

bcftools: $(HTSLIB) $(OBJS)
		$(CC) $(CFLAGS) -o $@ $(OBJS) $(HTSLIB) -lpthread -lz -lm -ldl


install: $(PROG)
		mkdir -p $(DESTDIR)$(bindir) $(DESTDIR)$(man1dir)
		$(INSTALL_PROGRAM) $(PROG) plot-vcfstats $(DESTDIR)$(bindir)
		$(INSTALL_DATA) bcftools.1 $(DESTDIR)$(man1dir)

cleanlocal: cleantest
		rm -fr gmon.out *.o a.out *.dSYM *~ $(PROG) version.h plugins/*.so

cleantest:
		rm -fr test/*.o test/*~ $(TEST_PROG)

clean:cleanlocal clean-htslib

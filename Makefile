# Makefile for bcftools, utilities for Variant Call Format VCF/BCF files.
#
#   Copyright (C) 2012-2014 Genome Research Ltd.
#
#   Author: Petr Danecek <pd3@sanger.ac.uk>
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
# THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
# DEALINGS IN THE SOFTWARE.

PROG=       bcftools
TEST_PROG=  test/test-rbuf


all: $(PROG) $(TEST_PROG)

# Adjust $(HTSDIR) to point to your top-level htslib directory
HTSDIR = ../htslib
include $(HTSDIR)/htslib.mk
HTSLIB = $(HTSDIR)/libhts.a
BGZIP  = $(HTSDIR)/bgzip
TABIX  = $(HTSDIR)/tabix

CC       = gcc
CPPFLAGS =
CFLAGS   = -g -Wall -Wc++-compat -O2
LDFLAGS  =
LIBS     =

OBJS     = main.o vcfindex.o tabix.o \
           vcfstats.o vcfisec.o vcfmerge.o vcfquery.o vcffilter.o filter.o vcfsom.o \
           vcfnorm.o vcfgtcheck.o vcfview.o vcfannotate.o vcfroh.o vcfconcat.o \
           vcfcall.o mcall.o vcmp.o gvcf.o reheader.o convert.o vcfconvert.o tsv2vcf.o \
           vcfcnv.o HMM.o vcfplugin.o consensus.o ploidy.o version.o \
           ccall.o em.o prob1.o kmin.o # the original samtools calling

EXTRA_CPPFLAGS = -I. -I$(HTSDIR) -DPLUGINPATH=\"$(pluginpath)\"
GSL_LIBS       =

# The polysomy command is not compiled by default because it brings dependency
# on libgsl. The command can be compiled wth `make USE_GPL=1`. See the INSTALL
# and LICENSE documents to understand license implications.
ifdef USE_GPL
    EXTRA_CPPFLAGS += -DUSE_GPL
    OBJS += polysomy.o
    GSL_LIBS = -lgsl -lcblas
endif

prefix      = /usr/local
exec_prefix = $(prefix)
bindir      = $(exec_prefix)/bin
libdir      = $(exec_prefix)/lib
libexecdir  = $(exec_prefix)/libexec
mandir      = $(prefix)/share/man
man1dir     = $(mandir)/man1

plugindir   = $(libexecdir)/bcftools
pluginpath  = $(plugindir)

MKDIR_P = mkdir -p
INSTALL = install -p
INSTALL_PROGRAM = $(INSTALL)
INSTALL_DATA    = $(INSTALL) -m 644
INSTALL_DIR     = $(MKDIR_P) -m 755


all:$(PROG) plugins

# See htslib/Makefile
PACKAGE_VERSION = 1.2
ifneq "$(wildcard .git)" ""
PACKAGE_VERSION := $(shell git describe --always --dirty)
DOC_VERSION :=  $(shell git describe --always)+
DOC_DATE := $(shell date +'%Y-%m-%d %R %Z')
version.h: $(if $(wildcard version.h),$(if $(findstring "$(PACKAGE_VERSION)",$(shell cat version.h)),,force))
endif
version.h:
	echo '#define BCFTOOLS_VERSION "$(PACKAGE_VERSION)"' > $@


.SUFFIXES:.c .o
.PHONY:all clean clean-all clean-plugins distclean install lib tags test testclean force plugins docs

force:

.c.o:
	$(CC) $(CFLAGS) $(EXTRA_CPPFLAGS) $(CPPFLAGS) -c -o $@ $<

test: $(PROG) plugins test/test-rbuf $(BGZIP) $(TABIX)
	./test/test.pl --exec bgzip=$(BGZIP) --exec tabix=$(TABIX)

test-plugins: $(PROG) plugins test/test-rbuf $(BGZIP) $(TABIX)
	./test/test.pl --plugins --exec bgzip=$(BGZIP) --exec tabix=$(TABIX)


# Plugin rules
PLUGINC = $(foreach dir, plugins, $(wildcard $(dir)/*.c))
PLUGINS = $(PLUGINC:.c=.so)
PLUGINM = $(PLUGINC:.c=.mk)

%.so: %.c version.h version.c $(HTSDIR)/libhts.so
	$(CC) -fPIC -shared $(CFLAGS) $(EXTRA_CPPFLAGS) $(CPPFLAGS) -L$(HTSDIR) $(LDFLAGS) -o $@ version.c $< -lhts $(LIBS)

-include $(PLUGINM)

plugins: $(PLUGINS)


bcftools_h = bcftools.h $(htslib_vcf_h)
call_h = call.h $(htslib_vcf_h) $(htslib_synced_bcf_reader_h) vcmp.h
convert_h = convert.h $(htslib_vcf_h)
tsv2vcf_h = tsv2vcf.h $(htslib_vcf_h)
filter_h = filter.h $(htslib_vcf_h)
prob1_h = prob1.h $(htslib_vcf_h) $(call_h)
roh_h = HMM.h $(htslib_vcf_h) $(htslib_synced_bcf_reader_h) $(HTSDIR)/htslib/kstring.h $(HTSDIR)/htslib/kseq.h $(bcftools_h)
cnv_h = HMM.h $(htslib_vcf_h) $(htslib_synced_bcf_reader_h)

main.o: main.c $(htslib_hts_h) version.h $(bcftools_h)
vcfannotate.o: vcfannotate.c $(htslib_vcf_h) $(htslib_synced_bcf_reader_h) $(HTSDIR)/htslib/kseq.h $(bcftools_h) vcmp.h $(filter_h)
vcfplugin.o: vcfplugin.c $(htslib_vcf_h) $(htslib_synced_bcf_reader_h) $(HTSDIR)/htslib/kseq.h $(bcftools_h) vcmp.h $(filter_h)
vcfcall.o: vcfcall.c $(htslib_vcf_h) $(HTSDIR)/htslib/kfunc.h $(htslib_synced_bcf_reader_h) $(bcftools_h) $(call_h) $(prob1_h)
vcfconcat.o: vcfconcat.c $(htslib_vcf_h) $(htslib_synced_bcf_reader_h) $(HTSDIR)/htslib/kseq.h $(bcftools_h)
vcfconvert.o: vcfconvert.c $(htslib_vcf_h) $(htslib_bgzf_h) $(htslib_synced_bcf_reader_h) $(htslib_vcfutils_h) $(bcftools_h) $(filter_h) $(convert_h) $(tsv2vcf_h)
vcffilter.o: vcffilter.c $(htslib_vcf_h) $(htslib_synced_bcf_reader_h) $(htslib_vcfutils_h) $(bcftools_h) $(filter_h) rbuf.h
vcfgtcheck.o: vcfgtcheck.c $(htslib_vcf_h) $(htslib_synced_bcf_reader_h) $(htslib_vcfutils_h) $(bcftools_h)
vcfindex.o: vcfindex.c $(htslib_vcf_h) $(htslib_tbx_h)
vcfisec.o: vcfisec.c $(htslib_vcf_h) $(htslib_synced_bcf_reader_h) $(htslib_vcfutils_h) $(bcftools_h) $(filter_h)
vcfmerge.o: vcfmerge.c $(htslib_vcf_h) $(htslib_synced_bcf_reader_h) $(htslib_vcfutils_h) $(bcftools_h) vcmp.h $(HTSDIR)/htslib/khash.h
vcfnorm.o: vcfnorm.c $(htslib_vcf_h) $(htslib_synced_bcf_reader_h) $(htslib_faidx_h) $(bcftools_h) rbuf.h
vcfquery.o: vcfquery.c $(htslib_vcf_h) $(htslib_synced_bcf_reader_h) $(htslib_vcfutils_h) $(bcftools_h) $(filter_h) $(convert_h)
vcfroh.o: vcfroh.c $(roh_h)
vcfcnv.o: vcfcnv.c $(cnv_h)
vcfsom.o: vcfsom.c $(htslib_vcf_h) $(htslib_synced_bcf_reader_h) $(htslib_vcfutils_h) $(bcftools_h)
vcfstats.o: vcfstats.c $(htslib_vcf_h) $(htslib_synced_bcf_reader_h) $(htslib_vcfutils_h) $(htslib_faidx_h) $(bcftools_h)
vcfview.o: vcfview.c $(htslib_vcf_h) $(htslib_synced_bcf_reader_h) $(htslib_vcfutils_h) $(bcftools_h) $(filter_h)
reheader.o: reheader.c $(htslib_vcf_h) $(htslib_bgzf_h) $(HTSDIR)/htslib/kseq.h $(bcftools_h)
tabix.o: tabix.c $(htslib_bgzf_h) $(htslib_tbx_h)
ccall.o: ccall.c $(HTSDIR)/htslib/kfunc.h $(call_h) kmin.h $(prob1_h)
convert.o: convert.c $(htslib_vcf_h) $(htslib_synced_bcf_reader_h) $(htslib_vcfutils_h) $(bcftools_h) $(convert_h)
tsv2vcf.o: tsv2vcf.c $(tsv2vcf_h)
em.o: em.c $(htslib_vcf_h) kmin.h $(call_h)
filter.o: filter.c $(HTSDIR)/htslib/khash_str2int.h $(filter_h) $(bcftools_h) $(htslib_hts_defs_h) $(htslib_vcfutils_h)
gvcf.o: gvcf.c $(call_h)
kmin.o: kmin.c kmin.h
mcall.o: mcall.c $(HTSDIR)/htslib/kfunc.h $(call_h)
prob1.o: prob1.c $(prob1_h)
vcmp.o: vcmp.c $(htslib_hts_h) vcmp.h
polysomy.o: polysomy.c $(htslib_hts_h)
consensus.o: consensus.c $(htslib_hts_h) $(HTSDIR)/htslib/kseq.h rbuf.h $(bcftools_h) $(HTSDIR)/htslib/regidx.h
version.o: version.h version.c

test/test-rbuf.o: test/test-rbuf.c rbuf.h

test/test-rbuf: test/test-rbuf.o
	$(CC) $(LDFLAGS) -o $@ $^ -lm -ldl $(LIBS)

bcftools: $(HTSLIB) $(OBJS)
	$(CC) $(LDFLAGS) -o $@ $(OBJS) $(HTSLIB) -lpthread -lz -lm -ldl $(GSL_LIBS) $(LIBS)

doc/bcftools.1: doc/bcftools.txt
	cd doc && a2x -adate="$(DOC_DATE)" -aversion=$(DOC_VERSION) --doctype manpage --format manpage bcftools.txt

doc/bcftools.html: doc/bcftools.txt
	cd doc && a2x -adate="$(DOC_DATE)" -aversion=$(DOC_VERSION) --doctype manpage --format xhtml bcftools.txt

docs: doc/bcftools.1 doc/bcftools.html

install: $(PROG) doc/bcftools.1
	$(INSTALL_DIR) $(DESTDIR)$(bindir) $(DESTDIR)$(man1dir) $(DESTDIR)$(plugindir)
	$(INSTALL_PROGRAM) $(PROG) plot-vcfstats vcfutils.pl $(DESTDIR)$(bindir)
	$(INSTALL_DATA) doc/bcftools.1 $(DESTDIR)$(man1dir)
	$(INSTALL_PROGRAM) plugins/*.so $(DESTDIR)$(plugindir)

clean: testclean clean-plugins
	-rm -f gmon.out *.o *~ $(PROG) version.h plugins/*.so plugins/*.P
	-rm -rf *.dSYM plugins/*.dSYM test/*.dSYM

clean-plugins:
	-rm -f plugins/*.so plugins/*.P
	-rm -rf plugins/*.dSYM

testclean:
	-rm -f test/*.o test/*~ $(TEST_PROG)

distclean: clean
	-rm -f TAGS

clean-all: clean clean-htslib

tags:
	ctags -f TAGS *.[ch] plugins/*.[ch]

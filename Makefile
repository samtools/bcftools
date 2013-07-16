# The default version string in bam.h and bcftools/bcf.h can be overriden directly
#   make VERSION="-DVERSION='\\\"my-version\\\"'"
# or using the git-stamp rule
#   make git-stamp
VERSION=

# Adjust $(HTSDIR) to point to your top-level htslib directory
HTSDIR = ../htslib
HTSLIB = $(HTSDIR)/libhts.a

CC=			gcc
CFLAGS=		-g -Wall -Wc++-compat -O2 $(VERSION)
DFLAGS=
OBJS=		main.o vcfview.o bcfidx.o tabix.o \
			vcfcheck.o vcfisec.o vcfmerge.o vcfquery.o vcffilter.o \
            vcfnorm.o vcfgtcheck.o vcfsubset.o
INCLUDES=	-I. -I$(HTSDIR)
PROG=		bcftools

.SUFFIXES:.c .o
.PHONY:all lib test

.c.o:
		$(CC) -c $(CFLAGS) $(DFLAGS) $(INCLUDES) $< -o $@

all:$(PROG)

git-stamp:
		make VERSION="-DVERSION='\"`git describe --always --dirty`\"'"

test:
		./test/test.pl

bcftools:$(OBJS)
		$(CC) $(CFLAGS) -o $@ $(OBJS) $(HTSLIB) -lpthread -lz -lm

clean:
		rm -fr gmon.out *.o a.out *.dSYM *~ $(PROG)

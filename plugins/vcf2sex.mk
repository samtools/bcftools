plugins/vcf2sex.so: plugins/vcf2sex.c version.h version.c ploidy.h ploidy.c $(HTSDIR)/libhts.so
	$(CC) $(CFLAGS) $(INCLUDES) -fPIC -shared -o $@ ploidy.c version.c $< -L$(HTSDIR) -lhts

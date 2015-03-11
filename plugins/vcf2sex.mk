plugins/vcf2sex.so: plugins/vcf2sex.c version.h version.c ploidy.h ploidy.c $(HTSDIR)/libhts.so
	$(CC) -fPIC -shared $(CFLAGS) $(EXTRA_CPPFLAGS) $(CPPFLAGS) -L$(HTSDIR) $(LDFLAGS) -o $@ ploidy.c version.c $< -lhts $(LIBS)

plugins/fixploidy.so: plugins/fixploidy.c version.h version.c ploidy.h ploidy.c $(HTSDIR)/libhts.so
	$(CC) $(CFLAGS) $(INCLUDES) -fPIC -shared -o $@ ploidy.c version.c $< -L$(HTSDIR) -lhts

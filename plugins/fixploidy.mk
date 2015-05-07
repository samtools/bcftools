plugins/fixploidy.so: plugins/fixploidy.c version.h version.c ploidy.h ploidy.c $(htslib_shared)
	$(CC) -fPIC -shared $(CFLAGS) $(EXTRA_CPPFLAGS) $(CPPFLAGS) -L$(HTSDIR) $(LDFLAGS) -o $@ ploidy.c version.c $< -lhts $(LIBS)

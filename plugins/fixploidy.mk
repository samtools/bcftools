plugins/fixploidy.so: plugins/fixploidy.c version.h version.c ploidy.h ploidy.c
	$(CC) $(CFLAGS) $(INCLUDES) -fPIC -shared -o $@ ploidy.c version.c $< -L$(HTSDIR) -lhts

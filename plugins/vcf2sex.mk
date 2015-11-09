plugins/vcf2sex.so: plugins/vcf2sex.c version.h version.c ploidy.h ploidy.c
	$(CC) $(PLUGIN_FLAGS) $(CFLAGS) $(EXTRA_CPPFLAGS) $(CPPFLAGS) $(LDFLAGS) -o $@ ploidy.c version.c $< $(LIBS)

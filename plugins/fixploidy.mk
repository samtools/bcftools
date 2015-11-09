plugins/fixploidy.so: plugins/fixploidy.c version.h version.c ploidy.h ploidy.c
	$(CC) $(PLUGIN_FLAGS) $(CFLAGS) $(EXTRA_CPPFLAGS) $(CPPFLAGS) $(LDFLAGS) -o $@ ploidy.c version.c $< $(LIBS)

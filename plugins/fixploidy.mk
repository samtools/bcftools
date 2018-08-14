plugins/fixploidy.so: plugins/fixploidy.c version.h version.c ploidy.h ploidy.c
	$(CC) $(PLUGIN_FLAGS) $(CFLAGS) $(ALL_CPPFLAGS) $(EXTRA_CPPFLAGS) $(LDFLAGS) -o $@ ploidy.c version.c $< $(LIBS)

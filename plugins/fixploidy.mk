plugins/fixploidy.so: plugins/fixploidy.c version.h version.c regidx.h regidx.c ploidy.h ploidy.c
	$(CC) $(PLUGIN_FLAGS) $(CFLAGS) $(ALL_CPPFLAGS) $(EXTRA_CPPFLAGS) $(LDFLAGS) -o $@ ploidy.c regidx.c version.c $< $(LIBS)

plugins/af-dist.so: plugins/af-dist.c version.h version.c bin.h bin.c
	$(CC) $(PLUGIN_FLAGS) $(CFLAGS) $(ALL_CPPFLAGS) $(EXTRA_CPPFLAGS) $(LDFLAGS) -o $@ bin.c version.c $< $(LIBS)

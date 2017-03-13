plugins/af-dist.so: plugins/af-dist.c version.h version.c bin.h bin.c
	$(CC) $(PLUGIN_FLAGS) $(CFLAGS) $(EXTRA_CPPFLAGS) $(CPPFLAGS) $(LDFLAGS) -o $@ bin.c version.c $< $(LIBS)

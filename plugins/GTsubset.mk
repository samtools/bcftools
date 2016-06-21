plugins/GTsubset.so: plugins/GTsubset.c version.h version.c 
	$(CC) $(PLUGIN_FLAGS) $(CFLAGS) $(EXTRA_CPPFLAGS) $(CPPFLAGS) $(LDFLAGS) -o $@ version.c $< $(LIBS)

plugins/GTisec.so: plugins/GTisec.c version.h version.c 
	$(CC) $(PLUGIN_FLAGS) $(CFLAGS) $(EXTRA_CPPFLAGS) $(CPPFLAGS) $(LDFLAGS) -o $@ version.c $< $(LIBS)

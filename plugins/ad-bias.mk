plugins/ad-bias.so: plugins/ad-bias.c version.h version.c convert.h convert.c
	$(CC) $(PLUGIN_FLAGS) $(CFLAGS) $(EXTRA_CPPFLAGS) $(CPPFLAGS) $(LDFLAGS) -o $@ convert.c version.c $< $(LIBS)

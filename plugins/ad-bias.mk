plugins/ad-bias.so: plugins/ad-bias.c version.h version.c convert.h convert.c
	$(CC) $(PLUGIN_FLAGS) $(CFLAGS) $(ALL_CPPFLAGS) $(EXTRA_CPPFLAGS) $(LDFLAGS) -o $@ convert.c version.c $< $(PLUGIN_LIBS) $(LIBS)

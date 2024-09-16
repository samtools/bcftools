plugins/noise-profile.so: plugins/noise-profile.c version.h regidx.h mpileup2/mpileup.o
	$(CC) $(PLUGIN_FLAGS) $(CFLAGS) $(ALL_CPPFLAGS) $(EXTRA_CPPFLAGS) $(LDFLAGS) -o $@ version.c $< $(PLUGIN_LIBS) $(LIBS)

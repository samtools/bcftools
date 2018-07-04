plugins/color-chrs.so: plugins/color-chrs.c version.h version.c HMM.h HMM.c
	$(CC) $(PLUGIN_FLAGS) $(CFLAGS) $(ALL_CPPFLAGS) $(EXTRA_CPPFLAGS) $(LDFLAGS) -o $@ HMM.c version.c $< $(LIBS)

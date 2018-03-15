plugins/setGT.so: plugins/setGT.c version.h version.c filter.h filter.c
	$(CC) $(PLUGIN_FLAGS) $(CFLAGS) $(ALL_CPPFLAGS) $(EXTRA_CPPFLAGS) $(PERL_CFLAGS) $(LDFLAGS) -o $@ filter.c version.c $< $(LIBS)

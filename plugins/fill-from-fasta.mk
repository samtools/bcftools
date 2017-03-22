plugins/fill-from-fasta.so: plugins/fill-from-fasta.c version.h version.c filter.h filter.c
	$(CC) $(PLUGIN_FLAGS) $(CFLAGS) $(ALL_CPPFLAGS) $(EXTRA_CPPFLAGS) $(LDFLAGS) -o $@ filter.c version.c $< $(LIBS)

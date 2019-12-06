plugins/isecGT.so: plugins/isecGT.c version.h version.c smpl_ilist.h
	$(CC) $(PLUGIN_FLAGS) $(CFLAGS) $(ALL_CPPFLAGS) $(EXTRA_CPPFLAGS) $(LDFLAGS) -o $@ smpl_ilist.c version.c $< $(PLUGIN_LIBS) $(LIBS)

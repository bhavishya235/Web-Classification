.c.o:
	$(CC) $(ALL_CPPFLAGS) $(ALL_CFLAGS) -c $< -o $@
.c.d:
	@echo "making $@ from $<"
	@$(CC) -MM $(ALL_CPPFLAGS) $< > $@
.m.d:
	@echo > $@

#
# ${R_HOME}/src/extra/bzip2/Makefile


srcdir = .
top_srcdir = ../../..

top_builddir = ../../..
subdir = src/extra/bzip2
R_HOME = $(top_builddir)

include $(top_builddir)/Makeconf

BZIP2_CPPFLAGS = -I$(srcdir)
ALL_CPPFLAGS = $(BZIP2_CPPFLAGS) $(R_XTRA_CPPFLAGS) $(CPPFLAGS) $(DEFS)

SOURCES = blocksort.c bzlib.c bzcompress.c crctable.c decompress.c \
  huffman.c randtable.c
HEADERS = bzlib.h bzlib_private.h
DEPENDS = $(SOURCES:.c=.d)
OBJECTS = $(SOURCES:.c=.o)
#ALL_CFLAGS = $(ALL_CFLAGS_LO) -fvisibility=hidden

distdir = $(top_builddir)/$(PACKAGE)-$(VERSION)/$(subdir)
DISTFILES = LICENSE Makefile.in Makefile.win R_changes \
  $(SOURCES) $(HEADERS)
noinst_LIBRARIES = libbz2.a
#noinst_LIBRARIES =  libbz2.a
libbz2_a_SOURCES = $(SOURCES)
libbz2_a_OBJECTS = $(OBJECTS)


all: Makefile Makedeps R

Makefile: $(srcdir)/Makefile.in \
  $(top_builddir)/config.status \
  $(SOURCES)
	@cd $(top_builddir) && $(SHELL) ./config.status $(subdir)/$@

Makedeps: Makefile $(DEPENDS)
	@cat $(DEPENDS) >> Makefile
	@touch $@

R: Makefile
	@$(MAKE) Makedeps
	@$(MAKE) $(noinst_LIBRARIES)

libbz2.a: $(libbz2_a_OBJECTS)
	rm -f $@
	$(AR) cr $@ $(libbz2_a_OBJECTS)
	$(RANLIB) $@


mostlyclean: clean
clean:
	@-rm -f Makedeps *.d *.o *.lo *.a
distclean: clean
	@-rm -f Makefile
maintainer-clean: distclean

install install-strip uninstall TAGS info dvi check:

distdir: $(DISTFILES)
	@for f in $(DISTFILES); do \
	  test -f $(distdir)/$${f} \
	    || ln $(srcdir)/$${f} $(distdir)/$${f} 2>/dev/null \
	    || cp -p $(srcdir)/$${f} $(distdir)/$${f}; \
	done

## Automagically generated dependencies:

blocksort.o: blocksort.c bzlib_private.h ../../../src/include/config.h \
 bzlib.h
bzlib.o: bzlib.c bzlib_private.h ../../../src/include/config.h bzlib.h
bzcompress.o: bzcompress.c bzlib_private.h ../../../src/include/config.h \
 bzlib.h
crctable.o: crctable.c bzlib_private.h ../../../src/include/config.h \
 bzlib.h
decompress.o: decompress.c bzlib_private.h ../../../src/include/config.h \
 bzlib.h
huffman.o: huffman.c bzlib_private.h ../../../src/include/config.h \
 bzlib.h
randtable.o: randtable.c bzlib_private.h ../../../src/include/config.h \
 bzlib.h

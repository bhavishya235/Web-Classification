#
# ${R_HOME}/src/modules/X11/Makefile


srcdir = .
top_srcdir = ../../..

top_builddir = ../../..
subdir = src/modules/X11
R_HOME = $(top_builddir)

include $(top_builddir)/Makeconf

ALL_CPPFLAGS = $(R_XTRA_CPPFLAGS) $(X_CFLAGS) -pthread -I/usr/include/pango-1.0 -I/usr/include/cairo -I/usr/include/glib-2.0 -I/usr/lib/i386-linux-gnu/glib-2.0/include -I/usr/include/pixman-1 -I/usr/include/freetype2 -I/usr/include/libpng12   \
  -I$(top_srcdir)/src/library/grDevices/src/cairo $(CPPFLAGS) $(DEFS)
ALL_CFLAGS = $(ALL_CFLAGS_LO)
ALL_X_LIBS = $(X_PRE_LIBS) -lpangocairo-1.0 -lpango-1.0 -lgobject-2.0 -lglib-2.0 -lpng12 -lz -lcairo -lX11   $(X_LIBS) $(X_EXTRA_LIBS)

SOURCES = devX11.c rotated.c rbitmap.c 
DEPENDS = $(SOURCES:.c=.d)
OBJECTS = $(SOURCES:.c=.o)
DE_SOURCES = dataentry.c
DE_OBJECTS = $(DE_SOURCES:.c=.o)
HEADERS = bitmap.h devX11.h rotated.h rlogo_icon.h

distdir = $(top_builddir)/$(PACKAGE)-$(VERSION)/$(subdir)
DISTFILES = \
	Makefile.in $(SOURCES) $(HEADERS) $(DE_SOURCES)

R_X11_la = R_X11$(SHLIB_EXT)
## AIX needs LIBM
R_X11_la_LIBADD = $(ALL_X_LIBS) -ljpeg -lpng -lz $(LIBR) $(LIBM) # $(LIBINTL)

R_de_la = R_de$(SHLIB_EXT)
R_de_la_LIBADD = $(ALL_X_LIBS) $(LIBR) $(LIBM) # $(LIBINTL)

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
	@$(MAKE) $(R_X11_la) $(R_de_la)
	@$(MAKE) rhome="$(abs_top_builddir)" install

$(R_X11_la): $(OBJECTS)
	$(SHLIB_LINK) -o $@ $(OBJECTS) $(R_X11_la_LIBADD)

$(R_de_la): $(DE_OBJECTS)
	$(SHLIB_LINK) -o $@ $(DE_OBJECTS) $(R_de_la_LIBADD)

install: installdirs
	@$(SHELL) $(top_srcdir)/tools/copy-if-change $(R_X11_la) "$(DESTDIR)$(Rexecmodulesdir)/$(R_X11_la)"
	@$(SHELL) $(top_srcdir)/tools/copy-if-change $(R_de_la) "$(DESTDIR)$(Rexecmodulesdir)/$(R_de_la)"
installdirs:
	@$(MKINSTALLDIRS) "$(DESTDIR)$(Rexecmodulesdir)"
install-strip: installdirs
	@$(SHELL) $(top_srcdir)/tools/copy-if-change $(R_X11_la) "$(DESTDIR)$(Rexecmodulesdir)/$(R_X11_la)"
	@$(SHELL) $(top_srcdir)/tools/copy-if-change $(R_X11_la) "$(DESTDIR)$(Rexecmodulesdir)/$(R_de_la)"
	@if test -n "$(STRIP_LIBS)"; then \
	  $(STRIP_LIBS) "$(DESTDIR)$(Rexecmodulesdir)/$(R_X11_la)" ;\
	  $(STRIP_LIBS) "$(DESTDIR)$(Rexecmodulesdir)/$(R_de_la)" ;\
	fi
uninstall:
	@rm -f "$(Rexecmodulesdir)/R_X11$(SHLIB_EXT)"\
	  "$(Rexecmodulesdir)/R_de(SHLIB_EXT)"

mostlyclean: clean
clean:
	@-rm -rf .libs _libs
	@-rm -f Makedeps *.d *.o *.a *.lo *.la *$(SHLIB_EXT)
distclean: clean
	@-rm -f Makefile
maintainer-clean: distclean

TAGS info dvi check:

distdir: $(DISTFILES)
	@for f in $(DISTFILES); do \
	  test -f $(distdir)/$${f} \
	    || ln $(srcdir)/$${f} $(distdir)/$${f} 2>/dev/null \
	    || cp -p $(srcdir)/$${f} $(distdir)/$${f}; \
	done

## Automagically generated dependencies:

devX11.o: devX11.c ../../../src/include/config.h \
 ../../../src/include/Defn.h ../../../src/include/R_ext/Complex.h \
 ../../../src/include/Rinternals.h ../../../src/include/R_ext/Arith.h \
 ../../../src/include/R_ext/libextern.h \
 ../../../src/include/R_ext/Boolean.h ../../../src/include/R_ext/Error.h \
 ../../../src/include/R_ext/Memory.h ../../../src/include/R_ext/Utils.h \
 ../../../src/include/R_ext/Print.h ../../../src/include/Rinlinedfuns.h \
 ../../../src/include/Errormsg.h \
 ../../../src/include/R_ext/GraphicsEngine.h \
 ../../../src/include/R_ext/GraphicsDevice.h \
 ../../../src/include/Fileio.h rotated.h \
 ../../../src/include/R_ext/eventloop.h devX11.h \
 /usr/include/pango-1.0/pango/pango.h \
 /usr/include/pango-1.0/pango/pango-attributes.h \
 /usr/include/pango-1.0/pango/pango-font.h \
 /usr/include/pango-1.0/pango/pango-coverage.h \
 /usr/include/glib-2.0/glib.h /usr/include/glib-2.0/glib/galloca.h \
 /usr/include/glib-2.0/glib/gtypes.h \
 /usr/lib/i386-linux-gnu/glib-2.0/include/glibconfig.h \
 /usr/include/glib-2.0/glib/gmacros.h \
 /usr/include/glib-2.0/glib/gversionmacros.h \
 /usr/include/glib-2.0/glib/garray.h \
 /usr/include/glib-2.0/glib/gasyncqueue.h \
 /usr/include/glib-2.0/glib/gthread.h \
 /usr/include/glib-2.0/glib/gatomic.h /usr/include/glib-2.0/glib/gerror.h \
 /usr/include/glib-2.0/glib/gquark.h \
 /usr/include/glib-2.0/glib/gbacktrace.h \
 /usr/include/glib-2.0/glib/gbase64.h \
 /usr/include/glib-2.0/glib/gbitlock.h \
 /usr/include/glib-2.0/glib/gbookmarkfile.h \
 /usr/include/glib-2.0/glib/gbytes.h \
 /usr/include/glib-2.0/glib/gcharset.h \
 /usr/include/glib-2.0/glib/gchecksum.h \
 /usr/include/glib-2.0/glib/gconvert.h \
 /usr/include/glib-2.0/glib/gdataset.h /usr/include/glib-2.0/glib/gdate.h \
 /usr/include/glib-2.0/glib/gdatetime.h \
 /usr/include/glib-2.0/glib/gtimezone.h /usr/include/glib-2.0/glib/gdir.h \
 /usr/include/glib-2.0/glib/genviron.h \
 /usr/include/glib-2.0/glib/gfileutils.h \
 /usr/include/glib-2.0/glib/ggettext.h /usr/include/glib-2.0/glib/ghash.h \
 /usr/include/glib-2.0/glib/glist.h /usr/include/glib-2.0/glib/gmem.h \
 /usr/include/glib-2.0/glib/ghmac.h \
 /usr/include/glib-2.0/glib/gchecksum.h \
 /usr/include/glib-2.0/glib/ghook.h \
 /usr/include/glib-2.0/glib/ghostutils.h \
 /usr/include/glib-2.0/glib/giochannel.h \
 /usr/include/glib-2.0/glib/gmain.h /usr/include/glib-2.0/glib/gpoll.h \
 /usr/include/glib-2.0/glib/gslist.h /usr/include/glib-2.0/glib/gstring.h \
 /usr/include/glib-2.0/glib/gunicode.h \
 /usr/include/glib-2.0/glib/gutils.h \
 /usr/include/glib-2.0/glib/gkeyfile.h \
 /usr/include/glib-2.0/glib/gmappedfile.h \
 /usr/include/glib-2.0/glib/gmarkup.h \
 /usr/include/glib-2.0/glib/gmessages.h \
 /usr/include/glib-2.0/glib/gnode.h /usr/include/glib-2.0/glib/goption.h \
 /usr/include/glib-2.0/glib/gpattern.h \
 /usr/include/glib-2.0/glib/gprimes.h /usr/include/glib-2.0/glib/gqsort.h \
 /usr/include/glib-2.0/glib/gqueue.h /usr/include/glib-2.0/glib/grand.h \
 /usr/include/glib-2.0/glib/gregex.h \
 /usr/include/glib-2.0/glib/gscanner.h \
 /usr/include/glib-2.0/glib/gsequence.h \
 /usr/include/glib-2.0/glib/gshell.h /usr/include/glib-2.0/glib/gslice.h \
 /usr/include/glib-2.0/glib/gspawn.h \
 /usr/include/glib-2.0/glib/gstrfuncs.h \
 /usr/include/glib-2.0/glib/gstringchunk.h \
 /usr/include/glib-2.0/glib/gtestutils.h \
 /usr/include/glib-2.0/glib/gthreadpool.h \
 /usr/include/glib-2.0/glib/gtimer.h \
 /usr/include/glib-2.0/glib/gtrashstack.h \
 /usr/include/glib-2.0/glib/gtree.h \
 /usr/include/glib-2.0/glib/gurifuncs.h \
 /usr/include/glib-2.0/glib/gvarianttype.h \
 /usr/include/glib-2.0/glib/gvariant.h \
 /usr/include/glib-2.0/glib/gversion.h \
 /usr/include/glib-2.0/glib/deprecated/gallocator.h \
 /usr/include/glib-2.0/glib/deprecated/gcache.h \
 /usr/include/glib-2.0/glib/deprecated/gcompletion.h \
 /usr/include/glib-2.0/glib/deprecated/gmain.h \
 /usr/include/glib-2.0/glib/deprecated/grel.h \
 /usr/include/glib-2.0/glib/deprecated/gthread.h \
 /usr/include/pango-1.0/pango/pango-types.h \
 /usr/include/glib-2.0/glib-object.h \
 /usr/include/glib-2.0/gobject/gbinding.h \
 /usr/include/glib-2.0/gobject/gobject.h \
 /usr/include/glib-2.0/gobject/gtype.h \
 /usr/include/glib-2.0/gobject/gvalue.h \
 /usr/include/glib-2.0/gobject/gparam.h \
 /usr/include/glib-2.0/gobject/gclosure.h \
 /usr/include/glib-2.0/gobject/gsignal.h \
 /usr/include/glib-2.0/gobject/gmarshal.h \
 /usr/include/glib-2.0/gobject/gboxed.h \
 /usr/include/glib-2.0/gobject/glib-types.h \
 /usr/include/glib-2.0/gobject/genums.h \
 /usr/include/glib-2.0/gobject/gparamspecs.h \
 /usr/include/glib-2.0/gobject/gsourceclosure.h \
 /usr/include/glib-2.0/gobject/gtypemodule.h \
 /usr/include/glib-2.0/gobject/gtypeplugin.h \
 /usr/include/glib-2.0/gobject/gvaluearray.h \
 /usr/include/glib-2.0/gobject/gvaluetypes.h \
 /usr/include/pango-1.0/pango/pango-gravity.h \
 /usr/include/pango-1.0/pango/pango-matrix.h \
 /usr/include/pango-1.0/pango/pango-script.h \
 /usr/include/pango-1.0/pango/pango-language.h \
 /usr/include/pango-1.0/pango/pango-bidi-type.h \
 /usr/include/pango-1.0/pango/pango-break.h \
 /usr/include/pango-1.0/pango/pango-item.h \
 /usr/include/pango-1.0/pango/pango-context.h \
 /usr/include/pango-1.0/pango/pango-fontmap.h \
 /usr/include/pango-1.0/pango/pango-fontset.h \
 /usr/include/pango-1.0/pango/pango-engine.h \
 /usr/include/pango-1.0/pango/pango-glyph.h \
 /usr/include/pango-1.0/pango/pango-enum-types.h \
 /usr/include/pango-1.0/pango/pango-features.h \
 /usr/include/pango-1.0/pango/pango-glyph-item.h \
 /usr/include/pango-1.0/pango/pango-layout.h \
 /usr/include/pango-1.0/pango/pango-tabs.h \
 /usr/include/pango-1.0/pango/pango-renderer.h \
 /usr/include/pango-1.0/pango/pango-utils.h \
 /usr/include/pango-1.0/pango/pangocairo.h /usr/include/cairo/cairo.h \
 /usr/include/cairo/cairo-version.h /usr/include/cairo/cairo-features.h \
 /usr/include/cairo/cairo-deprecated.h /usr/include/cairo/cairo-xlib.h \
 /usr/include/cairo/cairo.h rlogo_icon.h \
 ../../../src/include/Rmodules/RX11.h ../../../src/include/Rconnections.h \
 ../../../src/include/R_ext/Connections.h \
 ../../../src/library/grDevices/src/cairo/cairoFns.c bitmap.h \
 ../../../src/include/R_ext/Rdynload.h
rotated.o: rotated.c ../../../src/include/config.h rotated.h
rbitmap.o: rbitmap.c ../../../src/include/config.h bitmap.h \
 /usr/include/libpng12/png.h /usr/include/libpng12/pngconf.h \
 ../../../src/include/R_ext/Error.h

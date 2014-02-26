##
## ${R_HOME}/src/appl/Makefile


srcdir = .
top_srcdir = ../..

top_builddir = ../..
subdir = src/appl

include $(top_builddir)/Makeconf

SOURCES_C = \
	integrate.c interv.c maxcol.c optim.c pretty.c uncmin.c
SOURCES_F = \
	dchdc.f dpbfa.f dpbsl.f dpoco.f dpodi.f dpofa.f dposl.f dqrdc.f \
	dqrdc2.f dqrls.f dqrsl.f dqrutl.f dsvdc.f dtrco.f dtrsl.f
DEPENDS = $(SOURCES_C:.c=.d)
OBJECTS_BLAS =  blas.o  cmplxblas.o
OBJECTS = $(SOURCES_C:.c=.o) $(SOURCES_F:.f=.o) \
  # $(OBJECTS_BLAS)
HEADERS = lbfgsb.c

distdir = $(top_builddir)/$(PACKAGE)-$(VERSION)/$(subdir)
DISTFILES = \
	Makefile.in Makefile.win \
	$(HEADERS) \
	$(SOURCES_C) \
	$(SOURCES_F)

#ALL_CFLAGS = $(ALL_CFLAGS_LO)
#ALL_FFLAGS = $(ALL_FFLAGS_LO)

Rblas_la = libRblas$(DYLIB_EXT)
Rblas_la_LIBADD = # -flat_namespace -undefined suppress

all: Makefile Makedeps R

Makefile: $(srcdir)/Makefile.in $(top_builddir)/config.status $(SOURCES_C)
	@cd $(top_builddir) && $(SHELL) ./config.status $(subdir)/$@

Makedeps: Makefile $(DEPENDS)
	@cat $(DEPENDS) >> Makefile
	@touch $@

R: Makefile
	@$(MAKE) Makedeps
	@$(MAKE) libappl.a

libappl.a: $(OBJECTS)
	rm -f $@
	$(AR) cr $@ $(OBJECTS)
	$(RANLIB) $@

blas.o: $(top_srcdir)/src/extra/blas/blas.f
	$(F77) $(ALL_FFLAGS) -c $< -o $@
cmplxblas.o: $(top_srcdir)/src/extra/blas/cmplxblas.f
	$(F77) $(ALL_FFLAGS) -c $< -o $@

mostlyclean: clean
clean:
	@-rm -rf .libs _libs
	@-rm -f Makedeps *.d *.o *.a *.lo *.la stamp-lo $(Rblas_la)

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

integrate.o: integrate.c ../../src/include/config.h \
 ../../src/include/Rmath.h ../../src/include/R_ext/Boolean.h \
 ../../src/include/R_ext/Applic.h ../../src/include/R_ext/RS.h \
 ../../src/include/Rconfig.h ../../src/include/R_ext/Error.h \
 ../../src/include/R_ext/BLAS.h ../../src/include/R_ext/Complex.h
interv.o: interv.c ../../src/include/R_ext/Applic.h \
 ../../src/include/R_ext/Boolean.h ../../src/include/R_ext/RS.h \
 ../../src/include/Rconfig.h ../../src/include/R_ext/Error.h \
 ../../src/include/R_ext/BLAS.h ../../src/include/R_ext/Complex.h \
 ../../src/include/R_ext/Utils.h
maxcol.o: maxcol.c ../../src/include/config.h \
 ../../src/include/R_ext/Arith.h ../../src/include/R_ext/libextern.h \
 ../../src/include/Rmath.h ../../src/include/R_ext/Boolean.h \
 ../../src/include/R_ext/Random.h ../../src/include/R_ext/Applic.h \
 ../../src/include/R_ext/RS.h ../../src/include/Rconfig.h \
 ../../src/include/R_ext/Error.h ../../src/include/R_ext/BLAS.h \
 ../../src/include/R_ext/Complex.h
optim.o: optim.c ../../src/include/config.h ../../src/include/Defn.h \
 ../../src/include/R_ext/Complex.h ../../src/include/Rinternals.h \
 ../../src/include/R_ext/Arith.h ../../src/include/R_ext/libextern.h \
 ../../src/include/R_ext/Boolean.h ../../src/include/R_ext/Error.h \
 ../../src/include/R_ext/Memory.h ../../src/include/R_ext/Utils.h \
 ../../src/include/R_ext/Print.h ../../src/include/Rinlinedfuns.h \
 ../../src/include/Errormsg.h ../../src/include/R_ext/Random.h \
 ../../src/include/R_ext/Applic.h ../../src/include/R_ext/RS.h \
 ../../src/include/Rconfig.h ../../src/include/R_ext/BLAS.h lbfgsb.c \
 ../../src/include/R_ext/Linpack.h
pretty.o: pretty.c ../../src/include/config.h \
 ../../src/include/R_ext/Applic.h ../../src/include/R_ext/Boolean.h \
 ../../src/include/R_ext/RS.h ../../src/include/Rconfig.h \
 ../../src/include/R_ext/Error.h ../../src/include/R_ext/BLAS.h \
 ../../src/include/R_ext/Complex.h ../../src/include/Rmath.h
uncmin.o: uncmin.c ../../src/include/R_ext/Applic.h \
 ../../src/include/R_ext/Boolean.h ../../src/include/R_ext/RS.h \
 ../../src/include/Rconfig.h ../../src/include/R_ext/Error.h \
 ../../src/include/R_ext/BLAS.h ../../src/include/R_ext/Complex.h \
 ../../src/include/R_ext/Print.h ../../src/include/R_ext/PrtUtil.h \
 ../../src/include/Rinternals.h ../../src/include/R_ext/Arith.h \
 ../../src/include/R_ext/libextern.h ../../src/include/R_ext/Memory.h \
 ../../src/include/R_ext/Utils.h ../../src/include/R_ext/Linpack.h \
 ../../src/include/Rmath.h

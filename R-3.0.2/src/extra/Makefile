#
# ${R_HOME}/src/extra/Makefile


srcdir = .
top_srcdir = ../..

top_builddir = ../..
subdir = src/extra

include $(top_builddir)/Makeconf

distdir = $(top_builddir)/$(PACKAGE)-$(VERSION)/$(subdir)
DISTFILES = Makefile.in

SUBDIRS = blas bzip2 pcre tre xdr xz zlib  intl
SUBDIRS_WITH_NO_BUILD = graphapp trio tzone win_iconv # intl

RBLAS = make.blas
LIBBZ2 = make.bzip2
LIBPCRE = make.pcre
LIBZ = make.zlib
LIBXDR =# make.xdr
LIBINTL =# make.intl
LIBTRE = make.tre
LIBXZ = make.xz
EXTRAS = $(RBLAS) $(LIBBZ2) $(LIBPCRE) $(LIBZ) $(LIBXDR) $(LIBINTL) $(LIBTRE) $(LIBXZ)

all: Makefile R

Makefile: $(srcdir)/Makefile.in $(top_builddir)/config.status
	@cd $(top_builddir) && $(SHELL) ./config.status $(subdir)/$@

## can be done in parallel
R: Makefile $(EXTRAS)

make.blas: Makefile
	@(cd blas; $(MAKE))

make.bzip2: Makefile
	@(cd bzip2; $(MAKE))

make.pcre: Makefile
	@(cd pcre; $(MAKE))

make.zlib: Makefile
	@(cd zlib; $(MAKE))

make.tre: Makefile
	@(cd tre; $(MAKE))

make.xdr: Makefile
	@(cd xdr; $(MAKE))

make.xz: Makefile
	@(cd xz; $(MAKE))

make.intl: Makefile
	@(cd intl; $(MAKE))


mostlyclean: clean
clean:
	@for d in $(SUBDIRS); do \
	  (cd $${d} && $(MAKE) $@) || exit 1; \
	done
distclean: clean
	@for d in $(SUBDIRS); do \
	  (cd $${d} && $(MAKE) $@) || exit 1; \
	done
	-@rm -f Makefile
maintainer-clean: distclean

install install-strip:
	@for d in $(SUBDIRS); do \
	  (cd $${d} && $(MAKE) $@) || exit 1; \
	done

uninstall TAGS info dvi check:

distdir: $(DISTFILES)
	@for f in $(DISTFILES); do \
	  test -f $(distdir)/$${f} \
	    || ln $(srcdir)/$${f} $(distdir)/$${f} 2>/dev/null \
	    || cp -p $(srcdir)/$${f} $(distdir)/$${f}; \
	done
	@for d in $(SUBDIRS); do \
	  test -d $(distdir)/$${d} \
	    || mkdir $(distdir)/$${d} \
	    || exit 1; \
	  chmod 755 $(distdir)/$${d}; \
	  (cd $${d} && $(MAKE) distdir) \
	    || exit 1; \
	done
	@for d in $(SUBDIRS_WITH_NO_BUILD); do \
	  ((cd $(srcdir); $(TAR) -c -f - $(DISTDIR_TAR_EXCLUDE) $${d}) \
	      | (cd $(distdir); $(TAR) -x -f -)) \
	    || exit 1; \
	done

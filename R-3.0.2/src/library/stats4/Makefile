#
# ${R_HOME}/src/library/stats4/Makefile


srcdir = .
top_srcdir = ../../..

top_builddir = ../../..
subdir = src/library/stats4

include $(top_builddir)/Makeconf

distdir = $(top_builddir)/$(PACKAGE)-$(VERSION)/$(subdir)
DISTFILES = DESCRIPTION.in Makefile.in Makefile.win NAMESPACE

pkg = stats4

R_EXE = $(top_builddir)/bin/R --vanilla --slave
DEFPKGS = "methods,graphics,stats"

all: Makefile DESCRIPTION
	@$(ECHO) "building package '$(pkg)'"
	@$(MKINSTALLDIRS) $(top_builddir)/library/$(pkg)
	@$(MAKE) mkR2 mkdesc
#	@$(MAKE) mklazy
	@$(MAKE) mklazycomp

mkR: mkR2

include $(top_srcdir)/share/make/basepkg.mk
include $(top_srcdir)/share/make/lazycomp.mk

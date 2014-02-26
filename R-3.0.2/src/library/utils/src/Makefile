#
# ${R_HOME}/src/library/utils/src/Makefile

pkg = utils


srcdir = .
top_srcdir = ../../../..

top_builddir = ../../../..
subdir = src/library/$(pkg)/src
R_HOME = $(top_builddir)
R_SHARE_DIR = $(R_HOME)/share
R_INCLUDE_DIR = $(R_HOME)/include

SOURCES_C = init.c io.c size.c sock.c stubs.c utils.c

DEPENDS = $(SOURCES_C:.c=.d)
OBJECTS = $(SOURCES_C:.c=.o)

PKG_CFLAGS = $(C_VISIBILITY)

SHLIB = $(pkg).so

XZ_CPPFLAGS =   -I$(top_srcdir)/src/extra/xz/api
# need Defn.h etc, and config.h
PKG_CPPFLAGS =-I../../../include -I$(top_srcdir)/src/include -DHAVE_CONFIG_H \
  -I$(top_srcdir)/src/main $(XZ_CPPFLAGS)


all: Makefile Makedeps
	@$(MAKE) Makedeps
	@$(MAKE) shlib

Makefile: $(srcdir)/Makefile.in $(top_builddir)/config.status
	@cd $(top_builddir) && $(SHELL) ./config.status $(subdir)/$@

Makedeps: Makefile $(DEPENDS)
	@cat $(DEPENDS) >> Makefile
	@touch $@

shlib: $(SHLIB)
	@$(MKINSTALLDIRS) \
	  $(top_builddir)/library/$(pkg)/libs
	@$(SHELL) $(top_srcdir)/tools/copy-if-change $(SHLIB) \
	  $(top_builddir)/library/$(pkg)/libs/$(SHLIB)

include $(R_HOME)/etc${R_ARCH}/Makeconf
include $(top_srcdir)/share/make/shlib.mk
LTO = 

mostlyclean: clean
clean:
	@-rm -rf .libs _libs
	@-rm -f Makedeps *.d *.o *$(SHLIB_EXT)
distclean: clean
	@-rm -f Makefile
maintainer-clean: distclean

## Automagically generated dependencies:
init.o: init.c ../../../include/config.h ../../../../include/R.h \
 ../../../../include/Rconfig.h ../../../../include/R_ext/Arith.h \
 ../../../../include/R_ext/libextern.h \
 ../../../../include/R_ext/Boolean.h ../../../../include/R_ext/Complex.h \
 ../../../../include/R_ext/Constants.h ../../../../include/R_ext/Error.h \
 ../../../../include/R_ext/Memory.h ../../../../include/R_ext/Print.h \
 ../../../../include/R_ext/Random.h ../../../../include/R_ext/Utils.h \
 ../../../../include/R_ext/RS.h ../../../../include/Rinternals.h utils.h \
 ../../../../include/R_ext/Rdynload.h
io.o: io.c ../../../include/config.h ../../../include/Defn.h \
 ../../../../include/R_ext/Complex.h ../../../../include/Rinternals.h \
 ../../../../include/R_ext/Arith.h ../../../../include/R_ext/libextern.h \
 ../../../../include/R_ext/Boolean.h ../../../../include/R_ext/Error.h \
 ../../../../include/R_ext/Memory.h ../../../../include/R_ext/Utils.h \
 ../../../../include/R_ext/Print.h ../../../include/Rinlinedfuns.h \
 ../../../include/Errormsg.h ../../../include/Fileio.h \
 ../../../include/Rconnections.h ../../../../include/R_ext/Connections.h \
 ../../../include/Print.h ../../../include/Defn.h \
 ../../../../include/R_ext/PrtUtil.h ../../../include/rlocale.h \
 ../../../../src/main/RBufferUtils.h
size.o: size.c ../../../include/config.h ../../../include/Defn.h \
 ../../../../include/R_ext/Complex.h ../../../../include/Rinternals.h \
 ../../../../include/R_ext/Arith.h ../../../../include/R_ext/libextern.h \
 ../../../../include/R_ext/Boolean.h ../../../../include/R_ext/Error.h \
 ../../../../include/R_ext/Memory.h ../../../../include/R_ext/Utils.h \
 ../../../../include/R_ext/Print.h ../../../include/Rinlinedfuns.h \
 ../../../include/Errormsg.h
sock.o: sock.c ../../../../include/R.h ../../../../include/Rconfig.h \
 ../../../../include/R_ext/Arith.h ../../../../include/R_ext/libextern.h \
 ../../../../include/R_ext/Boolean.h ../../../../include/R_ext/Complex.h \
 ../../../../include/R_ext/Constants.h ../../../../include/R_ext/Error.h \
 ../../../../include/R_ext/Memory.h ../../../../include/R_ext/Print.h \
 ../../../../include/R_ext/Random.h ../../../../include/R_ext/Utils.h \
 ../../../../include/R_ext/RS.h ../../../../include/Rinternals.h
stubs.o: stubs.c ../../../include/config.h ../../../include/Defn.h \
 ../../../../include/R_ext/Complex.h ../../../../include/Rinternals.h \
 ../../../../include/R_ext/Arith.h ../../../../include/R_ext/libextern.h \
 ../../../../include/R_ext/Boolean.h ../../../../include/R_ext/Error.h \
 ../../../../include/R_ext/Memory.h ../../../../include/R_ext/Utils.h \
 ../../../../include/R_ext/Print.h ../../../include/Rinlinedfuns.h \
 ../../../include/Errormsg.h ../../../include/Internal.h \
 ../../../../include/Rinterface.h ../../../../include/R_ext/RStartup.h \
 ../../../include/Rdynpriv.h ../../../../include/R_ext/Rdynload.h \
 ../../../include/Rmodules/RX11.h ../../../include/Rconnections.h \
 ../../../../include/R_ext/Connections.h
utils.o: utils.c ../../../include/config.h ../../../../include/R.h \
 ../../../../include/Rconfig.h ../../../../include/R_ext/Arith.h \
 ../../../../include/R_ext/libextern.h \
 ../../../../include/R_ext/Boolean.h ../../../../include/R_ext/Complex.h \
 ../../../../include/R_ext/Constants.h ../../../../include/R_ext/Error.h \
 ../../../../include/R_ext/Memory.h ../../../../include/R_ext/Print.h \
 ../../../../include/R_ext/Random.h ../../../../include/R_ext/Utils.h \
 ../../../../include/R_ext/RS.h ../../../../include/Rinternals.h utils.h \
 ../../../../src/extra/xz/api/lzma.h \
 ../../../../src/extra/xz/api/lzma/version.h \
 ../../../../src/extra/xz/api/lzma/base.h \
 ../../../../src/extra/xz/api/lzma/vli.h \
 ../../../../src/extra/xz/api/lzma/check.h \
 ../../../../src/extra/xz/api/lzma/filter.h \
 ../../../../src/extra/xz/api/lzma/bcj.h \
 ../../../../src/extra/xz/api/lzma/delta.h \
 ../../../../src/extra/xz/api/lzma/lzma.h \
 ../../../../src/extra/xz/api/lzma/container.h \
 ../../../../src/extra/xz/api/lzma/stream_flags.h \
 ../../../../src/extra/xz/api/lzma/block.h \
 ../../../../src/extra/xz/api/lzma/index.h \
 ../../../../src/extra/xz/api/lzma/index_hash.h \
 ../../../../src/extra/xz/api/lzma/hardware.h

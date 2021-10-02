PKG_CFLAGS = $(C_VISIBILITY) @CFLAGS@ @OPENMP_CFLAGS@
PKG_CPPFLAGS = @CPPFLAGS@ @GHMM_CPPFLAGS@
PKG_LIBS = @OPENMP_CFLAGS@ @LDFLAGS@ @GHMM_LDFLAGS@ @LIBS@ @GHMM_LIBS@

if BUNDLED_GHMM
SUBDIRS = ghmm-0.9-rc3 .
endif

lib_LTLIBRARIES    = tsshmm.la
tsshmm_la_CFLAGS   = $(PKG_CFLAGS)
tsshmm_la_CPPFLAGS = $(PKG_CPPFLAGS)
tsshmm_la_LDFLAGS  = -module $(PKG_LIBS)
tsshmm_la_SOURCES  = R_init_tsshmm.c R_wrap_tsshmm.c models.c parameters.c \
train.c tss.c viterbi.c

ACLOCAL_AMFLAGS = -I tools

# Hook that runs after the default "all" rule.
if BUNDLED_GHMM
all-local : tsshmm.so libghmm.so
else
all-local : tsshmm.so
endif

# One of the limitations with POSIX-compliant `make` is not being able to
# specify multiple outputs from a single rule.  Therefore, even though libtool
# links and produces multiple files, the only explicit rule we can depend on is
# the .la output file.
tsshmm.so : tsshmm.la
	cp -av .libs/tsshmm.so.0.0.0 $@
	chmod -x $@

if BUNDLED_GHMM
libghmm.so : @GHMM_ROOT@/libghmm.la
	cp -av @GHMM_ROOT@/.libs/libghmm.so.1.0.0 $@
	chmod -x $@
endif

clean-local :
	rm -f tsshmm.so libghmm.so

distclean-local : clean-local

.PHONY : all-local clean-local distclean-local

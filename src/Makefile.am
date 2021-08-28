PKG_CFLAGS = $(C_VISIBILITY) @CFLAGS@ @OPENMP_CFLAGS@
PKG_CPPFLAGS = @CPPFLAGS@ @GHMM_CPPFLAGS@
PKG_LIBS = @OPENMP_CFLAGS@ @LDFLAGS@ @GHMM_LDFLAGS@ @LIBS@ @GHMM_LIBS@
SUBDIRS = ghmm-0.9-rc3 .

lib_LTLIBRARIES    = tsshmm.la
tsshmm_la_CFLAGS   = $(PKG_CFLAGS)
tsshmm_la_CPPFLAGS = $(PKG_CPPFLAGS)
tsshmm_la_LDFLAGS  = -module
tsshmm_la_LIBS     = $(PKG_LIBS)
tsshmm_la_SOURCES  = R_init_tsshmm.c R_wrap_tsshmm.c tss.c viterbi.c

# Hook that runs after the default "all" rule.
all-local :
	cp -av .libs/tsshmm.so.0.0.0 tsshmm.so
	chmod -x tsshmm.so
if BUNDLED_GHMM
	cp -av @GHMM_ROOT@/.libs/libghmm.so.1.0.0 libghmm.so
	chmod -x libghmm.so
endif

clean-local :
	rm -f tsshmm.so libghmm.so

distclean-local : clean-local

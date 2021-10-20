if BUNDLED_GHMM
SUBDIRS = ghmm-0.9-rc3 .
endif

lib_LTLIBRARIES = tsshmm.la
tsshmm_la_CFLAGS = @CFLAGS@ @OPENMP_CFLAGS@ $(C_VISIBILITY)
tsshmm_la_CPPFLAGS = @CPPFLAGS@ @GHMM_CPPFLAGS@
if BUNDLED_GHMM
tsshmm_la_LIBADD = @LIBS@ @GHMM_ROOT@/libghmm.la @LIBXML2_LIBS@
else
tsshmm_la_LIBADD = @LIBS@ @GHMM_LIBS@
endif
tsshmm_la_LDFLAGS = @LDFLAGS@ @OPENMP_CFLAGS@ @GHMM_LDFLAGS@ -module
tsshmm_la_SOURCES = R_init_tsshmm.c R_wrap_tsshmm.c models.c simulate.c	\
train.c tss.c viterbi.c

ACLOCAL_AMFLAGS = -I tools

# Hook that runs after the default "all" rule.
all-local : tsshmm.so

# One of the limitations with POSIX-compliant `make` is not being able to
# specify multiple outputs from a single rule.  Therefore, even though libtool
# links and produces multiple files, the only explicit rule we can depend on is
# the .la output file.
tsshmm.so : tsshmm.la
	cp -av .libs/tsshmm.so.0.0.0 $@
	chmod -x $@

clean-local :
	rm -f tsshmm.so

distclean-local : clean-local

.PHONY : all-local clean-local distclean-local

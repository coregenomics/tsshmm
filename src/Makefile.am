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

# We cannot use autoreconf because we want to store all autotools output files
# inside tools/
#
# libtoolize generates: tools/{ltmain.sh,libtool.m4,ltoptions.m4,ltsugar.m4}
#                       tools/{ltversion.m4,lt~obsolete.m4}
libtoolize --quiet
# aclocal generates: tools/aclocal.m4
aclocal --output=tools/aclocal.m4
# autoconf generates: configure
autoconf
# automake generates: tools/{compile,install-sh,missing}
automake --add-missing

export R_LIBS_USER=$PWD/library
export R_PROFILE_USER=$PWD/.Rprofile
mkdir -vp "$R_LIBS_USER"
cat <<EOF > "$R_PROFILE_USER"
local({
    .libPaths("$R_LIBS_USER")
    options(repos = BiocManager::repositories(),
            Ncpus = parallel::detectCores())
})
EOF
# Don't track the GHMM source files in version control, because running
# configure modifies many tracked files which creates version control noise.
# Additionally, the GHMM tarball guarantees the integrity of the source, and
# makes us disciplined in applying patches making for easier upstream
# contributions.
pushd src/ || exit
tar -xf ghmm-0.9-rc3.tar.gz
popd || exit
# We cannot use autoreconf because we don't want to dirty our parent directory
# and fail R CMD check, and instead want to store all autotools output files
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

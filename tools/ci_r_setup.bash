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

export R_LIBS_USER=$PWD/library
export R_PROFILE_USER=$PWD/.Rprofile
mkdir -vp "$R_LIBS_USER"
cat <<EOF > "$R_PROFILE_USER"
local({
    .libPaths("$R_LIBS_USER")
    options(repos = BiocManager::repositories())
})
EOF

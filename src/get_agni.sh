#!/bin/bash
# Update and compile AGNI

set -e

# User arguments
if [ -z "$1" ]; then
    testsuite="fast" # default
else
    testsuite="$1" # requested test suite
fi

# Update from GitHub (if we are in an updateable state)
if git pull 2>/dev/null; then
    echo "Updated from GitHub"
else
    echo "Not updating AGNI from GitHub. When installing AGNI through the PROTEUS installer this is expected, since it checks out a specific commit. If you want to update AGNI, you should install it manually or request an update through PROTEUS."
fi
root=$(dirname $(realpath $0))
root=$(realpath "$root/..")

# SOCRATES?
if [ -n "$RAD_DIR" ]; then
    echo "Found SOCRATES path: yes"
else
    echo "Found SOCRATES path: no"
    echo "    You need to install SOCRATES and set the RAD_DIR environment variable"
    echo "    This requires modifying your shell rc file to export RAD_DIR"
    echo "    Check the docs: https://www.h-nicholls.space/AGNI/"
    exit 1
fi

# FastChem?
if [ -n "$FC_DIR" ]; then
    echo "Found FastChem path: yes"
else
    echo "Found FastChem path: no"
    echo "    It is highly recommended that you install FastChem and set the FC_DIR environment variable"
    echo "    Check the docs: https://www.h-nicholls.space/AGNI/"
fi

# Install
echo "Installing AGNI..."
rm -f Manifest.toml
julia "$root/deps/build.jl"
julia -e 'using Pkg; Pkg.activate("."); Pkg.instantiate()'

# Run tests (suite set by CLI argument)
echo "Running tests..."
dir=$(pwd)
cd $root/test/
julia runtests.jl "$testsuite"
cd $dir

echo "Done!"
exit 0

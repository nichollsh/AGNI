#!/bin/bash
# Update and compile AGNI

set -e

# User arguments
if [ -z "$1" ]; then
    testsuite="all" # default
else
    testsuite="$1" # requested test suite
fi

# Update from GitHub
echo "Updating from GitHub..."
git pull
root=$(dirname $(realpath $0))
root=$(realpath "$root/..")

# SOCRATES?
if [ -n "$RAD_DIR" ]; then
    echo "Found SOCRATES path: yes"
else
    echo "Found SOCRATES path: no"
    echo "You need to install SOCRATES AND set the RAD_DIR environment variable"
    echo "Check the docs: https://www.h-nicholls.space/AGNI/"
    exit 1
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

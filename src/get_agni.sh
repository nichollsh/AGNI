#!/bin/bash
# Update and compile AGNI

set -e

# Update from GitHub
echo "Updating from GitHub..."
git pull
root=$(dirname $(realpath $0))
root=$(realpath "$root/..")

# SOCRATES?
if [ -n "$RAD_DIR" ]; then
    echo "SOCRATES is already installed"
else
    echo "You need to install SOCRATES and set RAD_DIR first. Check the docs."
    exit 1
fi

# Install
echo "Installing AGNI..."
rm -f Manifest.toml
julia -e 'using Pkg; Pkg.activate("."); Pkg.instantiate(); Pkg.build()'

# Run tests
echo "Running tests..."
julia $root/test/runtests.jl

echo "Done!"
exit 0

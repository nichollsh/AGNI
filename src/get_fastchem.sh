#!/bin/bash
# Download and compile fastchem

tag="v3.1.3"

set -e

root=$(dirname $(realpath $0))
root=$(realpath "$root/..")

# Download via HTTPS only
fcpath="$root/fastchem"
rm -rf "$fcpath"

git clone --depth 1 --branch "$tag" https://github.com/NewStrangeWorlds/FastChem.git "$fcpath"

# Compile FastChem
cd "$fcpath"
mkdir build
cd build

cmake ".."
make

cd $root
export FC_DIR=$fcpath

echo "FastChem has been installed"
echo "It is recommended that you add the following line to your shell rc file"
echo "export FC_DIR='$fcpath'"
exit 0

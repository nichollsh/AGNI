#!/bin/bash
# Download and compile fastchem

# Download via HTTPS only
fcpath="fastchem"
rm -rf "$fcpath"

git clone https://github.com/NewStrangeWorlds/FastChem.git "$fcpath"

# Compile FastChem
cd "$fcpath"
mkdir build
cd build 

cmake ".."
make

cd ..
export FC_DIR="$(realpath $fcpath)"

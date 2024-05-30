#!/bin/bash

# Download
rm -rf socrates
git clone git@github.com:nichollsh/SOCRATES.git socrates

# Environment
export LD_LIBRARY_PATH=""

# Install 
cd socrates
./configure
./build_code

# Environment
source ./set_rad_env
export LD_LIBRARY_PATH=""

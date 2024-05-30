#!/bin/bash


rm -rf socrates

# Download
git clone git@github.com:nichollsh/SOCRATES.git socrates

# Environment
export LD_LIBRARY_PATH=""

# Install 
cd socrates
./configure
./build_code

# Environment
source ./set_rad_env

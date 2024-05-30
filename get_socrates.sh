#!/bin/bash

# Remove old repo
rm -rf socrates

# Download
git clone git@github.com:nichollsh/SOCRATES.git socrates

# Install 
cd socrates
./configure
./build_code

# Environment
source ./set_rad_env

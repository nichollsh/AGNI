#!/bin/bash

# Check ssh access to GitHub
ssh -T git@github.com
if [ $? -eq 1 ]; then 
    use_ssh=true
else 
    use_ssh=false
fi

# Download
socpath="socrates"
rm -rf "$socpath"
if [ "$use_ssh" = true]; then
    git clone git@github.com:nichollsh/SOCRATES.git "$socpath"
else
    git clone https://github.com/nichollsh/SOCRATES.git "$socpath"
fi 

# Environment
export LD_LIBRARY_PATH=""

# Install 
cd "$socpath"
./configure
./build_code

# Environment
source ./set_rad_env
export LD_LIBRARY_PATH=""

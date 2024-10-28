#!/bin/bash
# Download and compile socrates

# Already setup?
if [ -n "$RAD_DIR" ]; then
    echo "WARNING: You already have SOCRATES installed"
    echo "         RAD_DIR=$RAD_DIR"
    exit 1
fi

# Check SSH access to GitHub
ssh -T git@github.com
if [ $? -eq 1 ]; then
    use_ssh=true
else
    use_ssh=false
fi

# Download (using SSH if possible)
socpath="$(realpath .)/socrates"
rm -rf "$socpath"
if [ "$use_ssh" = true ]; then
    git clone git@github.com:nichollsh/SOCRATES.git "$socpath"
else
    git clone https://github.com/nichollsh/SOCRATES.git "$socpath"
fi

# Compile SOCRATES
cd "$socpath"
./configure
./build_code

# Environment
export RAD_DIR=$socpath
cd ..

# Inform user
echo "SOCRATES has been installed"
echo "It is recommended that you add the following line to your shell rc file"
echo "export RAD_DIR='$socpath'"
exit 0

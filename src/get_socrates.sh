#!/bin/bash
# Download and compile socrates

auto_yes=false
while getopts ":y" opt; do
    case "$opt" in
        y) auto_yes=true ;;
        *)
            echo "Usage: $0 [-y]" >&2
            exit 1
            ;;
    esac
done
shift $((OPTIND - 1))

if [ $# -ne 0 ]; then
    echo "Usage: $0 [-y]" >&2
    exit 1
fi

confirm_yes() {
    local prompt="$1"
    if [ "$auto_yes" = true ]; then
        echo "$prompt [y/N] y (-y)"
        return 0
    fi

    read -r -p "$prompt [y/N] " choice
    case "$choice" in
        [Yy]|[Yy][Ee][Ss]) return 0 ;;
        *) return 1 ;;
    esac
}

warn_duplicate_exports() {
    local var_name="$1"
    local shell_rc="$2"
    local export_count

    if [ ! -f "$shell_rc" ]; then
        return
    fi

    export_count=$(grep -Ec "^[[:space:]]*export[[:space:]]+${var_name}=" "$shell_rc")
    if [ "$export_count" -gt 1 ]; then
        echo "WARNING: Found ${export_count} export entries for ${var_name} in '$shell_rc'."
        echo "         Consider keeping only one to avoid ambiguity."
    fi
}

# Do we have NetCDF?
if ! [ -x "$(command -v nc-config)" ]; then
  echo 'ERROR: NetCDF is not installed.' >&2
  exit 1
fi
if ! [ -x "$(command -v nf-config)" ]; then
  echo 'ERROR: NetCDF-Fortran library is not installed.' >&2
  exit 1
fi

# Do we have gfortran?
if ! [ -x "$(command -v gfortran)" ]; then
  echo 'ERROR: gfortran compiler is not installed.' >&2
  exit 1
fi

# Root paths
root=$(dirname "$(realpath "$0")")
root=$(realpath "$root/..")
default_socpath="$root/socrates"

# Installation path
if [ -n "$RAD_DIR" ]; then
    socpath="$RAD_DIR"
    echo "RAD_DIR is already set: '$socpath'"
    if ! confirm_yes "Install SOCRATES there and overwrite existing contents if needed?"; then
        echo "Exiting without changes."
        exit 0
    fi
else
    socpath="$default_socpath"
fi

if [ -z "$socpath" ] || [ "$socpath" = "/" ]; then
    echo "ERROR: Invalid installation path: '$socpath'" >&2
    exit 1
fi

# Check SSH access to GitHub
ssh -T git@github.com
if [ $? -eq 1 ]; then
    use_ssh=true
else
    use_ssh=false
fi

# Disable SSH (uncomment to allow SSH clone of SOCRATES)
use_ssh=false

# Download
mkdir -p "$(dirname "$socpath")"
rm -rf "$socpath"
if [ "$use_ssh" = true ]; then
    git clone git@github.com:FormingWorlds/SOCRATES.git "$socpath"
else
    git clone https://github.com/FormingWorlds/SOCRATES.git "$socpath"
fi

# Compile SOCRATES
cd "$socpath"
./configure
./build_code

# Environment
export RAD_DIR="$socpath"
cd "$root"

# Check radlib exists
radlib="$socpath/bin/radlib.a"
if [ -f "$radlib" ]; then
    echo "SOCRATES has been installed"
    echo ""
else
    echo "Could not find compiled SOCRATES binaries - failed to compile"
    exit 1
fi

detect_shell_rc() {
    case "$(basename "${SHELL:-}")" in
        bash) echo "$HOME/.bashrc" ;;
        zsh)  echo "$HOME/.zshrc" ;;
        ksh)  echo "$HOME/.kshrc" ;;
        *)    echo "$HOME/.profile" ;;
    esac
}

# Inform user
echo "You must now run the following command:"
echo "    export RAD_DIR='$socpath'"

rcfile="$(detect_shell_rc)"
if confirm_yes "Add RAD_DIR to '$rcfile'?"; then
    export_line="export RAD_DIR='$socpath'"
    touch "$rcfile"
    {
        echo ""
        echo "# SOCRATES installation"
        echo "$export_line"
    } >> "$rcfile"
    echo "Added RAD_DIR to '$rcfile'."
    echo "Please restart your terminal or run 'source $rcfile' to apply the changes."
    echo " "
else
    echo "Skipped editing shell rc file."
fi
warn_duplicate_exports "RAD_DIR" "$rcfile"

exit 0

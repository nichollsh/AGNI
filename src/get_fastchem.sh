#!/bin/bash
# Download and compile fastchem

tag="main" # Tag or branch

set -e

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

detect_shell_rc() {
    case "$(basename "${SHELL:-}")" in
        bash) echo "$HOME/.bashrc" ;;
        zsh)  echo "$HOME/.zshrc" ;;
        ksh)  echo "$HOME/.kshrc" ;;
        *)    echo "$HOME/.profile" ;;
    esac
}

root=$(dirname "$(realpath "$0")")
root=$(realpath "$root/..")

# Installation path
default_fcpath="$root/fastchem"
if [ -n "$FC_DIR" ]; then
    fcpath="$FC_DIR"
    echo "FC_DIR is already set: '$fcpath'"
    if ! confirm_yes "Install FastChem there and overwrite existing contents if needed?"; then
        echo "Exiting without changes."
        exit 0
    fi
else
    fcpath="$default_fcpath"
fi

if [ -z "$fcpath" ] || [ "$fcpath" = "/" ]; then
    echo "ERROR: Invalid installation path: '$fcpath'" >&2
    exit 1
fi

# Download via HTTPS only
mkdir -p "$(dirname "$fcpath")"
rm -rf "$fcpath"

git clone --depth 1 --branch "$tag" https://github.com/FormingWorlds/FastChem.git "$fcpath"

# Compile FastChem
cd "$fcpath"
mkdir build
cd build

cmake ".."
make

cd "$root"
export FC_DIR="$fcpath"

# Check fastchem executable exists
fcexec="$fcpath/fastchem"
if [ -f "$fcexec" ]; then
    echo "FastChem has been installed"
    echo ""
else
    echo "Could not find FastChem executable - failed to compile"
    exit 1
fi

rcfile="$(detect_shell_rc)"
if confirm_yes "Add FC_DIR to '$rcfile'?"; then
    export_line="export FC_DIR='$fcpath'"
    touch "$rcfile"
    {
        echo ""
        echo "# FastChem installation"
        echo "$export_line"
    } >> "$rcfile"
    echo "Added FC_DIR to '$rcfile'."
    echo "Please restart your terminal or run 'source $rcfile' to apply the changes."
    echo " "
else
    echo "Skipped editing shell rc file."
fi
warn_duplicate_exports "FC_DIR" "$rcfile"

exit 0
